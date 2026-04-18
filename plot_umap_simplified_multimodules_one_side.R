#!/usr/bin/env Rscript
# Single-population UMAP / LocusZoom / Coloc-Z-score plot.
# Per module the figure row is: [LocusZoom] -> Beta UMAP -> [Coloc Z-scores].

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(data.table)
  library(optparse)
  library(patchwork)
  library(ggrepel)
  library(ragg)
  library(scattermore)
})

# --- CLI Options ---
option_list <- list(
  make_option("--colors", type = "character", default = NULL, help = "Path to palette TSV (only required if --show_ref is used)"),
  
  # Required
  make_option("--umap", type = "character", help = "Path to UMAP TSV"),
  make_option("--master", type = "character", help = "Path to Master Annotations TSV"),
  make_option("--modules", type = "character", help = "Comma-separated Module IDs"),
  
  # Optional per-module side panels
  make_option("--summary_table", type = "character", default = NULL, help = "Path to Summary Table (for disease titles)"),
  make_option("--z_files", type = "character", default = NULL, help = "Comma-separated paths to Z-Score CSVs (use NA for missing)"),
  make_option("--lz_files", type = "character", default = NULL, help = "Comma-separated paths to LocusZoom CSVs (use NA for missing). Expected columns: CHR,CELL,GENE,POS,P"),
  make_option("--name", type = "character", default = "Pop", help = "Display name used in panel titles"),
  
  # General Settings
  make_option("--join_col", type = "character", default = "celltype_2", help = "Column name in UMAP for celltypes [default: %default]"),
  make_option("--anno_join_col", type = "character", default = NULL, help = "Column name in Master Annotations (if different from join_col)"),
  make_option("--gene_col", type = "character", default = "eGene_symbol", help = "Column name for gene symbol [default: %default]"),
  make_option("--gtf", type = "character", default = NULL, help = "Optional GENCODE GTF (e.g. gencode.v49.annotation.gtf). When supplied, the LocusZoom gene track shows every protein-coding gene in the window. A cached TSV is created beside the GTF on first use."),
  make_option("--annotations", type = "character", default = NULL, help = "Optional comma-separated annotation file(s) to show in a zoom panel around the lead SNP. Each file may be simple 4-col (chrom, feature_type, start, end) or full 9-col GFF. Feature types become separate tracks."),
  make_option("--zoom_window", type = "numeric", default = 5000, help = "Half-width in bp for the zoom annotation panel around the lead SNP [default: %default]"),
  
  make_option("--out", type = "character", default = "umap_plot", help = "Output filename prefix"),
  make_option("--pt_size", type = "numeric", default = 0.25, help = "Point size"),
  make_option("--max_cells", type = "numeric", default = 250000, help = "Max cells to plot"),
  
  # Flipped Toggle Flags
  make_option("--png", action="store_true", default=FALSE, help="Save as PNG instead of PDF (PDF is default)"),
  make_option("--no_raster", action="store_true", default=FALSE, help="Disable rasterizing points (Raster is default)"),
  make_option("--show_ref", action="store_true", default=FALSE, help="Show the reference UMAP plot (Hidden by default)"),
  make_option("--no_labels", action="store_true", default=FALSE, help="Hide active cell type labels (Shown by default)")
)

opt <- parse_args(OptionParser(option_list = option_list))
anno_col <- if (!is.null(opt$anno_join_col)) opt$anno_join_col else opt$join_col

save_pdf   <- !opt$png
use_raster <- !opt$no_raster
show_ref   <- opt$show_ref
show_labels <- !opt$no_labels

# --- Validation & Parsing ---
required_args <- c("umap", "master", "modules")
missing_args <- setdiff(required_args, names(opt))
if (length(missing_args) > 0) stop(paste("Missing required arguments:", paste(missing_args, collapse = ", ")))
if (show_ref && is.null(opt$colors)) stop("You must provide the --colors file if you enable --show_ref.")

mods    <- trimws(unlist(strsplit(opt$modules, ",")))
has_z   <- !is.null(opt$z_files)
has_lz  <- !is.null(opt$lz_files)

if (has_z) {
  z_files <- trimws(unlist(strsplit(opt$z_files, ",")))
  if (length(z_files) != length(mods)) stop("--z_files must have the same number of items as --modules")
}
if (has_lz) {
  lz_files <- trimws(unlist(strsplit(opt$lz_files, ",")))
  if (length(lz_files) != length(mods)) stop("--lz_files must have the same number of items as --modules")
}

# --- Helper Functions ---

# plot_ref(): reference UMAP coloured by celltype palette; labels celltype
# centroids with ggrepel. Used when --show_ref is enabled.
#
# Example:
#   plot_ref(merged, "GH Reference", "celltype_2", 0.25, TRUE)
#   # -> ggplot of the full UMAP in palette colours, labelled.
plot_ref <- function(df, title_text, join_col, pt_size, use_raster) {
  centroids <- df %>% group_by(.data[[join_col]]) %>%
    summarise(UMAP_1 = median(UMAP_1, na.rm = TRUE),
              UMAP_2 = median(UMAP_2, na.rm = TRUE),
              .groups = "drop")
  p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2))
  if (use_raster) p <- p + geom_scattermore(aes(color = hex_color), pointsize = pt_size + 1.5, pixels = c(2048, 2048)) else p <- p + geom_point(aes(color = hex_color), size = pt_size, stroke = 0)
  p + scale_color_identity() +
    geom_text_repel(data = centroids, aes(label = .data[[join_col]]),
                    size = 3.5, fontface = "bold", bg.color = "white",
                    bg.r = 0.1, min.segment.length = 0) +
    coord_equal() + labs(title = title_text, x = NULL, y = NULL) +
    theme_minimal() +
    theme(panel.grid = element_blank(), axis.text = element_blank(),
          plot.title = element_text(face = "bold"))
}

# plot_beta(): UMAP coloured by a per-cell-type beta_val, with a symmetric
# blue-yellow-red diverging palette anchored at 0.
#
# Example:
#   plot_beta(plot_df, "GH | RBFA\n[M_18448]", 0.25, "celltype_2",
#             show_legend=TRUE, use_raster=TRUE, show_labels=TRUE)
plot_beta <- function(df, title, pt_size, join_col, show_legend = TRUE,
                      use_raster = FALSE, show_labels = FALSE) {
  b_max <- max(abs(df$beta_val), na.rm = TRUE)
  if (!is.finite(b_max) || b_max <= 0) b_max <- 1e-9
  colors_beta <- c("#00008B", "#1E90FF", "#90CAF9", "#E8E8E8", "#FFEB3B", "#FF9800", "#E53935")
  
  p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = beta_val))
  if (use_raster) p <- p + geom_scattermore(pointsize = pt_size + 1.5, pixels = c(2048, 2048)) else p <- p + geom_point(size = pt_size, stroke = 0)
  
  p <- p + scale_color_gradientn(colors = colors_beta,
                                 limits = c(-b_max, b_max), na.value = "grey50",
                                 name = "Beta",
                                 values = c(0.0, 0.20, 0.45, 0.5, 0.55, 0.80, 1.0)) +
    coord_equal() + labs(title = title, x = NULL, y = NULL) + theme_minimal() +
    theme(panel.grid = element_blank(), axis.text = element_blank(),
          plot.title  = element_text(size = 8.5, face = "bold"),
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 10))
  
  if (!show_legend) p <- p + theme(legend.position = "none")
  if (show_labels) {
    centroids <- df %>% group_by(.data[[join_col]]) %>%
      summarise(UMAP_1 = median(UMAP_1, na.rm = TRUE),
                UMAP_2 = median(UMAP_2, na.rm = TRUE),
                mod_id = if ("mod_id_col" %in% names(.)) { v <- na.omit(mod_id_col); if (length(v) > 0) v[1] else NA_character_ } else { NA_character_ },
                .groups = "drop") %>%
      filter(!is.na(mod_id) & mod_id != "") %>%
      mutate(label_text = .data[[join_col]])
    p <- p + geom_text_repel(data = centroids,
                             aes(x = UMAP_1, y = UMAP_2, label = label_text),
                             inherit.aes = FALSE, size = 4.5,
                             fontface = "bold", bg.color = "white", bg.r = 0.15,
                             color = "black", min.segment.length = 0)
  }
  return(p)
}

# plot_zscore(): QTL vs disease Z-score scatter for a single module.
#
# Example:
#   plot_zscore("icd10_j15_zscore_merge_notnull.csv", "Coloc Z-Scores\n[M_18448]")
#   # -> ggplot showing points with diagonals and dashed zero lines.
plot_zscore <- function(z_file, title = "Coloc Z-Scores") {
  if (is.na(z_file) || z_file == "NA" || !file.exists(z_file)) {
    if (!is.na(z_file) && z_file != "NA") cat(sprintf("Warning: Z-score file not found at %s\n", z_file))
    return(plot_spacer())
  }
  df <- fread(z_file)
  qtl_col <- names(df)[grepl("z_qtl", names(df), ignore.case = TRUE)][1]
  if (is.na(qtl_col)) qtl_col <- names(df)[grepl("qtl", names(df), ignore.case = TRUE)][1]
  dis_col <- names(df)[grepl("z_icd10|z_dis", names(df), ignore.case = TRUE)][1]
  if (is.na(dis_col)) dis_col <- names(df)[!names(df) %in% c("snp", qtl_col)][1]
  if (is.na(qtl_col) || is.na(dis_col)) {
    cat(sprintf("Warning: Could not identify Z-score columns in %s\n", z_file))
    return(plot_spacer())
  }
  max_val   <- max(abs(c(df[[qtl_col]], df[[dis_col]])), na.rm = TRUE)
  limit_val <- max_val * 1.1
  ggplot(df, aes(x = .data[[qtl_col]], y = .data[[dis_col]])) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = 0, color = "black", linetype = "dashed", alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, color = "gray80", linetype = "dotted") +
    geom_abline(slope = -1, intercept = 0, color = "gray80", linetype = "dotted") +
    geom_point(fill = "#1E90FF", color = "white", shape = 21,
               size = 2.5, alpha = 0.8, stroke = 0.3) +
    coord_fixed(xlim = c(-limit_val, limit_val), ylim = c(-limit_val, limit_val)) +
    theme_minimal() +
    labs(title = title, x = "QTL Z-Score", y = "Disease Z-Score") +
    theme(plot.title  = element_text(size = 10, face = "bold"),
          axis.title  = element_text(size = 10, face = "bold"),
          axis.text   = element_text(size = 9),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "gray90", fill = NA, linewidth = 1))
}

# load_gene_table(): return a data.table of protein-coding gene bodies ready
# to use for the LocusZoom gene track. Accepts either:
#   (a) a pre-built cache TSV with columns chrom, start, end, strand,
#       gene_name, gene_id -- loaded directly; or
#   (b) a raw GENCODE GTF -- in which case a cache TSV is built next to the
#       GTF (re-built only when the GTF is newer than the cache).
#
# Cache header is detected by reading the first non-empty line of the file.
#
# Example:
#   load_gene_table("gencode.v49.annotation.gtf")
#   # -> builds "gencode.v49.annotation.gtf.coding_genes.tsv" (~20k rows)
#   load_gene_table("gencode.v49.annotation.gtf.coding_genes.tsv")
#   # -> loads the cached TSV directly, no awk step.
load_gene_table <- function(gtf_path) {
  if (!file.exists(gtf_path)) stop(sprintf("GTF / gene cache not found: %s", gtf_path))
  expected_header <- c("chrom", "start", "end", "strand", "gene_name", "gene_id")
  first_line <- readLines(gtf_path, n = 1, warn = FALSE)
  if (length(first_line) > 0) {
    tokens <- strsplit(first_line, "\t", fixed = TRUE)[[1]]
    if (identical(tokens, expected_header)) {
      cat(sprintf("Loading pre-built gene cache %s ...\n", gtf_path))
      return(fread(gtf_path))
    }
  }
  cache_path <- paste0(gtf_path, ".coding_genes.tsv")
  if (file.exists(cache_path) &&
      file.info(cache_path)$mtime >= file.info(gtf_path)$mtime) {
    cat(sprintf("Loading existing gene cache %s ...\n", cache_path))
    return(fread(cache_path))
  }
  cat(sprintf("Building protein-coding gene cache from %s (one-off, ~1 min) ...\n", gtf_path))
  tmp <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp), add = TRUE)
  cmd <- sprintf("awk -F'\\t' '$0 !~ /^#/ && $3 == \"gene\"' %s > %s",
                 shQuote(gtf_path), shQuote(tmp))
  if (system(cmd) != 0) stop("Failed to pre-filter GTF with awk.")
  raw <- fread(tmp, sep = "\t", header = FALSE,
               col.names = c("chrom", "source", "feature", "start", "end",
                             "score", "strand", "frame", "attrs"))
  raw[, gene_type := sub('.*gene_type "([^"]+)".*', "\\1", attrs)]
  raw <- raw[gene_type == "protein_coding"]
  raw[, gene_name := sub('.*gene_name "([^"]+)".*', "\\1", attrs)]
  raw[, gene_id   := sub('.*gene_id "([^"]+)".*',   "\\1", attrs)]
  out <- raw[, .(chrom, start, end, strand, gene_name, gene_id)]
  fwrite(out, cache_path, sep = "\t")
  cat(sprintf("Cached %d protein-coding genes to %s\n", nrow(out), cache_path))
  out
}

# get_genes_in_region(): protein-coding genes overlapping [win_start,win_end]
# on chromosome chr.
#
# Example:
#   get_genes_in_region(gtf_tbl, "chr18", 79800000, 80200000)
#   # -> data.table with RBFA, ATP9B, NFATC1 ... rows.
get_genes_in_region <- function(gtf_tbl, chr, win_start, win_end) {
  if (is.null(gtf_tbl) || nrow(gtf_tbl) == 0) return(NULL)
  chr_no <- sub("^chr", "", chr, ignore.case = TRUE)
  targets <- unique(c(chr, chr_no, paste0("chr", chr_no)))
  hit <- gtf_tbl[chrom %in% targets & start <= win_end & end >= win_start]
  if (nrow(hit) == 0) return(NULL)
  hit
}

# assign_gene_lanes(): greedy lane packing so overlapping gene arrows stack
# in separate lanes. Adds a `lane` column (1 = bottom).
#
# Example:
#   assign_gene_lanes(data.frame(start=c(1,5,20), end=c(10,15,30)))
#   # -> lanes = c(1, 2, 1)
assign_gene_lanes <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  df <- df[order(df$start), , drop = FALSE]
  lanes <- integer(nrow(df))
  lane_ends <- numeric(0)
  for (i in seq_len(nrow(df))) {
    placed <- FALSE
    for (k in seq_along(lane_ends)) {
      if (df$start[i] > lane_ends[k]) {
        lanes[i] <- k; lane_ends[k] <- df$end[i]; placed <- TRUE; break
      }
    }
    if (!placed) { lanes[i] <- length(lane_ends) + 1L; lane_ends <- c(lane_ends, df$end[i]) }
  }
  df$lane <- lanes
  df
}

# load_annotation_file(): read a regulatory annotation file into a tidy
# data.table with columns chrom, feature_type, start, end. Accepts either
# the 4-column simple format (chrom, feature_type, start, end) or a full
# 9-column GFF/GTF. Chromosome names are normalised by stripping any
# leading "chr".
#
# Example:
#   load_annotation_file("homo_sapiens...promoter.gff")
#   # -> data.table chrom=18, feature_type="promoter_flanking_region", ...
load_annotation_file <- function(path) {
  if (!file.exists(path)) stop(sprintf("Annotation file not found: %s", path))
  con <- file(path, "r"); on.exit(close(con), add = TRUE)
  first_data_line <- NULL
  while (length(line <- readLines(con, n = 1)) > 0) {
    if (!grepl("^#", line) && nchar(line) > 0) { first_data_line <- line; break }
  }
  if (is.null(first_data_line)) return(NULL)
  n_cols <- length(strsplit(first_data_line, "\t", fixed = TRUE)[[1]])
  
  if (n_cols >= 9) {
    raw <- fread(path, sep = "\t", header = FALSE, skip = "#",
                 col.names = c("chrom", "source", "feature_type", "start", "end",
                               "score", "strand", "frame", "attrs"),
                 select = 1:9, fill = TRUE)
    out <- raw[, .(chrom = as.character(chrom),
                   feature_type = as.character(feature_type),
                   start = as.integer(start),
                   end   = as.integer(end))]
  } else if (n_cols == 4) {
    raw <- fread(path, sep = "\t", header = FALSE,
                 col.names = c("chrom", "feature_type", "start", "end"))
    out <- raw[, .(chrom = as.character(chrom),
                   feature_type = as.character(feature_type),
                   start = as.integer(start),
                   end   = as.integer(end))]
  } else {
    stop(sprintf("Unsupported annotation format (%d columns) in %s", n_cols, path))
  }
  out[, chrom := sub("^chr", "", chrom, ignore.case = TRUE)]
  out
}

# filter_annotations_window(): subset an annotation data.table to features
# overlapping [win_start, win_end] on chromosome chr.
#
# Example:
#   filter_annotations_window(anno_tbl, "chr18", 80000273, 80010273)
#   # -> data.table containing the promoter_flanking_region at
#   #    80009802-80011199.
filter_annotations_window <- function(tbl, chr, win_start, win_end) {
  if (is.null(tbl) || nrow(tbl) == 0) return(tbl[0])
  chr_no <- sub("^chr", "", chr, ignore.case = TRUE)
  tbl[chrom == chr_no & start <= win_end & end >= win_start]
}

# plot_zoom_annotations(): narrow zoom panel around the lead SNP with one
# lane per feature_type, filled rectangles for overlapping features and a
# dashed purple vertical line at the lead SNP. Empty panels still render
# with a small "No annotations in window" note.
#
# Example:
#   plot_zoom_annotations(filter_annotations_window(anno_tbl, "18",
#                                                    80000273, 80010273),
#                         "18", 80005273, 5000,
#                         c("promoter", "promoter_flanking_region"))
plot_zoom_annotations <- function(anno_win, chr_label, lead_pos, win_half,
                                  track_order) {
  win_s <- lead_pos - win_half
  win_e <- lead_pos + win_half
  n_tracks <- max(length(track_order), 1L)
  pal <- c("#1565C0", "#2E7D32", "#EF6C00", "#6A1B9A", "#00838F",
           "#AD1457", "#558B2F", "#4E342E", "#283593", "#F57F17")
  track_colors <- setNames(rep(pal, length.out = n_tracks), track_order)
  
  p <- ggplot() +
    scale_y_continuous(breaks = seq_along(track_order), labels = track_order,
                       limits = c(0.3, n_tracks + 0.7)) +
    scale_fill_manual(values = track_colors, guide = "none", limits = track_order) +
    scale_x_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
    coord_cartesian(xlim = c(win_s, win_e), expand = FALSE) +
    labs(title = sprintf("Zoom +/- %s bp around lead SNP (chr%s:%s)",
                         format(win_half, big.mark = ",", scientific = FALSE),
                         chr_label,
                         format(lead_pos, big.mark = ",", scientific = FALSE)),
         x = paste0("Chromosome ", chr_label, " position (bp)"), y = NULL) +
    theme_minimal() +
    theme(plot.title         = element_text(size = 9, face = "bold"),
          axis.title.x       = element_text(size = 9, face = "bold"),
          axis.text          = element_text(size = 8),
          panel.grid.minor   = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.border       = element_rect(color = "gray80", fill = NA, linewidth = 0.8))
  
  if (!is.null(anno_win) && nrow(anno_win) > 0) {
    anno_win <- as.data.frame(anno_win)
    anno_win$lane    <- match(anno_win$feature_type, track_order)
    anno_win$x_start <- pmax(anno_win$start, win_s)
    anno_win$x_end   <- pmin(anno_win$end,   win_e)
    anno_win <- anno_win[!is.na(anno_win$lane), , drop = FALSE]
    p <- p + geom_rect(data = anno_win,
                       aes(xmin = x_start, xmax = x_end,
                           ymin = lane - 0.35, ymax = lane + 0.35,
                           fill = feature_type),
                       color = NA, alpha = 0.85)
  } else {
    p <- p + annotate("text", x = lead_pos, y = (n_tracks + 1) / 2,
                      label = "No annotations in window",
                      size = 3.2, color = "gray55", fontface = "italic")
  }
  p + geom_vline(xintercept = lead_pos, linetype = "dashed",
                 color = "#8E24AA", alpha = 0.8, linewidth = 0.5)
}

# plot_locuszoom(): composite LocusZoom panel.
#   - main: -log10(P) vs position (Mb) with lead SNP highlighted
#   - gene track: all protein-coding genes in the window (focal eGene red)
#   - optional zoom: regulatory annotations within +/- zoom_window bp
#
# Example:
#   plot_locuszoom("j15_locuszoom_ENSG00000101546.csv",
#                  locus_info = list(chrom="chr18", gene_start=79850000,
#                                    gene_end=79960000, gene_strand="+",
#                                    gene_symbol="RBFA", lead_pos=80005273),
#                  title = "LocusZoom | GH | RBFA\n[M_18448]",
#                  genes_df = get_genes_in_region(gtf_tbl,"chr18",
#                                                 79800000, 80100000))
plot_locuszoom <- function(lz_file, locus_info = NULL, title = "LocusZoom",
                           genes_df = NULL,
                           annotations_tbl = NULL,
                           annotation_tracks = NULL,
                           zoom_window = 5000) {
  if (is.na(lz_file) || lz_file == "NA" || !file.exists(lz_file)) {
    if (!is.na(lz_file) && lz_file != "NA") cat(sprintf("Warning: LocusZoom file not found at %s\n", lz_file))
    return(plot_spacer())
  }
  df <- fread(lz_file)
  needed <- c("CHR", "POS", "P")
  if (!all(needed %in% names(df))) {
    cat(sprintf("Warning: %s is missing required columns (CHR, POS, P)\n", lz_file))
    return(plot_spacer())
  }
  df <- df[!is.na(P) & is.finite(as.numeric(P))]
  if (nrow(df) == 0) return(plot_spacer())
  df[, logp := -log10(pmax(as.numeric(P), 1e-300))]
  
  chr_label <- sub("^chr", "", as.character(df$CHR[1]), ignore.case = TRUE)
  
  lead_pos <- if (!is.null(locus_info) && !is.null(locus_info$lead_pos)) locus_info$lead_pos else NA_real_
  if (is.na(lead_pos) || !is.finite(as.numeric(lead_pos))) {
    top_idx <- which.max(df$logp); lead_pos <- df$POS[top_idx]
  } else {
    lead_pos <- as.numeric(lead_pos)
    if (!any(df$POS == lead_pos)) {
      lead_pos <- df$POS[which.min(abs(df$POS - lead_pos))]
    }
  }
  df[, is_lead := POS == lead_pos]
  lead_row <- df[is_lead == TRUE][1]
  lead_label <- sprintf("chr%s:%s\n-log10(P) = %.1f",
                        chr_label,
                        format(lead_row$POS, big.mark = ",", scientific = FALSE),
                        lead_row$logp)
  
  has_focal <- !is.null(locus_info) &&
               !is.null(locus_info$gene_start) && !is.na(locus_info$gene_start) &&
               !is.null(locus_info$gene_end)   && !is.na(locus_info$gene_end)
  has_many  <- !is.null(genes_df) && nrow(genes_df) > 0
  
  x_min_bp <- min(df$POS); x_max_bp <- max(df$POS)
  if (has_focal) {
    x_min_bp <- min(x_min_bp, as.numeric(locus_info$gene_start))
    x_max_bp <- max(x_max_bp, as.numeric(locus_info$gene_end))
  }
  if (has_many) {
    x_min_bp <- min(x_min_bp, min(genes_df$start))
    x_max_bp <- max(x_max_bp, max(genes_df$end))
  }
  pad <- max((x_max_bp - x_min_bp) * 0.02, 1)
  x_lims_mb <- c((x_min_bp - pad) / 1e6, (x_max_bp + pad) / 1e6)
  y_max <- max(df$logp, na.rm = TRUE) * 1.15
  
  p_main <- ggplot(df[is_lead == FALSE], aes(x = POS / 1e6, y = logp)) +
    geom_point(shape = 21, fill = "#4F8EDC", color = "white",
               size = 1.6, alpha = 0.85, stroke = 0.25) +
    geom_point(data = df[is_lead == TRUE], aes(x = POS / 1e6, y = logp),
               shape = 23, fill = "#8E24AA", color = "black",
               size = 3.5, stroke = 0.6) +
    geom_text_repel(data = df[is_lead == TRUE],
                    aes(x = POS / 1e6, y = logp, label = lead_label),
                    size = 2.9, fontface = "bold", color = "#4A148C",
                    nudge_y = y_max * 0.08, segment.color = "#8E24AA",
                    min.segment.length = 0, box.padding = 0.3) +
    coord_cartesian(xlim = x_lims_mb, ylim = c(0, y_max)) +
    labs(title = title, x = NULL, y = expression(-log[10](italic(P)))) +
    theme_minimal() +
    theme(plot.title = element_text(size = 10, face = "bold"),
          axis.title.y = element_text(size = 10, face = "bold"),
          axis.text = element_text(size = 9),
          axis.text.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.8),
          plot.margin = margin(5, 5, 0, 5))
  
  build_zoom_panel <- function() {
    if (is.null(annotations_tbl) || is.null(annotation_tracks) ||
        length(annotation_tracks) == 0) return(NULL)
    anno_win <- filter_annotations_window(annotations_tbl, chr_label,
                                          lead_pos - zoom_window,
                                          lead_pos + zoom_window)
    plot_zoom_annotations(anno_win, chr_label, lead_pos, zoom_window,
                          annotation_tracks)
  }
  
  if (!has_focal && !has_many) {
    p_main <- p_main +
      labs(x = paste0("Chromosome ", chr_label, " position (Mb)")) +
      theme(axis.text.x = element_text(size = 9),
            axis.title.x = element_text(size = 10, face = "bold"))
    p_zoom <- build_zoom_panel()
    if (is.null(p_zoom)) return(p_main)
    zoom_h <- max(1.2, 0.6 * length(annotation_tracks) + 0.6)
    return(p_main / p_zoom + plot_layout(heights = c(5, zoom_h)))
  }
  
  focal_sym <- if (!is.null(locus_info) &&
                   !is.null(locus_info$gene_symbol) &&
                   !is.na(locus_info$gene_symbol)) as.character(locus_info$gene_symbol) else NA_character_
  
  if (has_many) {
    track_df <- data.frame(
      start = as.numeric(genes_df$start),
      end   = as.numeric(genes_df$end),
      strand = as.character(genes_df$strand),
      gene_name = as.character(genes_df$gene_name),
      stringsAsFactors = FALSE)
  } else {
    track_df <- data.frame(
      start = as.numeric(locus_info$gene_start),
      end   = as.numeric(locus_info$gene_end),
      strand = if (!is.null(locus_info$gene_strand)) as.character(locus_info$gene_strand) else "+",
      gene_name = if (!is.na(focal_sym)) focal_sym else "gene",
      stringsAsFactors = FALSE)
  }
  track_df$is_focal <- !is.na(focal_sym) & track_df$gene_name == focal_sym
  track_df <- assign_gene_lanes(track_df)
  n_lanes <- max(track_df$lane, 1L)
  track_df$seg_x    <- ifelse(track_df$strand == "-", track_df$end,   track_df$start) / 1e6
  track_df$seg_xend <- ifelse(track_df$strand == "-", track_df$start, track_df$end)   / 1e6
  track_df$mid_mb   <- (track_df$start + track_df$end) / 2 / 1e6
  track_df$col      <- ifelse(track_df$is_focal, "#C62828", "#1A237E")
  
  p_gene <- ggplot(track_df) +
    geom_segment(aes(x = seg_x, xend = seg_xend, y = lane, yend = lane,
                     color = I(col)),
                 linewidth = 1.1,
                 arrow = arrow(ends = "last", type = "closed",
                               length = unit(0.07, "inches"))) +
    geom_text(aes(x = mid_mb, y = lane + 0.35, label = gene_name,
                  color = I(col),
                  fontface = ifelse(is_focal, "bold.italic", "italic")),
              size = 2.9) +
    scale_y_continuous(limits = c(0.3, n_lanes + 0.9), breaks = NULL) +
    coord_cartesian(xlim = x_lims_mb, clip = "off") +
    labs(x = paste0("Chromosome ", chr_label, " position (Mb)"), y = NULL) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x = element_text(size = 10, face = "bold"),
          axis.text.x = element_text(size = 9),
          panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.8),
          plot.margin = margin(0, 5, 5, 5))
  
  # Row heights scale with the number of stacked gene lanes and annotation
  # tracks so labels stay readable when the locus is gene-dense.
  gene_weight <- max(1, min(n_lanes, 5))
  p_zoom <- build_zoom_panel()
  if (is.null(p_zoom)) {
    return(p_main / p_gene + plot_layout(heights = c(5, gene_weight)))
  }
  zoom_h <- max(1.2, 0.6 * length(annotation_tracks) + 0.6)
  p_main / p_gene / p_zoom + plot_layout(heights = c(5, gene_weight, zoom_h))
}

# get_module_locus_info(): locus coordinates + lead SNP for a module, pulled
# from the master annotation table. Returns a list with chrom, gene_start,
# gene_end, gene_strand, gene_symbol, lead_pos (any field may be NA).
#
# Example:
#   get_module_locus_info(meta, "M_18448", gene_col = "eGene_symbol")
#   # -> list(chrom="chr18", gene_start=79850000, gene_end=79960000,
#   #         gene_strand="+", gene_symbol="RBFA", lead_pos=80005273)
get_module_locus_info <- function(df, mod_id, gene_col = "eGene_symbol") {
  sub_df <- df %>% filter(module == mod_id)
  if (nrow(sub_df) == 0) return(NULL)
  if ("cs_max_pip" %in% names(sub_df)) sub_df <- sub_df %>% arrange(desc(as.numeric(cs_max_pip)))
  row1 <- sub_df[1, ]
  pick <- function(col) if (col %in% names(row1)) row1[[col]][1] else NA
  list(
    chrom       = pick("chrom"),
    gene_start  = suppressWarnings(as.numeric(pick("eGene_start"))),
    gene_end    = suppressWarnings(as.numeric(pick("eGene_end"))),
    gene_strand = pick("eGene_strand"),
    gene_symbol = if (gene_col %in% names(row1)) row1[[gene_col]][1] else pick("eGene_symbol"),
    lead_pos    = suppressWarnings(as.numeric(pick("most_likely_snp_pos")))
  )
}

# get_module_betas(): per-celltype beta table for a module, with top gene
# symbol, top SNP rsID, and an exact chi-square p-value where available.
#
# Example:
#   get_module_betas(meta, "M_18448", "cell", "celltype_2", "eGene_symbol")
#   # -> tibble with one row per celltype: beta_val, gene_sym_col, snp_id, ...
get_module_betas <- function(df, mod_id, anno_col, join_col, gene_col) {
  sub_df <- df %>% filter(module == mod_id)
  if (nrow(sub_df) == 0) return(NULL)
  if ("cs_max_pip" %in% names(sub_df)) sub_df <- sub_df %>% arrange(desc(as.numeric(cs_max_pip)))
  sub_df %>%
    distinct(.data[[anno_col]], .keep_all = TRUE) %>%
    mutate(
      final_beta   = if ("cs_top_snp_beta" %in% names(.)) cs_top_snp_beta else if ("most_likely_snp_beta" %in% names(.)) most_likely_snp_beta else beta,
      gene_sym_col = if (gene_col %in% names(.)) .data[[gene_col]] else if ("eGene" %in% names(.)) eGene else "Unknown_Gene",
      mod_id_col   = mod_id,
      snp_id       = if ("most_likely_snp_rsID" %in% names(.)) most_likely_snp_rsID else NA,
      pval_exact   = if ("most_likely_snp_chisq" %in% names(.)) pchisq(as.numeric(most_likely_snp_chisq), df = 1, lower.tail = FALSE) else NA
    ) %>%
    select(all_of(anno_col), beta_val = final_beta, gene_sym_col, mod_id_col, snp_id, pval_exact) %>%
    rename(!!join_col := all_of(anno_col))
}

# --- Main Logic ---
cat("Loading data...\n")
if (show_ref) colors_df <- fread(opt$colors)

gtf_tbl <- NULL
if (!is.null(opt$gtf)) {
  gtf_tbl <- load_gene_table(opt$gtf)
}

annotations_tbl   <- NULL
annotation_tracks <- NULL
if (!is.null(opt$annotations)) {
  anno_files <- trimws(unlist(strsplit(opt$annotations, ",")))
  annotations_tbl <- rbindlist(lapply(anno_files, load_annotation_file),
                               use.names = TRUE, fill = TRUE)
  annotation_tracks <- unique(annotations_tbl$feature_type)
  cat(sprintf("Loaded %d annotation rows across %d track(s): %s\n",
              nrow(annotations_tbl), length(annotation_tracks),
              paste(annotation_tracks, collapse = ", ")))
}

sum_tbl <- if (!is.null(opt$summary_table)) fread(opt$summary_table) else NULL

umap <- fread(opt$umap)
if ("UMAP_0" %in% names(umap) && "UMAP_1" %in% names(umap)) umap <- umap %>% rename(UMAP_2 = UMAP_1, UMAP_1 = UMAP_0)
merged <- if (show_ref) umap %>% left_join(colors_df %>% select(all_of(opt$join_col), hex_color = color_ct2), by = opt$join_col) else umap
if (nrow(merged) > opt$max_cells) { set.seed(42); merged <- merged %>% slice_sample(n = opt$max_cells) }
meta <- fread(opt$master)

cat("Generating plot grid...\n")
plot_list <- list()
n_items <- length(mods)

for (i in seq_len(n_items)) {
  skip_mod <- FALSE
  b_dis <- NA; se_dis <- NA; trt <- NA; snp_dis <- NA
  valid_b <- FALSE; valid_se <- FALSE
  
  if (!is.null(sum_tbl) && mods[i] != "NA") {
    s_match <- sum_tbl %>% filter(module == mods[i])
    if (nrow(s_match) > 0) {
      b_col <- names(s_match)[grepl("^most_likely_beta_disease$", names(s_match), ignore.case = TRUE)][1]
      if (is.na(b_col)) b_col <- names(s_match)[grepl("beta.*disease|disease.*beta", names(s_match), ignore.case = TRUE)][1]
      if (is.na(b_col)) b_col <- names(s_match)[grepl("most_likely_b|most_likely_beta", names(s_match), ignore.case = TRUE) & !grepl("snp", names(s_match), ignore.case = TRUE)][1]
      
      se_col <- names(s_match)[grepl("^most_likely_se_disease$", names(s_match), ignore.case = TRUE)][1]
      if (is.na(se_col)) se_col <- names(s_match)[grepl("se.*disease|disease.*se|se_snp_disease", names(s_match), ignore.case = TRUE)][1]
      
      trt_col <- names(s_match)[grepl("coloc_trait", names(s_match), ignore.case = TRUE)][1]
      snp_dis_col <- names(s_match)[grepl("most_likely_snp", names(s_match), ignore.case = TRUE)][1]
      
      if (is.na(b_col) || is.na(se_col)) {
        cat(sprintf("\n[DIAGNOSTIC] Missing Beta/SE for module %s.\nAvailable columns in your CSV are: %s\n\n",
                    mods[i], paste(names(s_match), collapse = ", ")))
      } else {
        b_dis  <- as.numeric(s_match[[b_col]][1])
        se_dis <- as.numeric(s_match[[se_col]][1])
        trt    <- if (!is.na(trt_col)) as.character(s_match[[trt_col]][1]) else "Unknown_Trait"
        snp_dis <- if (!is.na(snp_dis_col)) as.character(s_match[[snp_dis_col]][1]) else NA
        valid_b  <- length(b_dis) == 1 && !is.na(b_dis)
        valid_se <- length(se_dis) == 1 && !is.na(se_dis)
        if (valid_b && valid_se && b_dis == 0 && se_dis == 0) skip_mod <- TRUE
      }
    }
  }
  
  if (skip_mod) {
    cat(sprintf("Notice: Skipping module %s because Beta and SE are exactly 0.\n", mods[i]))
    next
  }
  if (mods[i] == "NA") next
  
  beta_tbl <- get_module_betas(meta, mods[i], anno_col, opt$join_col, opt$gene_col)
  if (is.null(beta_tbl)) {
    cat(paste("Warning: Module", mods[i], "not found in annotations - skipping panel.\n"))
    next
  }
  
  plot_df <- merged %>% left_join(beta_tbl, by = opt$join_col) %>%
    mutate(beta_val = coalesce(as.numeric(beta_val), 0)) %>%
    arrange(abs(beta_val))
  
  gene_name <- head(na.omit(plot_df$gene_sym_col), 1)
  if (length(gene_name) == 0) gene_name <- "Unknown_Gene"
  
  snp_rsid <- head(na.omit(plot_df$snp_id), 1)
  disp_snp <- if (length(snp_rsid) > 0 && !is.na(snp_rsid) && snp_rsid != "") snp_rsid else snp_dis
  if (is.na(disp_snp) || disp_snp == "") disp_snp <- "Unknown_SNP"
  
  snp_str <- ""; pval_str <- ""
  if (valid_b && valid_se) {
    p_dis <- if (se_dis != 0) 2 * pnorm(-abs(b_dis / se_dis)) else NA
    p_str <- if (!is.na(p_dis)) sprintf(" | P = %.2e", p_dis) else " | P = NA"
    snp_str  <- sprintf("\n%s | %s", trt, disp_snp)
    pval_str <- sprintf(" | Beta = %.3f%s", b_dis, p_str)
  } else {
    snp_str <- sprintf("\n%s", disp_snp)
  }
  
  beta_title <- paste0(opt$name, " | ", gene_name, "\n[", mods[i], "]", snp_str, pval_str)
  
  # Panel order per module: LocusZoom -> Beta -> Coloc Z-scores
  if (has_lz) {
    locus_info <- get_module_locus_info(meta, mods[i], opt$gene_col)
    lz_title <- paste0("LocusZoom | ", opt$name, " | ", gene_name, "\n[", mods[i], "]")
    genes_region <- NULL
    if (!is.null(gtf_tbl) && !is.null(locus_info) &&
        !is.null(locus_info$chrom) && !is.na(locus_info$chrom)) {
      lz_tmp <- fread(lz_files[i])
      win_s <- min(lz_tmp$POS, na.rm = TRUE)
      win_e <- max(lz_tmp$POS, na.rm = TRUE)
      if (!is.na(locus_info$gene_start)) win_s <- min(win_s, locus_info$gene_start)
      if (!is.na(locus_info$gene_end))   win_e <- max(win_e, locus_info$gene_end)
      genes_region <- get_genes_in_region(gtf_tbl, locus_info$chrom, win_s, win_e)
    }
    plot_list[[length(plot_list) + 1]] <- plot_locuszoom(
      lz_files[i], locus_info,
      title             = lz_title,
      genes_df          = genes_region,
      annotations_tbl   = annotations_tbl,
      annotation_tracks = annotation_tracks,
      zoom_window       = opt$zoom_window)
  }
  plot_list[[length(plot_list) + 1]] <- plot_beta(
    plot_df, beta_title, opt$pt_size, opt$join_col,
    show_legend = TRUE, use_raster = use_raster, show_labels = show_labels)
  if (has_z) {
    plot_list[[length(plot_list) + 1]] <- plot_zscore(
      z_files[i], title = paste0("Coloc Z-Scores\n[", mods[i], "]"))
  }
}

# --- Layout Logic ---
panels_per_module <- 1L + as.integer(has_z) + as.integer(has_lz)
cols <- panels_per_module
n_rows <- ceiling(length(plot_list) / cols)
grid_plot <- wrap_plots(plot_list, ncol = cols)

if (show_ref) {
  p_ref <- plot_ref(merged, paste(opt$name, "Reference"), opt$join_col, opt$pt_size, use_raster)
  final_plot <- p_ref | grid_plot
  final_plot <- final_plot + plot_layout(widths = c(1, cols))
  plot_width <- (cols + 1) * 6
} else {
  final_plot <- grid_plot
  plot_width <- cols * 6
}

base_height <- 4.5
plot_height <- max(base_height * n_rows, 5)

# --- Save Logic ---
out_base <- sub("\\.png$|\\.pdf$", "", opt$out, ignore.case = TRUE)
if (save_pdf) {
  out_file <- paste0(out_base, ".pdf")
  cat(paste("Saving as vector PDF with rasterized points to:", out_file, "\n"))
  ggsave(out_file, plot = final_plot, width = plot_width, height = plot_height,
         device = cairo_pdf, limitsize = FALSE)
} else {
  out_file <- paste0(out_base, ".png")
  cat(paste("Saving as fast PNG to:", out_file, "\n"))
  ggsave(out_file, plot = final_plot, width = plot_width, height = plot_height,
         dpi = 300, bg = "white", device = ragg::agg_png, limitsize = FALSE)
}
