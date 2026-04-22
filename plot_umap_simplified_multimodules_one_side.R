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

# Load modularised helpers from R/ (same directory as this script).
# Works under Rscript (uses --file=...), source() (uses sys.frame(1)$ofile),
# or interactive runs (falls back to getwd()).
script_dir <- tryCatch(
  dirname(normalizePath(sys.frame(1)$ofile, mustWork = FALSE)),
  error = function(e) getwd()
)
if (!dir.exists(file.path(script_dir, "R"))) {
  .cli_args <- commandArgs(trailingOnly = FALSE)
  .file_arg <- sub("--file=", "", .cli_args[grepl("^--file=", .cli_args)])
  if (length(.file_arg) == 1 && nzchar(.file_arg))
    script_dir <- dirname(normalizePath(.file_arg))
}
for (f in sort(list.files(file.path(script_dir, "R"),
                          pattern = "\\.R$", full.names = TRUE))) {
  source(f, local = FALSE)
}

# --- CLI Options ---
option_list <- list(
  make_option("--colors", type = "character", default = NULL, help = "Path to palette TSV (only required if --show_ref is used)"),
  
  # All optional. At minimum you need EITHER --lz_files alone (in which case
  # the labelled / anchor SNP becomes the top-P SNP per file), OR the
  # UMAP + master + modules triple (classic per-module plotting).
  make_option("--umap", type = "character", default = NULL, help = "Path to UMAP TSV (optional). Omit to skip the Beta UMAP panel. Pair with --master to enable it."),
  make_option("--master", type = "character", default = NULL, help = "Path to Master Annotations TSV (optional). Used for per-(cell, gene) betas, gene-track symbol lookups and the module-level most-likely SNP."),
  make_option("--modules", type = "character", default = NULL, help = "Comma-separated Module IDs (optional). When omitted the script treats each --lz_files entry as one implicit module and uses the top -log10(P) SNP as the anchor."),
  
  # Optional per-module side panels
  make_option("--summary_table", type = "character", default = NULL, help = "Path to Summary Table (for disease titles)"),
  make_option("--z_files", type = "character", default = NULL, help = "Comma-separated paths to Z-Score CSVs (use NA for missing)"),
  make_option("--lz_files", type = "character", default = NULL, help = "Comma-separated paths to LocusZoom CSVs (use NA for missing). Expected columns: CHR,CELL,GENE,POS,P"),
  make_option("--ld_files", type = "character", default = NULL, help = "Comma-separated paths to plink2 --export A .raw genotype files (use NA for missing). One per module; drives the r^2-based LD colouring of the LocusZoom points."),
  make_option("--name", type = "character", default = "Pop", help = "Display name used in panel titles"),
  
  # General Settings
  make_option("--join_col", type = "character", default = "celltype_2", help = "Column name in UMAP for celltypes [default: %default]"),
  make_option("--anno_join_col", type = "character", default = NULL, help = "Column name in Master Annotations (if different from join_col)"),
  make_option("--gene_col", type = "character", default = "eGene_symbol", help = "Column name for gene symbol [default: %default]"),
  make_option("--gtf", type = "character", default = NULL, help = "Optional GENCODE GTF (e.g. gencode.v49.annotation.gtf) or a pre-built .coding_genes.tsv. When supplied, the LocusZoom gene track shows every protein-coding gene in the window."),
  make_option("--annotations", type = "character", default = NULL, help = "Optional comma-separated annotation file(s) to show in a zoom panel around the lead SNP. BED-style (chrom, start, end, feature_type [, score]) or full 9-col GFF. A numeric 5th column turns the lane into a continuous profile."),
  make_option("--zoom_window", type = "numeric", default = 5000, help = "Half-width in bp for the zoom annotation panel around the lead SNP [default: %default]"),
  make_option("--cell", type = "character", default = NULL, help = "Optional comma-separated list of cell types to keep. Default: plot every (cell, gene) pair a module has."),
  make_option("--gene", type = "character", default = NULL, help = "Optional comma-separated list of gene identifiers to keep (matched against eGene_symbol or bare eGene ENSG)."),
  
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
has_umap    <- !is.null(opt$umap)
has_master  <- !is.null(opt$master)
has_modules <- !is.null(opt$modules)
has_z       <- !is.null(opt$z_files)
has_lz      <- !is.null(opt$lz_files)
has_ld      <- !is.null(opt$ld_files)

# Beta UMAP needs both the UMAP coords and the master annotation.
has_beta_panel <- has_umap && has_master
if (has_umap && !has_master)
  cat("Note: --umap provided without --master; the Beta UMAP panel is skipped (master annotations hold the per-(cell, gene) betas).\n")
if (has_master && !has_umap)
  cat("Note: --master provided without --umap; the Beta UMAP panel is skipped (no embedding to plot on).\n")

# Need at least something to plot.
if (!has_lz && !has_beta_panel)
  stop("Nothing to plot: please provide either --lz_files, or all of --umap + --master + --modules.")

if (show_ref && is.null(opt$colors))
  stop("You must provide the --colors file if you enable --show_ref.")
if (show_ref && !has_umap)
  stop("--show_ref requires --umap (there is no UMAP to draw the reference on).")

cell_filter <- if (!is.null(opt[["cell"]])) trimws(unlist(strsplit(opt[["cell"]], ","))) else NULL
gene_filter <- if (!is.null(opt[["gene"]])) trimws(unlist(strsplit(opt[["gene"]], ","))) else NULL

# Establish the module list (either explicit via --modules or synthesised
# from --lz_files, one implicit module per file). When only --modules is
# given and no --lz_files, we fall back to the master-annotation flow.
if (has_modules) {
  mods <- trimws(unlist(strsplit(opt$modules, ",")))
} else if (has_lz) {
  tmp_lz <- trimws(unlist(strsplit(opt$lz_files, ",")))
  mods <- sprintf("auto_%s",
                  tools::file_path_sans_ext(basename(tmp_lz)))
  cat(sprintf("--modules not provided; treating each of the %d LZ file(s) as one implicit module: %s\n",
              length(mods), paste(mods, collapse = ", ")))
} else {
  stop("Provide --modules or --lz_files; cannot derive modules otherwise.")
}
n_items <- length(mods)

if (has_z) {
  z_files <- trimws(unlist(strsplit(opt$z_files, ",")))
  if (length(z_files) != n_items) stop("--z_files must have the same number of items as --modules (or --lz_files when --modules is omitted).")
}
if (has_lz) {
  lz_files <- trimws(unlist(strsplit(opt$lz_files, ",")))
  if (length(lz_files) != n_items) stop("--lz_files must have the same number of items as --modules (or --lz_files itself defines the module count when --modules is omitted).")
}
if (has_ld) {
  ld_files <- trimws(unlist(strsplit(opt$ld_files, ",")))
  if (length(ld_files) != n_items) stop("--ld_files must have the same number of items as --modules / --lz_files.")
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

umap <- NULL; merged <- NULL; meta <- NULL
if (has_umap) {
  umap <- fread(opt$umap)
  if ("UMAP_0" %in% names(umap) && "UMAP_1" %in% names(umap))
    umap <- umap %>% rename(UMAP_2 = UMAP_1, UMAP_1 = UMAP_0)
  merged <- if (show_ref)
              umap %>% left_join(colors_df %>% select(all_of(opt$join_col),
                                                       hex_color = color_ct2),
                                 by = opt$join_col)
            else umap
  if (nrow(merged) > opt$max_cells) { set.seed(42); merged <- merged %>% slice_sample(n = opt$max_cells) }
}
if (has_master) meta <- fread(opt$master)

cat("Generating plot grid...\n")
module_composites <- list()  # one patchwork per module (CS grid + merged box)
total_cs_rows     <- 0L
# Column layout: LZ, optional Beta UMAP, optional Coloc Z.
cols <- as.integer(has_lz) + as.integer(has_beta_panel) + as.integer(has_z)
if (cols == 0)
  stop("Nothing to plot (no LZ, no Beta panel, no Z).")

for (i in seq_len(n_items)) {
  mod_id <- mods[i]
  if (mod_id == "NA") next
  
  # --- disease-level summary (same for all CS rows of this module)
  skip_mod <- FALSE
  b_dis <- NA; se_dis <- NA; trt <- NA; snp_dis <- NA
  valid_b <- FALSE; valid_se <- FALSE
  if (has_master && !is.null(sum_tbl)) {
    s_match <- sum_tbl %>% filter(module == mod_id)
    if (nrow(s_match) > 0) {
      b_col <- names(s_match)[grepl("^most_likely_beta_disease$", names(s_match), ignore.case = TRUE)][1]
      if (is.na(b_col)) b_col <- names(s_match)[grepl("beta.*disease|disease.*beta", names(s_match), ignore.case = TRUE)][1]
      if (is.na(b_col)) b_col <- names(s_match)[grepl("most_likely_b|most_likely_beta", names(s_match), ignore.case = TRUE) & !grepl("snp", names(s_match), ignore.case = TRUE)][1]
      
      se_col <- names(s_match)[grepl("^most_likely_se_disease$", names(s_match), ignore.case = TRUE)][1]
      if (is.na(se_col)) se_col <- names(s_match)[grepl("se.*disease|disease.*se|se_snp_disease", names(s_match), ignore.case = TRUE)][1]
      trt_col     <- names(s_match)[grepl("coloc_trait", names(s_match), ignore.case = TRUE)][1]
      snp_dis_col <- names(s_match)[grepl("most_likely_snp", names(s_match), ignore.case = TRUE)][1]
      
      if (is.na(b_col) || is.na(se_col)) {
        cat(sprintf("\n[DIAGNOSTIC] Missing Beta/SE for module %s.\nAvailable columns in your CSV are: %s\n\n",
                    mod_id, paste(names(s_match), collapse = ", ")))
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
    cat(sprintf("Notice: Skipping module %s because disease Beta and SE are exactly 0.\n", mod_id))
    next
  }
  
  # --- enumerate credible sets for this module (from LZ CS column when
  # available, otherwise fall back to (cell, gene) in master). Each CS
  # produces one grid row of [LocusZoom | Beta UMAP | optional Coloc Z].
  lz_path  <- if (has_lz) lz_files[i] else NULL
  cs_rows  <- module_cs_list(lz_path, meta, mod_id, anno_col,
                             cell_filter = cell_filter,
                             gene_filter = gene_filter)
  if (nrow(cs_rows) == 0) {
    cat(sprintf("Warning: Module %s has no credible sets matching the filters - skipping.\n", mod_id))
    next
  }
  
  # LD matrix (r^2) for this module; all CSs share the same genotype file.
  ld_mat <- if (has_ld) load_ld_matrix(ld_files[i]) else NULL
  
  module_plots <- list()
  lead_bag     <- numeric(0)   # lead SNP positions, for the merged box
  cs_labels    <- character(0)
  chr_for_mod  <- NA_character_
  n_cs_mod     <- nrow(cs_rows)  # total CSs kept, used to place the zoom
                                 # connector only under the last LocusZoom
  
  # Module-level most-likely SNP -- the shared anchor for the whole module.
  # Every LocusZoom diamond is drawn at this SNP and the LD r^2 colouring
  # is computed against it, so all LZ panels of a module share identical
  # LD patterns. Also re-used as the distinct marker on the merged box.
  #
  # Source priority:
  #   1. master annotation's `most_likely_snp` (picks the row with highest
  #      cs_max_pip for this module),
  #   2. top -log10(P) SNP across the module's LZ file when no master is
  #      provided.
  mod_master_snp_pos <- NA_real_
  mod_master_snp_lab <- NA_character_
  mod_master_snp_id  <- NA_character_   # "chr<N>:<pos>:<ref>:<alt>" for LD lookup
  if (has_master) {
    mrows <- meta %>% filter(module == mod_id)
    if ("cs_max_pip" %in% names(mrows))
      mrows <- mrows %>% arrange(desc(as.numeric(cs_max_pip)))
    if (nrow(mrows) > 0) {
      if ("most_likely_snp_pos" %in% names(mrows))
        mod_master_snp_pos <- suppressWarnings(as.numeric(mrows[["most_likely_snp_pos"]][1]))
      if ("most_likely_snp" %in% names(mrows) &&
          !is.na(mrows[["most_likely_snp"]][1]) &&
          nzchar(mrows[["most_likely_snp"]][1])) {
        mod_master_snp_lab <- as.character(mrows[["most_likely_snp"]][1])
        mod_master_snp_id  <- mod_master_snp_lab
      } else if (!is.na(mod_master_snp_pos)) {
        mod_master_snp_lab <- sprintf("most-likely SNP | %s bp",
                                      format(mod_master_snp_pos, big.mark = ",", scientific = FALSE))
      }
      # If we have a pos but no full chr:pos:ref:alt, fall back to a
      # "chr<N>:<pos>:N:N" token; the LD lookup matches by <chr>:<pos> only.
      if (is.na(mod_master_snp_id) && !is.na(mod_master_snp_pos)) {
        mod_chr <- if ("most_likely_snp_chrom" %in% names(mrows))
                     sub("^chr", "", as.character(mrows[["most_likely_snp_chrom"]][1]),
                         ignore.case = TRUE)
                   else if ("chrom" %in% names(mrows))
                     sub("^chr", "", as.character(mrows[["chrom"]][1]), ignore.case = TRUE)
                   else NA_character_
        if (!is.na(mod_chr))
          mod_master_snp_id <- sprintf("chr%s:%d:N:N", mod_chr,
                                        as.integer(mod_master_snp_pos))
      }
    }
  }
  # No master annotation -> derive the module anchor from the LZ file.
  # Pick the SNP with the highest -log10(P). When the LZ file includes a
  # CS column, try to reuse the chr:pos:ref:alt token embedded in that CS
  # for consistent labelling; otherwise synthesise "chr<N>:<pos>:N:N".
  if (is.na(mod_master_snp_pos) && has_lz) {
    lz_tmp <- load_lz_file(lz_files[i])
    if (!is.null(lz_tmp) && nrow(lz_tmp) > 0) {
      lz_tmp[, logp := -log10(pmax(as.numeric(P), 1e-300))]
      top_idx <- which.max(lz_tmp$logp)
      mod_master_snp_pos <- as.numeric(lz_tmp$POS[top_idx])
      top_chr <- sub("^chr","", as.character(lz_tmp$CHR[top_idx]), ignore.case = TRUE)
      if ("CS" %in% names(lz_tmp) && !is.na(lz_tmp$CS[top_idx])) {
        m <- regmatches(lz_tmp$CS[top_idx],
                        regexpr("chr[0-9XYM]+:\\d+:[ACGTN]+:[ACGTN]+",
                                lz_tmp$CS[top_idx], ignore.case = TRUE))
        if (length(m) == 1 && nzchar(m)) {
          parts <- strsplit(sub("^chr","",m), ":", fixed = TRUE)[[1]]
          if (length(parts) >= 2 &&
              as.integer(parts[2]) == as.integer(mod_master_snp_pos)) {
            mod_master_snp_lab <- m
            mod_master_snp_id  <- m
          }
        }
      }
      if (is.na(mod_master_snp_lab))
        mod_master_snp_lab <- sprintf("chr%s:%d", top_chr,
                                      as.integer(mod_master_snp_pos))
      if (is.na(mod_master_snp_id))
        mod_master_snp_id  <- sprintf("chr%s:%d:N:N", top_chr,
                                      as.integer(mod_master_snp_pos))
      cat(sprintf("  [%s] no master; using top-P SNP %s as module anchor.\n",
                  mod_id, mod_master_snp_lab))
    }
  }
  
  for (p_idx in seq_len(nrow(cs_rows))) {
    this_cs   <- cs_rows$cs[p_idx]
    this_cell <- cs_rows$cell[p_idx]
    this_gene <- cs_rows$eGene[p_idx]
    this_sym  <- cs_rows$eGene_symbol[p_idx]
    if (is.na(this_sym) || this_sym == "") this_sym <- this_gene
    
    # Beta UMAP is only produced when both --umap and --master are set.
    plot_df <- NULL
    if (has_beta_panel) {
      beta_tbl <- get_module_betas(meta, mod_id, focal_gene = this_gene,
                                   anno_col, opt$join_col, opt$gene_col)
      if (is.null(beta_tbl)) {
        cat(sprintf("Warning: No beta rows for module %s, gene %s - skipping Beta UMAP for this CS.\n",
                    mod_id, this_gene))
      } else {
        plot_df <- merged %>% left_join(beta_tbl, by = opt$join_col) %>%
          mutate(beta_val = coalesce(as.numeric(beta_val), 0)) %>%
          arrange(abs(beta_val))
      }
    }
    
    # Title shows the MODULE-level master "chr:pos:a1:a2" so it matches the
    # diamond/LD anchor drawn on every LocusZoom of this module. Fallbacks:
    # CS-encoded lead -> master rsID -> summary-table disease SNP.
    disp_snp <- if (!is.na(mod_master_snp_lab) && nzchar(mod_master_snp_lab))
                  mod_master_snp_lab
                else NA_character_
    if (is.na(disp_snp) && !is.na(this_cs)) {
      m <- regmatches(this_cs,
                      regexpr("chr[0-9XYM]+:\\d+:[ACGTN]+:[ACGTN]+",
                              this_cs, ignore.case = TRUE))
      if (length(m) == 1 && nzchar(m)) disp_snp <- m
    }
    if (is.na(disp_snp) && !is.null(plot_df)) {
      snp_rsid <- head(na.omit(plot_df$snp_id), 1)
      if (length(snp_rsid) > 0 && !is.na(snp_rsid) && snp_rsid != "")
        disp_snp <- snp_rsid else disp_snp <- snp_dis
    }
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
    beta_title <- paste0(opt$name, " | ", this_cell, " | ", this_sym,
                         "\n[", mod_id, "]", snp_str, pval_str)
    
    if (has_lz) {
      # Build locus info. If the master annotation is available, use it for
      # gene coords / eGene symbol / strand; otherwise start from an empty
      # list and rely on the CS-derived chrom + lead_pos alone.
      locus_info <- if (has_master)
                      get_module_locus_info(meta, mod_id,
                                            focal_cell = this_cell,
                                            focal_gene = this_gene,
                                            anno_col   = anno_col,
                                            gene_col   = opt$gene_col)
                    else NULL
      if (is.null(locus_info)) locus_info <- list()
      if (!is.na(cs_rows$lead_pos[p_idx])) {
        locus_info$lead_pos <- cs_rows$lead_pos[p_idx]
        if (is.null(locus_info$chrom) || is.na(locus_info$chrom))
          locus_info$chrom <- cs_rows$chrom[p_idx]
      }
      # Still no chrom? Fall back to the LZ file's first CHR entry.
      if (is.null(locus_info$chrom) || is.na(locus_info$chrom)) {
        lz_peek <- load_lz_file(lz_files[i])
        if (!is.null(lz_peek) && nrow(lz_peek) > 0)
          locus_info$chrom <- sub("^chr","", as.character(lz_peek$CHR[1]),
                                  ignore.case = TRUE)
      }
      lz_title <- paste0("LocusZoom | ", opt$name, " | ", this_cell, " | ",
                         this_sym, "\n[", mod_id, "]")
      
      # Region for the gene track: union of SNP extents (filtered) + eGene body.
      genes_region <- NULL
      if (!is.null(gtf_tbl) && !is.null(locus_info) &&
          !is.null(locus_info$chrom) && !is.na(locus_info$chrom)) {
        lz_tmp <- load_lz_file(lz_files[i])
        if (!is.null(lz_tmp)) {
          if ("CELL" %in% names(lz_tmp)) lz_tmp <- lz_tmp[as.character(CELL) == this_cell]
          if ("GENE" %in% names(lz_tmp)) lz_tmp <- lz_tmp[sub("\\.\\d+$","",as.character(GENE)) == this_gene]
          if ("CS"   %in% names(lz_tmp) && !is.na(this_cs)) lz_tmp <- lz_tmp[as.character(CS) == this_cs]
        } else {
          lz_tmp <- data.table(POS = integer(0))
        }
        if (nrow(lz_tmp) > 0) {
          win_s <- min(lz_tmp$POS, na.rm = TRUE)
          win_e <- max(lz_tmp$POS, na.rm = TRUE)
          if (!is.null(locus_info$gene_start) && !is.na(locus_info$gene_start))
            win_s <- min(win_s, locus_info$gene_start)
          if (!is.null(locus_info$gene_end) && !is.na(locus_info$gene_end))
            win_e <- max(win_e, locus_info$gene_end)
          genes_region <- get_genes_in_region(gtf_tbl, locus_info$chrom,
                                              win_s, win_e)
        }
      }
      
      # LD r^2 vector anchored on the MODULE-level most-likely SNP (same
      # across every CS of this module). This way both LocusZooms show
      # identical LD colouring; the diamond is also pinned to the same
      # position via `force_lead_pos` below.
      ld_vec_cs <- NULL
      if (!is.null(ld_mat)) {
        ld_anchor_id <- if (!is.na(mod_master_snp_id)) mod_master_snp_id else {
          # Fallback (only used when master annotation lacks a most-likely
          # SNP): CS-encoded lead, then locus_info, then constructed chr:pos.
          lead_id <- NA_character_
          if (!is.na(this_cs)) {
            m <- regmatches(this_cs,
                            regexpr("chr[0-9XYM]+:\\d+:[ACGTN]+:[ACGTN]+",
                                    this_cs, ignore.case = TRUE))
            if (length(m) == 1 && nzchar(m)) lead_id <- m
          }
          if (is.na(lead_id) && !is.null(locus_info) &&
              !is.null(locus_info$chrom) && !is.na(locus_info$lead_pos)) {
            lead_id <- sprintf("chr%s:%d:N:N",
                               sub("^chr","",locus_info$chrom,ignore.case=TRUE),
                               as.integer(locus_info$lead_pos))
          }
          lead_id
        }
        ld_vec_cs <- ld_vec_for_lead(ld_mat, ld_anchor_id)
        if (!is.null(ld_vec_cs)) {
          bin_counts <- table(bin_ld_r2(ld_vec_cs))
          cat(sprintf("  LD bins for %s (%s vs module anchor %s): %s\n",
                      mod_id,
                      if (!is.na(this_cs)) sub(".*::([^:]+)$", "\\1", this_cs) else paste0("CS", p_idx),
                      ld_anchor_id,
                      paste(sprintf("%s=%d", names(bin_counts), as.integer(bin_counts)),
                            collapse = ", ")))
        }
      }
      
      # Only the last CS row of a module draws the zoom connector diagonals;
      # earlier rows omit it so the module's figure has one zoom-out fan at
      # the bottom, flowing into the merged annotation box.
      is_last_cs <- (p_idx == n_cs_mod)
      
      module_plots[[length(module_plots) + 1]] <- plot_locuszoom(
        lz_files[i], locus_info,
        title                  = lz_title,
        genes_df               = genes_region,
        annotations_tbl        = annotations_tbl,
        annotation_tracks      = annotation_tracks,
        zoom_window            = opt$zoom_window,
        lz_filter_cell         = this_cell,
        lz_filter_gene         = this_gene,
        lz_filter_cs           = if (!is.na(this_cs)) this_cs else NULL,
        include_zoom_panel     = FALSE,  # zoom panel is now one merged box per module
        include_zoom_connector = is_last_cs,
        ld_vec                 = ld_vec_cs,
        force_lead_pos         = mod_master_snp_pos)
      
      if (!is.na(cs_rows$lead_pos[p_idx])) {
        lead_bag <- c(lead_bag, cs_rows$lead_pos[p_idx])
        # Build a compact, human-readable label "cell | symbol | L<n>" so the
        # merged box's dashed line points back at the right LocusZoom row.
        l_tok <- if (!is.na(this_cs)) sub(".*::([^:]+)$", "\\1", this_cs) else NA_character_
        short_sym <- if (!is.na(this_sym) && nchar(this_sym) > 0) this_sym else this_gene
        lab <- paste(c(this_cell, short_sym,
                       if (!is.na(l_tok) && nchar(l_tok) > 0 && l_tok != this_cs) l_tok else NULL),
                     collapse = " | ")
        cs_labels <- c(cs_labels, lab)
      }
      if (is.na(chr_for_mod) && !is.null(locus_info) &&
          !is.null(locus_info$chrom) && !is.na(locus_info$chrom)) {
        chr_for_mod <- sub("^chr", "", locus_info$chrom, ignore.case = TRUE)
      }
    }
    if (has_beta_panel && !is.null(plot_df)) {
      module_plots[[length(module_plots) + 1]] <- plot_beta(
        plot_df, beta_title, opt$pt_size, opt$join_col,
        show_legend = TRUE, use_raster = use_raster, show_labels = show_labels)
    } else if (has_beta_panel) {
      # Column must stay aligned -- insert a blank placeholder when this
      # CS had no master-derived beta rows.
      module_plots[[length(module_plots) + 1]] <- plot_spacer()
    }
    if (has_z) {
      module_plots[[length(module_plots) + 1]] <- plot_zscore(
        z_files[i],
        title = paste0("Coloc Z-Scores | ", this_cell, " | ", this_sym,
                       "\n[", mod_id, "]"))
    }
  }
  
  if (length(module_plots) == 0) next
  
  n_cs_for_mod <- length(module_plots) %/% cols
  total_cs_rows <- total_cs_rows + n_cs_for_mod
  cs_grid <- wrap_plots(module_plots, ncol = cols)
  
  # Merged, full-width annotation box for the whole module.
  if (!is.null(annotations_tbl) && !is.null(annotation_tracks) &&
      length(lead_bag) > 0 && !is.na(chr_for_mod)) {
    merged_box <- build_merged_zoom_box(
      chr_label         = chr_for_mod,
      lead_positions    = lead_bag,
      cs_labels         = cs_labels,
      annotations_tbl   = annotations_tbl,
      annotation_tracks = annotation_tracks,
      zoom_window       = opt$zoom_window,
      title = sprintf("Merged annotations | [%s] (%d credible set%s)",
                      mod_id, length(lead_bag),
                      if (length(lead_bag) > 1) "s" else ""),
      extra_lead_pos    = mod_master_snp_pos,
      extra_lead_label  = if (!is.na(mod_master_snp_lab))
                            paste0("module most-likely | ", mod_master_snp_lab)
                          else NA_character_)
    # Heights: each CS row ~5 units, merged box ~2 units. Keeps the box
    # readable without dominating the figure. Wrapping `cs_grid` with
    # `wrap_elements()` keeps it atomic so `/` stacks two blocks rather
    # than flattening the inner panels (which breaks when cols == 1).
    merged_h <- 2
    composite <- wrap_elements(cs_grid) / merged_box +
      plot_layout(heights = c(n_cs_for_mod * 5, merged_h))
    module_composites[[length(module_composites) + 1]] <- composite
    total_cs_rows <- total_cs_rows + (merged_h / 5)  # ~0.4 row-equivalents
  } else {
    module_composites[[length(module_composites) + 1]] <- cs_grid
  }
}

# --- Layout Logic ---
if (length(module_composites) == 0) stop("Nothing to plot - every module was skipped.")
grid_plot <- if (length(module_composites) == 1) {
               module_composites[[1]]
             } else {
               wrap_plots(module_composites, ncol = 1)
             }
n_rows <- max(1, total_cs_rows)

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
