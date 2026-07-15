#!/usr/bin/env Rscript
# Single-population UMAP / optional LocusZoom / optional Coloc-beta plot.
# Per credible set the grid row is: [optional LocusZoom] -> UMAP (beta or expression) -> [optional Coloc Beta].

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
# Rscript sets --file=<path>; sys.frame(1)$ofile is not reliable at top level,
# so always derive script_dir from --file= first (works when cwd is a parent folder).
.cli_args <- commandArgs(trailingOnly = FALSE)
.file_arg <- sub("^--file=", "", .cli_args[grepl("^--file=", .cli_args)])
script_dir <- if (length(.file_arg) >= 1L && nzchar(.file_arg[[1L]])) {
  dirname(normalizePath(.file_arg[[1L]], mustWork = FALSE))
} else {
  tryCatch(
    dirname(normalizePath(sys.frame(1)$ofile, mustWork = FALSE)),
    error = function(e) normalizePath(getwd(), mustWork = FALSE)
  )
}
if (!dir.exists(file.path(script_dir, "R"))) {
  stop("Cannot find R/ next to this script. script_dir=", script_dir,
       "\nRun from the project folder or use: Rscript path/to/plot_umap_simplified_multimodules_one_side.R")
}
for (f in sort(list.files(file.path(script_dir, "R"),
                          pattern = "\\.R$", full.names = TRUE))) {
  source(f, local = FALSE)
}
# Ensure plot_beta() is the version shipped with this script (e.g. label_cells_only).
source(file.path(script_dir, "R", "11_plot_beta.R"), local = FALSE)
# Same for beta-scatter helpers (build_beta_column, plot_beta_scatter_df).
source(file.path(script_dir, "R", "12_plot_zscore.R"), local = FALSE)
if (!exists("build_beta_column", mode = "function"))
  stop("R/12_plot_zscore.R must define build_beta_column(); refresh that file from the repo.")

# --- CLI Options ---
option_list <- list(
  make_option("--colors", type = "character", default = NULL, help = "Path to palette TSV (only required if --show_ref is used)"),
  
  # All optional. At minimum you need EITHER --lz_files alone (in which case
  # the labelled / anchor SNP becomes the top-P SNP per file), OR the
  # UMAP + master + modules triple (classic per-module plotting).
  # Register --umap_sample_n before --umap so getopt does not treat the
  # integer default -1 as the path for --umap (shared prefix).
  make_option("--max_cells", type = "numeric", default = 250000,
              help = "When --umap_sample_n is not set: subsample the UMAP to this many rows only if the input has more rows than this value [default: %default]"),
  make_option("--umap_sample_n", type = "integer", default = -1L,
              help = "If >= 1, randomly subsample UMAP rows to at most this many cells (set.seed 42) for smaller PDFs, even when total rows is below --max_cells. Use -1 to rely on --max_cells only [default: %default]"),
  make_option("--umap", type = "character", default = NULL, help = "Path to UMAP TSV/CSV (optional). Narrow table: coords + cell type (--umap_color_mode beta, needs --master). Wide table: coords + cell type + one gene expression column (--umap_color_mode expression)."),
  make_option("--umap_color_mode", type = "character", default = "beta",
              help = "Middle UMAP panel colouring: 'beta' (per-cell-type QTL beta from --master, default) or 'expression' (per-cell values from a gene column in --umap) [default: %default]"),
  make_option("--log1p", action = "store_true", default = FALSE,
              help = "Apply log1p to expression before colouring (expression mode only)."),
  make_option("--clip_quantile", type = "double", default = 0.99,
              help = "Clip expression colour scale at this upper quantile of non-zero cells; set to 1 to disable (expression mode) [default: %default]"),
  make_option("--expr_palette", type = "character", default = "yellowred",
              help = "Expression colour ramp: yellowred | brightyellowred | viridis | magma | inferno | plasma | turbo | cividis [default: %default]"),
  make_option("--master", type = "character", default = NULL, help = "Path to Master Annotations TSV (optional). Used for per-(cell, gene) betas, gene-track symbol lookups and the module-level most-likely SNP."),
  make_option("--modules", type = "character", default = NULL, help = "Comma-separated Module IDs (optional). When omitted the script treats each --lz_files entry as one implicit module and uses the top -log10(P) SNP as the anchor."),
  # Manual / single-CS mode: bypass --master and --modules by spelling out
  # the credible-set string + a single focal beta. Useful when you only have
  # an ad-hoc beta for one (cell, gene, snp) and don't want to materialise a
  # master annotation TSV just to plot it.
  make_option("--cs_name", type = "character", default = NULL,
              help = "Manual credible-set identifier, e.g. 'chr22::study::NK_CD16_adaptive:ENSG00000015475::chr22:17604981:C:T::L1'. When set together with --beta_cs the script runs in MANUAL mode: --master and --modules are ignored (a notice is printed if they are also set) and the cell / gene / chromosome / lead SNP are parsed from this string. --cell, when provided, overrides the cell parsed from --cs_name."),
  make_option("--beta_cs", type = "numeric", default = NA_real_,
              help = "Beta value used to colour the focal cell type on the Beta UMAP when --cs_name is set (manual mode). The diverging colour scale is centred at 0 and its extent is |--beta_cs|."),
  
  # Optional per-module side panels
  make_option("--summary_table", type = "character", default = NULL, help = "Path to Summary Table (for disease titles)"),
  make_option("--beta_files", type = "character", default = NULL, help = "Comma-separated paths to Coloc Beta CSVs (use NA for missing). Columns: beta_qtl vs beta_disease (legacy z_qtl / z_disease / z_icd10* still accepted); optional snp, cs_qtl. Multiple distinct cs_qtl values -> stacked beta panels; if a grid row's CS matches cs_qtl, only that subset is plotted for that row."),
  make_option("--z_files", type = "character", default = NULL, help = "Deprecated alias for --beta_files (kept for backward compatibility)."),
  make_option("--beta_axes", type = "character", default = NULL,
              help = "Coloc Beta panel axis scaling: 'linked' (default) uses the same symmetric limits on QTL-beta and disease-beta from max(|beta|) over both columns with coord_fixed; 'independent' sets each axis from max(|beta|) in that column only (coord_cartesian)."),
  make_option("--z_axes", type = "character", default = NULL,
              help = "Deprecated alias for --beta_axes (kept for backward compatibility)."),
  make_option("--lz_files", type = "character", default = NULL, help = "Optional comma-separated LocusZoom CSVs (use NA for missing) when you want LocusZoom panels. Omit when using --umap + --master + --modules: credible sets come from the master table only. Expected columns: CHR,CELL,GENE,POS,P. Coloc Beta tables belong in --beta_files."),
  make_option("--lz_layout", type = "character", default = "stacked",
              help = "LocusZoom layout: 'stacked' (one row per credible set) or 'merged' (single panel, all CS overlaid; LZ-only, no beta UMAP or Z). Merged adds a bottom strip (genomic span per CS, same Mb axis)."),
  make_option("--lz_xlim", type = "character", default = "context",
              help = "LocusZoom x-axis span: 'context' (default) widens to focal eGene + genes from --gtf; 'snp' uses only min/max association POS in the panel (+ 2%% pad). Gene bodies may be clipped when using 'snp'."),
  make_option("--ld_files", type = "character", default = NULL, help = "Comma-separated paths to plink2 --export A .raw genotype files (use NA for missing). One per module; drives the r^2-based LD colouring of the LocusZoom points."),
  make_option("--name", type = "character", default = "Pop", help = "Display name used in panel titles"),
  
  # General Settings
  make_option("--join_col", type = "character", default = "celltype_3", help = "Column name in UMAP for celltypes used for JOIN / FILTER logic: --cell matching, joining to --master, and identifying focal cells from --cs_name [default: %default]."),
  make_option("--anno_join_col", type = "character", default = NULL, help = "Column name in Master Annotations (if different from join_col)"),
  make_option("--label_col", type = "character", default = NULL, help = "Column name in UMAP used ONLY for the cell-type LABELS drawn on the Beta / Reference UMAPs (centroid text). Defaults to --join_col. Set this to e.g. 'celltype_3' to label cells at a different granularity from the one used to join with --master / --cell. When --label_col differs from --join_col, only label-column groups that contain at least one focal cell (per --join_col) are labelled."),
  make_option("--gene_col", type = "character", default = "eGene_symbol", help = "Column name for gene symbol [default: %default]"),
  make_option("--gtf", type = "character", default = NULL, help = "Optional GENCODE GTF (e.g. gencode.v49.annotation.gtf) or a pre-built .coding_genes.tsv. When supplied, the LocusZoom gene track shows every protein-coding gene in the window."),
  make_option("--annotations", type = "character", default = NULL, help = "Optional comma-separated annotation file(s) to show in a zoom panel around the lead SNP. BED-style (chrom, start, end, feature_type [, score]) or full 9-col GFF. A numeric 5th column turns the lane into a continuous profile."),
  make_option("--zoom_window", type = "numeric", default = 5000, help = "Half-width in bp for the zoom annotation panel around the lead SNP [default: %default]"),
  make_option("--cell", type = "character", default = NULL, help = "Optional comma-separated list of cell types to keep for credible-set rows (module_cs_list). On the beta UMAP, only these types use the beta colour scale; other cells stay grey (--join_col must match the UMAP). Default: all (cell, gene) pairs and all cell types coloured by beta."),
  make_option("--cells", type = "character", default = NULL, help = "Alias for --cell (same comma-separated list). If both are set, --cell wins."),
  make_option("--gene", type = "character", default = NULL, help = "Optional comma-separated list of gene identifiers to keep (matched against eGene_symbol or bare eGene ENSG)."),
  
  make_option("--out", type = "character", default = "umap_plot", help = "Output filename prefix"),
  make_option("--pt_size", type = "numeric", default = 0.25, help = "UMAP point scale (default %default). Raster PDFs use geom_scattermore: larger values map to visibly bigger pixels; try 0.4–1.2. With --no_raster, maps to ggplot point size (~1.2 + 5*pt_size)."),
  
  # Flipped Toggle Flags
  make_option("--png", action="store_true", default=FALSE, help="Save as PNG instead of PDF (PDF is default)"),
  make_option("--no_raster", action="store_true", default=FALSE, help="Disable rasterizing points (Raster is default)"),
  make_option("--show_ref", action="store_true", default=FALSE, help="Show the reference UMAP plot (Hidden by default)"),
  make_option("--no_labels", action="store_true", default=FALSE, help="Hide active cell type labels (Shown by default)"),
  make_option("--beta_umap_all_cells", action="store_true", default=FALSE,
              help="Colour beta on every cell type in each Beta UMAP (legacy). Default: only the focal cell type for that credible set is on-scale; others are grey.")
)

opt <- parse_args(OptionParser(option_list = option_list))
# --beta_files / --beta_axes are the current names; --z_files / --z_axes are
# kept as deprecated aliases. Prefer the beta-named flag, then fall back to the
# z-named one, then to the default. Downstream code keeps using opt$z_files /
# opt$z_axes as the internal names.
if (!is.null(opt$beta_files) && !(length(opt$beta_files) == 1L && is.na(opt$beta_files[[1L]]))) {
  opt$z_files <- opt$beta_files
}
beta_axes_set <- !is.null(opt$beta_axes) &&
  !(length(opt$beta_axes) == 1L && is.na(opt$beta_axes[[1L]]))
if (beta_axes_set) {
  opt$z_axes <- opt$beta_axes
} else if (is.null(opt$z_axes) || (length(opt$z_axes) == 1L && is.na(opt$z_axes[[1L]]))) {
  opt$z_axes <- "linked"
}
# optparse may leave unused character flags as NA (length-1) rather than NULL;
# treat those as "not provided" so LZ-only runs never hit fread(opt$umap).
opt$umap    <- if (is.null(opt$umap) || length(opt$umap) == 0L ||
                  (length(opt$umap) == 1L && is.na(opt$umap[[1L]]))) NULL else opt$umap
opt$master  <- if (is.null(opt$master) || length(opt$master) == 0L ||
                  (length(opt$master) == 1L && is.na(opt$master[[1L]]))) NULL else opt$master
opt$modules <- if (is.null(opt$modules) || length(opt$modules) == 0L ||
                  (length(opt$modules) == 1L && is.na(opt$modules[[1L]]))) NULL else opt$modules
opt$lz_files <- if (is.null(opt$lz_files) || length(opt$lz_files) == 0L ||
                   (length(opt$lz_files) == 1L && is.na(opt$lz_files[[1L]]))) NULL else opt$lz_files
opt$z_files <- if (is.null(opt$z_files) || length(opt$z_files) == 0L ||
                  (length(opt$z_files) == 1L && is.na(opt$z_files[[1L]]))) NULL else opt$z_files
opt$ld_files <- if (is.null(opt$ld_files) || length(opt$ld_files) == 0L ||
                   (length(opt$ld_files) == 1L && is.na(opt$ld_files[[1L]]))) NULL else opt$ld_files
# getopt can assign the integer default of --umap_sample_n into --umap when
# both flags share a prefix; drop any non-character so fread() is never called
# with -1.
if (!is.null(opt$umap) && !is.character(opt$umap)) opt$umap <- NULL
anno_col <- if (!is.null(opt$anno_join_col)) opt$anno_join_col else opt$join_col

save_pdf   <- !opt$png
use_raster <- !opt$no_raster
show_ref   <- opt$show_ref
show_labels <- !opt$no_labels

# optparse/getopt can mis-bind values when one long flag is a prefix of another
# (e.g. --umap vs --umap_sample_n). Treat paths only as "set" when they are
# non-empty character strings (not integers like the default -1).
has_nonempty_char <- function(x) {
  !is.null(x) && is.character(x) && length(x) >= 1L &&
    !is.na(x[[1L]]) && nzchar(x[[1L]])
}

# --- Validation & Parsing ---
has_umap    <- has_nonempty_char(opt$umap)
has_master  <- has_nonempty_char(opt$master)
has_modules <- has_nonempty_char(opt$modules)
has_z       <- has_nonempty_char(opt$z_files)
has_lz      <- has_nonempty_char(opt$lz_files)
has_ld      <- has_nonempty_char(opt$ld_files)

# --- MANUAL MODE: --cs_name + --beta_cs ---
# A single credible-set identifier (--cs_name) plus a single beta (--beta_cs)
# lets you plot one (cell, gene, snp) without supplying --master or --modules.
# Both flags must be set together. When active, parse_cs(--cs_name) provides
# the focal cell, gene, chromosome and lead SNP; --beta_cs is the beta used
# to colour the focal cell type on the Beta UMAP. --cell, if set, overrides
# the parsed cell.
cs_name_set <- has_nonempty_char(opt$cs_name)
beta_cs_set <- !is.null(opt$beta_cs) && length(opt$beta_cs) == 1L &&
               is.finite(suppressWarnings(as.numeric(opt$beta_cs)))
if (cs_name_set != beta_cs_set)
  stop("--cs_name and --beta_cs must be provided together (or neither).")
manual_mode <- cs_name_set && beta_cs_set
parsed_cs <- if (manual_mode) parse_cs(opt$cs_name) else NULL
manual_label <- NA_character_
if (manual_mode) {
  toks     <- strsplit(opt$cs_name, "::", fixed = TRUE)[[1]]
  last_tok <- if (length(toks) > 0L) tail(toks, 1L) else ""
  manual_label <- if (length(last_tok) == 1L && grepl("^L\\d+$", last_tok)) last_tok
                  else if (!is.na(parsed_cs$chrom) && !is.na(parsed_cs$lead_pos))
                    sprintf("chr%s:%d", parsed_cs$chrom, as.integer(parsed_cs$lead_pos))
                  else "manual_cs"
  if (has_master)
    cat("Note: --master is ignored in --cs_name manual mode.\n")
  if (has_modules)
    cat("Note: --modules is ignored in --cs_name manual mode.\n")
  opt$master  <- NULL
  opt$modules <- manual_label   # feeds the existing has_modules branch below
  has_master  <- FALSE
  has_modules <- TRUE
}

umap_color_mode <- tolower(trimws(as.character(opt$umap_color_mode)))
if (!umap_color_mode %in% c("beta", "expression"))
  stop("--umap_color_mode must be 'beta' or 'expression'.")

# Beta UMAP needs --umap + (--master or manual --cs_name/--beta_cs).
# Expression UMAP needs --umap only; per-cell values come from a gene column
# in the UMAP file (matched to each CS row's eGene).
has_expr_panel  <- has_umap && umap_color_mode == "expression"
has_beta_panel  <- has_umap && umap_color_mode == "beta" && (has_master || manual_mode)
has_umap_panel  <- has_expr_panel || has_beta_panel

if (has_umap && umap_color_mode == "beta" && !has_master && !manual_mode)
  cat("Note: --umap provided without --master; the Beta UMAP panel is skipped (master annotations hold the per-(cell, gene) betas).\n")
if (has_umap && umap_color_mode == "expression")
  cat("Note: --umap_color_mode expression; per-cell expression is read from the --umap file (not --master).\n")
if (has_master && !has_umap)
  cat("Note: --master provided without --umap; the UMAP panel is skipped (no embedding to plot on).\n")

lz_layout_val <- tolower(trimws(opt$lz_layout))
if (!lz_layout_val %in% c("stacked", "merged"))
  stop("--lz_layout must be 'stacked' or 'merged'.")
if (lz_layout_val == "merged" && (has_umap_panel || has_z)) {
  cat("Warning: --lz_layout merged applies only to LocusZoom-only runs (no UMAP panel, no --beta_files); using stacked.\n")
  lz_layout_val <- "stacked"
}
use_lz_merged <- has_lz && lz_layout_val == "merged" && !has_umap_panel && !has_z

lz_xlim_mode <- tolower(trimws(as.character(opt$lz_xlim)))
if (!lz_xlim_mode %in% c("context", "snp"))
  stop("--lz_xlim must be 'context' or 'snp'.")

z_axes_mode <- tolower(trimws(as.character(opt$z_axes)))
if (!z_axes_mode %in% c("linked", "independent"))
  stop("--z_axes must be 'linked' or 'independent'.")
z_linked_axes <- z_axes_mode == "linked"

# Need at least something to plot.
if (!has_lz && !has_umap_panel)
  stop("Nothing to plot: provide --lz_files (LocusZoom only), or --umap + --master (with --modules, beta mode), or --umap + --cs_name + --beta_cs (manual beta mode), or --umap + --lz_files/--master (expression mode).")

beta_umap_all_cells <- isTRUE(opt$beta_umap_all_cells)

if (show_ref && is.null(opt$colors))
  stop("You must provide the --colors file if you enable --show_ref.")
if (show_ref && !has_umap)
  stop("--show_ref requires --umap (there is no UMAP to draw the reference on).")

cell_filter <- {
  raw <- NULL
  if (!is.null(opt[["cell"]]) && nzchar(as.character(opt[["cell"]]))) {
    raw <- as.character(opt[["cell"]])
  } else if (!is.null(opt[["cells"]]) && nzchar(as.character(opt[["cells"]]))) {
    raw <- as.character(opt[["cells"]])
  }
  if (!is.null(raw)) trimws(unlist(strsplit(raw, ","))) else NULL
}
gene_filter <- if (!is.null(opt[["gene"]])) trimws(unlist(strsplit(opt[["gene"]], ","))) else NULL

meta <- NULL
if (has_master) meta <- fread(opt$master)

# Establish the module list (either explicit via --modules or synthesised
# from --lz_files, one implicit module per file). With --umap + --master and
# no --lz_files, all modules in the master table are used when --modules is omitted.
if (has_modules) {
  mods <- trimws(unlist(strsplit(opt$modules, ",")))
} else if (has_lz) {
  tmp_lz <- trimws(unlist(strsplit(opt$lz_files, ",")))
  mods <- sprintf("auto_%s",
                  tools::file_path_sans_ext(basename(tmp_lz)))
  cat(sprintf("--modules not provided; treating each of the %d LZ file(s) as one implicit module: %s\n",
              length(mods), paste(mods, collapse = ", ")))
} else if (has_beta_panel && !is.null(meta) && "module" %in% names(meta)) {
  mods <- sort(unique(as.character(stats::na.omit(meta$module))))
  if (length(mods) == 0L)
    stop("No modules found in --master (column 'module') for --modules omission.")
  cat(sprintf("--modules not provided; using all %d module(s) from --master.\n", length(mods)))
} else {
  stop("Provide --modules, --lz_files, or (--umap + --master) so module IDs can be determined.")
}
n_items <- length(mods)

if (has_z) {
  z_files <- trimws(unlist(strsplit(opt$z_files, ",")))
  if (length(z_files) != n_items) stop("--beta_files must have the same number of items as --modules (or --lz_files when --modules is omitted).")
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

umap <- NULL; merged <- NULL
if (has_umap) {
  umap <- fread(opt$umap)
  if ("UMAP_0" %in% names(umap) && "UMAP_1" %in% names(umap))
    umap <- umap %>% rename(UMAP_2 = UMAP_1, UMAP_1 = UMAP_0)
  merged <- if (show_ref)
              umap %>% left_join(colors_df %>% select(all_of(opt$join_col),
                                                       hex_color = color_ct2),
                                 by = opt$join_col)
            else umap
  n_umap_in <- nrow(merged)
  if (n_umap_in > 0) {
    cap <- NULL
    if (!is.na(opt$umap_sample_n) && opt$umap_sample_n >= 1L) {
      cap <- min(as.integer(opt$umap_sample_n), n_umap_in)
    } else if (n_umap_in > opt$max_cells) {
      cap <- as.integer(opt$max_cells)
    }
    if (!is.null(cap) && n_umap_in > cap) {
      set.seed(42L)
      merged <- merged %>% slice_sample(n = cap)
      cat(sprintf("UMAP subsampled to %d cells (from %d).\n", cap, n_umap_in))
    }
  }
}
# meta already loaded when has_master (see module list above).
module_composites <- list()  # one patchwork per module (CS grid + merged box)
total_cs_rows     <- 0L
# Column layout: LZ, optional UMAP (beta or expression), optional Coloc Beta.
cols <- as.integer(has_lz) + as.integer(has_umap_panel) + as.integer(has_z)
if (cols == 0)
  stop("Nothing to plot (no LZ, no UMAP panel, no Z).")

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
  # produces one grid row of [LocusZoom | Beta UMAP | optional Coloc Beta].
  lz_path  <- if (has_lz) lz_files[i] else NULL
  if (manual_mode) {
    # In manual mode we don't have a master annotation to enumerate
    # (cell, gene) pairs from, so synthesise exactly one row from
    # parse_cs(--cs_name). --cell (when provided) wins over the parsed cell.
    pc <- parsed_cs
    focal_cell_man <- if (!is.null(cell_filter) && length(cell_filter) > 0L)
                       cell_filter[1] else pc$cell
    cs_rows <- data.table::data.table(
      cs           = opt$cs_name,
      cell         = focal_cell_man,
      eGene        = pc$eGene,
      eGene_symbol = pc$eGene,  # no symbol lookup without master
      chrom        = pc$chrom,
      lead_pos     = pc$lead_pos
    )
  } else {
    cs_rows  <- module_cs_list(lz_path, meta, mod_id, anno_col,
                               cell_filter = cell_filter,
                               gene_filter = gene_filter)
  }
  if (nrow(cs_rows) == 0) {
    cat(sprintf("Warning: Module %s has no credible sets matching the filters - skipping.\n", mod_id))
    next
  }
  
  # LD matrix (r^2) for this module; skipped when LZ panels are merged (LD is
  # disabled there; points are coloured by credible set instead).
  ld_mat <- if (has_ld && !use_lz_merged) load_ld_matrix(ld_files[i]) else NULL
  
  z_tbl_mod <- NULL
  if (has_z) {
    zpath <- z_files[i]
    if (!is.na(zpath) && zpath != "NA" && file.exists(zpath)) {
      z_tbl_mod <- fread(zpath)
    } else if (!is.na(zpath) && zpath != "NA") {
      cat(sprintf("Warning: Z-score file not found at %s\n", zpath))
    }
  }
  
  module_plots   <- list()
  cs_row_titles  <- list()     # one title+subtitle per CS row (patchwork)
  lead_bag       <- numeric(0) # lead SNP positions, for the merged box
  cs_labels      <- character(0)
  chr_for_mod    <- NA_character_
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
  if (manual_mode) {
    # Anchor SNP comes from the parsed --cs_name. Try to use the full
    # chr:pos:ref:alt token (better LD lookup); fall back to chr:pos:N:N.
    mod_master_snp_pos <- parsed_cs$lead_pos
    pc_tok <- regmatches(opt$cs_name,
                         regexpr("chr[0-9XYM]+:\\d+:[ACGTN]+:[ACGTN]+",
                                 opt$cs_name, ignore.case = TRUE))
    if (length(pc_tok) == 1L && nzchar(pc_tok)) {
      mod_master_snp_lab <- pc_tok
      mod_master_snp_id  <- pc_tok
    } else if (!is.na(parsed_cs$chrom) && !is.na(parsed_cs$lead_pos)) {
      mod_master_snp_lab <- sprintf("chr%s:%d", parsed_cs$chrom,
                                    as.integer(parsed_cs$lead_pos))
      mod_master_snp_id  <- sprintf("chr%s:%d:N:N", parsed_cs$chrom,
                                    as.integer(parsed_cs$lead_pos))
    }
  } else if (has_master) {
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
  
  if (use_lz_merged && has_lz) {
    cat(sprintf("  [%s] LocusZoom layout: merged (all credible sets in one panel).\n", mod_id))
    p1 <- 1L
    this_cell <- cs_rows$cell[p1]
    this_gene <- cs_rows$eGene[p1]
    this_sym  <- cs_rows$eGene_symbol[p1]
    if (is.na(this_sym) || this_sym == "") this_sym <- this_gene
    
    disp_snp <- if (!is.na(mod_master_snp_lab) && nzchar(mod_master_snp_lab))
                  mod_master_snp_lab else NA_character_
    this_cs0 <- cs_rows$cs[p1]
    if (is.na(disp_snp) && !is.na(this_cs0)) {
      m <- regmatches(this_cs0,
                      regexpr("chr[0-9XYM]+:\\d+:[ACGTN]+:[ACGTN]+",
                              this_cs0, ignore.case = TRUE))
      if (length(m) == 1 && nzchar(m)) disp_snp <- m
    }
    if (is.na(disp_snp) || disp_snp == "") disp_snp <- "Unknown_SNP"
    beta_disp_m <- if (manual_mode) suppressWarnings(as.numeric(opt$beta_cs))
                   else if (valid_b) b_dis else NA_real_
    cs_row_titles[[length(cs_row_titles) + 1L]] <-
      build_cs_figure_title(mod_id, this_cell, this_sym, disp_snp, beta_disp_m)
    lz_title <- "LocusZoom (merged credible sets)"
    
    locus_info <- if (manual_mode) {
                    list(chrom       = parsed_cs$chrom,
                         gene_start  = NA_real_,
                         gene_end    = NA_real_,
                         gene_strand = NA_character_,
                         gene_symbol = parsed_cs$eGene,
                         lead_pos    = parsed_cs$lead_pos)
                  } else if (has_master) {
                    get_module_locus_info(meta, mod_id,
                                          focal_cell = this_cell,
                                          focal_gene = this_gene,
                                          anno_col   = anno_col,
                                          gene_col   = opt$gene_col)
                  } else NULL
    if (is.null(locus_info)) locus_info <- list()
    if (!is.na(cs_rows$lead_pos[p1])) {
      locus_info$lead_pos <- cs_rows$lead_pos[p1]
      if (is.null(locus_info$chrom) || is.na(locus_info$chrom))
        locus_info$chrom <- cs_rows$chrom[p1]
    }
    if (is.null(locus_info$chrom) || is.na(locus_info$chrom)) {
      lz_peek <- load_lz_file(lz_files[i])
      if (!is.null(lz_peek) && nrow(lz_peek) > 0)
        locus_info$chrom <- sub("^chr", "", as.character(lz_peek$CHR[1]),
                                ignore.case = TRUE)
    }
    
    lz_full <- load_lz_file(lz_files[i])
    genes_region <- NULL
    if (!is.null(gtf_tbl) && !is.null(locus_info) &&
        !is.null(locus_info$chrom) && !is.na(locus_info$chrom) &&
        !is.null(lz_full) && nrow(lz_full) > 0) {
      win_s <- min(lz_full$POS, na.rm = TRUE)
      win_e <- max(lz_full$POS, na.rm = TRUE)
      if (!is.null(locus_info$gene_start) && !is.na(locus_info$gene_start))
        win_s <- min(win_s, as.numeric(locus_info$gene_start))
      if (!is.null(locus_info$gene_end) && !is.na(locus_info$gene_end))
        win_e <- max(win_e, as.numeric(locus_info$gene_end))
      genes_region <- get_genes_in_region(gtf_tbl, locus_info$chrom,
                                          win_s, win_e)
    }
    
    for (p_idx in seq_len(nrow(cs_rows))) {
      this_cs   <- cs_rows$cs[p_idx]
      this_cell_i <- cs_rows$cell[p_idx]
      this_gene_i <- cs_rows$eGene[p_idx]
      this_sym_i  <- cs_rows$eGene_symbol[p_idx]
      if (is.na(this_sym_i) || this_sym_i == "") this_sym_i <- this_gene_i
      if (!is.na(cs_rows$lead_pos[p_idx])) {
        lead_bag <- c(lead_bag, cs_rows$lead_pos[p_idx])
        l_tok <- if (!is.na(this_cs)) sub(".*::([^:]+)$", "\\1", this_cs) else NA_character_
        short_sym <- if (!is.na(this_sym_i) && nchar(this_sym_i) > 0) this_sym_i else this_gene_i
        lab <- paste(c(this_cell_i, short_sym,
                       if (!is.na(l_tok) && nchar(l_tok) > 0 && l_tok != this_cs) l_tok else NULL),
                     collapse = " | ")
        cs_labels <- c(cs_labels, lab)
      }
      if (is.na(chr_for_mod) && !is.na(cs_rows$chrom[p_idx]) &&
          nzchar(as.character(cs_rows$chrom[p_idx])))
        chr_for_mod <- sub("^chr", "", as.character(cs_rows$chrom[p_idx]),
                           ignore.case = TRUE)
    }
    if (is.na(chr_for_mod) && !is.null(locus_info) &&
        !is.null(locus_info$chrom) && !is.na(locus_info$chrom))
      chr_for_mod <- sub("^chr", "", locus_info$chrom, ignore.case = TRUE)
    
    if (!is.null(lz_full) && !("CS" %in% names(lz_full)))
      cat(sprintf("Warning: --lz_layout merged expects a CS column in %s; plotting all rows together.\n",
                  lz_files[i]))
    
    module_plots[[length(module_plots) + 1]] <- plot_locuszoom(
      lz_files[i], locus_info,
      title                  = lz_title,
      genes_df               = genes_region,
      annotations_tbl        = NULL,
      annotation_tracks      = NULL,
      zoom_window            = opt$zoom_window,
      lz_filter_cell         = NULL,
      lz_filter_gene         = NULL,
      lz_filter_cs           = NULL,
      include_zoom_panel     = FALSE,
      include_zoom_connector = FALSE,
      ld_vec                 = NULL,
      force_lead_pos         = mod_master_snp_pos,
      merge_all_cs           = TRUE,
      xlim_mode              = lz_xlim_mode)
  } else {
  for (p_idx in seq_len(nrow(cs_rows))) {
    this_cs   <- cs_rows$cs[p_idx]
    this_cell <- cs_rows$cell[p_idx]
    this_gene <- cs_rows$eGene[p_idx]
    this_sym  <- cs_rows$eGene_symbol[p_idx]
    if (is.na(this_sym) || this_sym == "") this_sym <- this_gene
    
    # UMAP middle panel: beta (from --master) or expression (from --umap file).
    plot_df <- NULL
    expr_gene_col <- NULL
    if (has_expr_panel) {
      meta_cols <- expr_meta_columns(
        opt$join_col,
        if (has_nonempty_char(opt$label_col)) opt$label_col else NULL)
      tryCatch({
        expr_gene_col <- resolve_expr_gene_column(
          names(merged), this_gene, meta_cols, gene_symbol = this_sym)
        plot_df <- merged
      }, error = function(e) {
        cat(sprintf("Warning: Module %s, gene %s - %s\n",
                    mod_id, this_gene, conditionMessage(e)))
      })
    } else if (has_beta_panel) {
      beta_tbl <- if (manual_mode) {
        # One-row synthetic beta table: the focal cell (parsed from --cs_name
        # or overridden by --cell) gets --beta_cs, everyone else stays NA
        # (rendered grey by plot_beta() via na.value).
        tibble::tibble(
          !!opt$join_col := as.character(this_cell),
          beta_val      = as.numeric(opt$beta_cs),
          gene_sym_col  = as.character(this_sym),
          mod_id_col    = as.character(mod_id),
          snp_id        = NA_character_,
          pval_exact    = NA_real_
        )
      } else {
        get_module_betas(meta, mod_id, focal_gene = this_gene,
                         anno_col, opt$join_col, opt$gene_col)
      }
      if (is.null(beta_tbl)) {
        cat(sprintf("Warning: No beta rows for module %s, gene %s - skipping Beta UMAP for this CS.\n",
                    mod_id, this_gene))
      } else {
        jc_beta <- opt$join_col
        plot_df <- merged %>% left_join(beta_tbl, by = jc_beta) %>%
          mutate(beta_num = as.numeric(beta_val))
        if (!jc_beta %in% names(plot_df)) {
          stop("Beta UMAP: column ", jc_beta, " (--join_col) not found after joining master betas.")
        }
        # Default: one Beta UMAP per credible set colours only that CS's focal
        # cell type (others grey). --beta_umap_all_cells restores the old behaviour
        # (every cell type on the diverging scale). --cell still restricts which
        # CS rows exist; within a row, focal highlighting uses this_cell when known.
        focal_one_cs <- !beta_umap_all_cells && !is.na(this_cell) &&
          nzchar(as.character(this_cell))
        if (focal_one_cs) {
          focal <- as.character(plot_df[[jc_beta]]) == as.character(this_cell)
          plot_df <- plot_df %>%
            mutate(beta_val = if_else(focal, coalesce(beta_num, 0), NA_real_)) %>%
            # Draw non-focal (NA -> grey) cells FIRST so they sit in the
            # background; then focal cells in ascending |beta| so the most
            # extreme blues/reds end up on top and remain visible instead of
            # being painted over by the grey layer.
            arrange(desc(!focal), abs(coalesce(beta_val, 0))) %>%
            select(-beta_num)
        } else if (!is.null(cell_filter) && length(cell_filter) > 0) {
          focal <- as.character(plot_df[[jc_beta]]) %in% cell_filter
          plot_df <- plot_df %>%
            mutate(beta_val = if_else(focal, coalesce(beta_num, 0), NA_real_)) %>%
            arrange(desc(!focal), abs(coalesce(beta_val, 0))) %>%
            select(-beta_num)
        } else {
          plot_df <- plot_df %>%
            mutate(beta_val = coalesce(beta_num, 0)) %>%
            arrange(abs(beta_val)) %>%
            select(-beta_num)
        }
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

    # Beta for the figure subtitle (focal cell / manual --beta_cs).
    beta_disp <- NA_real_
    if (manual_mode) {
      beta_disp <- suppressWarnings(as.numeric(opt$beta_cs))
    } else if (!is.null(plot_df) && "beta_val" %in% names(plot_df)) {
      fb <- unique(stats::na.omit(plot_df$beta_val[!is.na(plot_df$beta_val)]))
      if (length(fb) >= 1L) beta_disp <- as.numeric(fb[1])
    }
    cs_row_titles[[length(cs_row_titles) + 1L]] <-
      build_cs_figure_title(mod_id, this_cell, this_sym, disp_snp, beta_disp,
                            color_mode = umap_color_mode)

    if (has_lz) {
      # Build locus info. If the master annotation is available, use it for
      # gene coords / eGene symbol / strand; otherwise start from an empty
      # list and rely on the CS-derived chrom + lead_pos alone.
      locus_info <- if (manual_mode) {
                      list(chrom       = parsed_cs$chrom,
                           gene_start  = NA_real_,
                           gene_end    = NA_real_,
                           gene_strand = NA_character_,
                           gene_symbol = parsed_cs$eGene,
                           lead_pos    = parsed_cs$lead_pos)
                    } else if (has_master) {
                      get_module_locus_info(meta, mod_id,
                                            focal_cell = this_cell,
                                            focal_gene = this_gene,
                                            anno_col   = anno_col,
                                            gene_col   = opt$gene_col)
                    } else NULL
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
      # Panel label only; module / cell / gene / SNP / beta live in the
      # row-level plot_annotation (Figure 5 panel-b style).
      lz_title <- "LocusZoom"
      
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
        force_lead_pos         = mod_master_snp_pos,
        xlim_mode              = lz_xlim_mode)
      
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
    label_cells_arg <- if (!is.na(this_cell) && nzchar(as.character(this_cell))) {
      as.character(this_cell)
    } else if (!is.null(cell_filter) && length(cell_filter) > 0) {
      cell_filter
    } else {
      NULL
    }
    label_col_arg <- if (has_nonempty_char(opt$label_col)) opt$label_col else opt$join_col

    if (has_expr_panel && !is.null(plot_df) && !is.null(expr_gene_col)) {
      pe_args <- list(
        plot_df,
        expr_gene_col,
        opt$pt_size,
        opt$join_col,
        show_legend = TRUE,
        use_raster = use_raster,
        show_labels = show_labels,
        label_cells_only = label_cells_arg,
        label_col = label_col_arg,
        log1p = isTRUE(opt$log1p),
        clip_q = opt$clip_quantile,
        palette = opt$expr_palette
      )
      module_plots[[length(module_plots) + 1]] <- do.call(plot_expression, pe_args)
    } else if (has_expr_panel) {
      module_plots[[length(module_plots) + 1]] <- plot_spacer()
    } else if (has_beta_panel && !is.null(plot_df)) {
      pb_args <- list(
        plot_df,
        NULL,
        opt$pt_size,
        opt$join_col,
        show_legend = TRUE,
        use_raster = use_raster,
        show_labels = show_labels
      )
      if ("label_col" %in% names(formals(plot_beta))) {
        pb_args$label_col <- label_col_arg
      }
      if ("label_cells_only" %in% names(formals(plot_beta))) {
        pb_args$label_cells_only <- label_cells_arg
      }
      module_plots[[length(module_plots) + 1]] <- do.call(plot_beta, pb_args)
    } else if (has_beta_panel) {
      module_plots[[length(module_plots) + 1]] <- plot_spacer()
    }
    if (has_z) {
      module_plots[[length(module_plots) + 1]] <- build_beta_column(
        z_tbl_mod, "Coloc Betas", this_cs = this_cs,
        linked_axes = z_linked_axes)
    }
  }
  }
  
  if (length(module_plots) == 0) next
  
  n_cs_for_mod <- length(module_plots) %/% cols
  total_cs_rows <- total_cs_rows + n_cs_for_mod
  # Per-CS-row column widths: LocusZoom (when present) is rendered narrower
  # than the other side panels (beta UMAP, Z-score) but with enough room for
  # its multi-line title and gene track to lay out cleanly. The figure-wide
  # plot_width below scales by sum(col_widths) so the other panels actually
  # close the gap instead of leaving empty horizontal space (beta UMAP + coloc-beta
  # use coord_equal()/coord_fixed() and can't grow into the slack).
  col_widths <- rep(1, cols)
  if (has_lz) col_widths[1] <- 1.2
  cs_grid <- wrap_cs_grid_with_titles(module_plots, cols, col_widths, cs_row_titles)
  # Apply the global Helvetica / large / plain text styling NOW, before
  # wrap_elements() makes cs_grid atomic. `& theme()` does NOT propagate
  # through wrap_elements(), so the final `final_plot & big_helvetica_theme()`
  # at the bottom never reaches the LZ / beta / Z panels -- only the merged
  # annotation box and (when shown) the reference UMAP.
  cs_grid <- cs_grid & big_helvetica_theme(base_size = 18)
  
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
      extra_lead_label  = {
        lab <- NA_character_
        if (!is.na(mod_master_snp_lab) && nzchar(mod_master_snp_lab) &&
            !startsWith(mod_master_snp_lab, "most-likely SNP |"))
          lab <- as.character(mod_master_snp_lab)
        else if (!is.na(mod_master_snp_id) && nzchar(as.character(mod_master_snp_id)))
          lab <- as.character(mod_master_snp_id)
        lab
      })
    # The merged annotation box was full-width before, which made it look
    # disconnected from the LocusZoom even though they share the same
    # genomic axis. Now we put it in a 1-row patchwork beside `plot_spacer()`
    # blocks for the beta / Z columns and reuse the SAME col_widths, so the
    # box visually anchors directly under the LocusZoom column.
    #
    # It also uses a slightly smaller base text size than the top panels
    # since its column is narrower (text wraps better at ~14 pt).
    merged_box <- merged_box + big_helvetica_theme(base_size = 14)
    if (cols > 1L) {
      merged_row_parts <- c(list(merged_box),
                            replicate(cols - 1L, patchwork::plot_spacer(),
                                      simplify = FALSE))
      merged_row <- patchwork::wrap_plots(merged_row_parts, ncol = cols) +
                    patchwork::plot_layout(widths = col_widths)
    } else {
      merged_row <- merged_box
    }
    # Heights: each CS row ~5 units, merged box ~1.5 units (lowered from 2
    # so the top row keeps more vertical real-estate, since its panels are
    # the larger / more information-dense ones). Wrapping `cs_grid` AND the
    # merged row with `wrap_elements()` keeps each block atomic so `/` stacks
    # them cleanly without flattening inner panels.
    # merged_h bumped slightly (1.5 -> 2) so the y-axis track labels
    # (promoter_flanking_region, open_chromatin_region, ...) don't end up
    # stacked on top of each other now that the merged box is also
    # narrower (LZ-column width only).
    merged_h <- 2
    composite <- wrap_elements(cs_grid) / wrap_elements(merged_row) +
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

# Effective per-row width is the sum of col_widths (the LZ slot is 0.5,
# others are 1) -- not `cols` -- so when the LZ is halved the row shrinks
# accordingly instead of leaving empty horizontal space.
row_w <- sum(col_widths)
if (show_ref) {
  ref_label_col <- if (has_nonempty_char(opt$label_col)) opt$label_col else opt$join_col
  pr_args <- list(merged, paste(opt$name, "Reference"), opt$join_col,
                  opt$pt_size, use_raster)
  if ("label_col" %in% names(formals(plot_ref))) pr_args$label_col <- ref_label_col
  p_ref <- do.call(plot_ref, pr_args)
  final_plot <- p_ref | grid_plot
  final_plot <- final_plot + plot_layout(widths = c(1, row_w))
  plot_width <- (row_w + 1) * 6
} else {
  final_plot <- grid_plot
  plot_width <- row_w * 6
}

# Bumped from 4.5 to 6 so the top row's panels (LZ, beta UMAP, Z-score) keep
# enough vertical room after the new bigger Helvetica titles + bigger axis text.
base_height <- 6
plot_height <- max(base_height * n_rows, 6)

# Helvetica / large / plain text on the outer patchwork. This reaches the
# merged annotation box and (when present) the reference UMAP. The inner
# LZ / beta / Z grid was already styled above before wrap_elements().
final_plot <- final_plot & big_helvetica_theme(base_size = 18)

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
