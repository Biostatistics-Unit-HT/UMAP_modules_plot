#!/usr/bin/env Rscript
# Fast UMAP plots coloured by betas (PDF, Raster, Labels, Layout Toggles)
# Dependencies: ggplot2, dplyr, data.table, optparse, patchwork, ggrepel, ragg, ggrastr

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
  make_option("--umap", type = "character", help = "Path to UMAP TSV"),
  make_option("--colors", type = "character", help = "Path to palette TSV"),
  make_option("--annotations", type = "character", help = "Path to Annotations TSV (containing modules/betas)"),
  
  make_option("--module_a", type = "character", help = "Primary Module ID (e.g., M_13916)"),
  make_option("--ensg_a", type = "character", help = "Primary eGene ID (e.g., ENSG0000...)"),
  
  make_option("--module_b", type = "character", default = NULL, help = "Secondary Module ID (optional)"),
  make_option("--ensg_b", type = "character", default = NULL, help = "Secondary eGene ID (optional)"),
  
  make_option("--join_col", type = "character", default = "celltype_2", help = "Shared column for mapping across all tables [default: %default]"),
  make_option("--out", type = "character", default = "umap_module_betas.png", help = "Output filename [default: %default]"),
  make_option("--pt_size", type = "numeric", default = 0.25, help = "Point size [default: %default]"),
  make_option("--max_cells", type = "numeric", default = 250000, help = "Max cells to plot to speed up rendering [default: %default]"),
  
  # New Flags
  make_option("--pdf", action="store_true", default=FALSE, help="Save as PDF instead of PNG"),
  make_option("--raster", action="store_true", default=FALSE, help="Rasterize points (Highly recommended for PDF to save file size)"),
  make_option("--no_ref", action="store_true", default=FALSE, help="Hide the reference UMAP plot on the left"),
  make_option("--label_betas", action="store_true", default=FALSE, help="Overlay cell type labels onto the Beta UMAPs")
)

opt <- parse_args(OptionParser(option_list = option_list))

# --- Validation ---
required_args <- c("umap", "colors", "annotations", "module_a", "ensg_a")
missing_args <- setdiff(required_args, names(opt))
if (length(missing_args) > 0) stop(paste("Missing required arguments:", paste(missing_args, collapse = ", ")))

# --- Helper Functions ---

# Plot reference UMAP
plot_ref <- function(df, join_col, pt_size, use_raster) {
  centroids <- df %>%
    group_by(.data[[join_col]]) %>%
    summarise(UMAP_1 = median(UMAP_1, na.rm = TRUE), UMAP_2 = median(UMAP_2, na.rm = TRUE), .groups = "drop")
  
  p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2))
  
  # Apply rasterization if requested
  if (use_raster) {
    # pixels sets the resolution of the seamless image. 2048 is ultra-crisp.
    # scattermore point sizes scale slightly differently, so we usually add +1 or +2
    p <- p + geom_scattermore(aes(color = hex_color), pointsize = pt_size + 1.5, pixels = c(2048, 2048))
  } else {
    p <- p + geom_point(aes(color = hex_color), size = pt_size, stroke = 0)
  }
  
  p + scale_color_identity() +
    geom_text_repel(data = centroids, aes(label = .data[[join_col]]), size = 2.5, fontface = "bold", bg.color = "white", bg.r = 0.1) +
    coord_equal() +
    labs(title = "UMAP Reference", x = NULL, y = NULL) +
    theme_minimal() +
    theme(panel.grid = element_blank(), axis.text = element_blank())
}

# Plot Beta UMAP
plot_beta <- function(df, title, pt_size, join_col, show_legend = TRUE, use_raster = FALSE, show_labels = FALSE) {
  b_max <- max(abs(df$beta_val), na.rm = TRUE)
  if (!is.finite(b_max) || b_max <= 0) b_max <- 1e-9
  
  # Swapped the pale yellow (#FFF9C4) for a vibrant, aggressive yellow (#FFEB3B)
  # and the amber (#FFC107) for a deeper orange (#FF9800)
  colors_beta <- c("#00008B", "#1E90FF", "#90CAF9", "#E8E8E8", "#FFEB3B", "#FF9800", "#E53935")
  
  p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = beta_val))
  
  if (use_raster) {
    p <- p + geom_scattermore(pointsize = pt_size + 1.5, pixels = c(2048, 2048))
  } else {
    p <- p + geom_point(size = pt_size, stroke = 0)
  }
  
  p <- p + scale_color_gradientn(
    colors = colors_beta, 
    limits = c(-b_max, b_max), 
    na.value = "grey50", 
    name = "Beta",
    # SQUEEZE THE WHITE: 
    # 0.5 is exactly Zero (White).
    # Setting the next color to 0.55 means yellow kicks in almost immediately!
    values = c(0.0, 0.20, 0.45, 0.5, 0.55, 0.80, 1.0)
  ) +    coord_equal() +
    labs(title = title, x = NULL, y = NULL) +
    theme_minimal() +
    theme(panel.grid = element_blank(), axis.text = element_blank(), plot.title = element_text(size = 9))
  
  if (!show_legend) p <- p + theme(legend.position = "none")
  
  # Add cell type labels to the beta plot if requested
  if (show_labels) {
    centroids <- df %>%
      group_by(.data[[join_col]]) %>%
      summarise(UMAP_1 = median(UMAP_1, na.rm = TRUE), UMAP_2 = median(UMAP_2, na.rm = TRUE), .groups = "drop")
    
    p <- p + geom_text_repel(
      data = centroids, 
      aes(x = UMAP_1, y = UMAP_2, label = .data[[join_col]]), 
      inherit.aes = FALSE, # Bypasses the color=beta_val requirement
      size = 2.5, fontface = "bold", bg.color = "white", bg.r = 0.1, color = "black"
    )
  }
  
  return(p)
}

# Filter and deduplicate beta values
get_module_betas <- function(df, mod_id, ensg_id, join_col) {
  sub_df <- df %>% filter(module == mod_id, eGene == ensg_id)
  if (nrow(sub_df) == 0) stop(paste("No rows found in annotations for module:", mod_id, "and eGene:", ensg_id))
  if ("cs_max_pip" %in% names(sub_df)) sub_df <- sub_df %>% arrange(desc(as.numeric(cs_max_pip)))
  
  sub_df %>% 
    distinct(.data[[join_col]], .keep_all = TRUE) %>%
    mutate(final_beta = if("most_likely_snp_beta" %in% names(.)) most_likely_snp_beta else beta) %>%
    select(all_of(join_col), beta_val = final_beta)
}

# --- Main Logic ---

cat("Loading data...\n")
umap_df <- fread(opt$umap)
colors_df <- fread(opt$colors)
meta_df <- fread(opt$annotations)

if ("UMAP_0" %in% names(umap_df) && "UMAP_1" %in% names(umap_df)) {
  umap_df <- umap_df %>% rename(UMAP_2 = UMAP_1, UMAP_1 = UMAP_0)
}

merged_base <- umap_df %>%
  left_join(colors_df %>% select(all_of(opt$join_col), hex_color = color_ct2), by = opt$join_col)

if (nrow(merged_base) > opt$max_cells) {
  cat(sprintf("Dataset is large. Subsampling to %s cells...\n", format(opt$max_cells, big.mark=",")))
  set.seed(42)
  merged_base <- merged_base %>% slice_sample(n = opt$max_cells)
}

# Generate Reference Plot (if not disabled)
if (!opt$no_ref) {
  p_ref <- plot_ref(merged_base, opt$join_col, opt$pt_size, opt$raster)
}

# Generate Beta A Plot
cat(paste("Processing Primary Module:", opt$module_a, "\n"))
beta_a <- get_module_betas(meta_df, opt$module_a, opt$ensg_a, opt$join_col)
plot_df_a <- merged_base %>% left_join(beta_a, by = opt$join_col) %>% mutate(beta_val = coalesce(as.numeric(beta_val), 0)) %>% arrange(abs(beta_val))
p_a <- plot_beta(plot_df_a, paste(opt$module_a, "|", opt$ensg_a), opt$pt_size, opt$join_col, show_legend = TRUE, use_raster = opt$raster, show_labels = opt$label_betas)

# Generate Beta B Plot (If requested)
has_b <- !is.null(opt$module_b) && !is.null(opt$ensg_b)
if (has_b) {
  cat(paste("Processing Secondary Module:", opt$module_b, "\n"))
  beta_b <- get_module_betas(meta_df, opt$module_b, opt$ensg_b, opt$join_col)
  plot_df_b <- merged_base %>% left_join(beta_b, by = opt$join_col) %>% mutate(beta_val = coalesce(as.numeric(beta_val), 0)) %>% arrange(abs(beta_val))
  p_b <- plot_beta(plot_df_b, paste(opt$module_b, "|", opt$ensg_b), opt$pt_size, opt$join_col, show_legend = FALSE, use_raster = opt$raster, show_labels = opt$label_betas)
}

# Layout Logic
if (opt$no_ref) {
  final_plot <- if (has_b) (p_a / p_b) else p_a
  plot_width <- 5
  plot_height <- if (has_b) 8 else 5
} else {
  final_plot <- if (has_b) p_ref | (p_a / p_b) else p_ref | p_a
  plot_width <- 11
  plot_height <- if (has_b) 8 else 5
}

# --- Save Logic ---
if (opt$pdf) {
  # Force .pdf extension if saving as PDF
  out_file <- sub("\\.png$", ".pdf", opt$out, ignore.case = TRUE)
  if (!grepl("\\.pdf$", out_file, ignore.case = TRUE)) out_file <- paste0(out_file, ".pdf")
  
  cat(paste("Saving as vector PDF with rasterized points to:", out_file, "\n"))
  ggsave(out_file, plot = final_plot, width = plot_width, height = plot_height, device = cairo_pdf)
} else {
  cat(paste("Saving as fast PNG to:", opt$out, "\n"))
  ggsave(opt$out, plot = final_plot, width = plot_width, height = plot_height, dpi = 300, bg = "white", device = ragg::agg_png)
}