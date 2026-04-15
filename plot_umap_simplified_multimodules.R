#!/usr/bin/env Rscript
# Ultra-Streamlined Dynamic UMAP plots (Smart defaults, Minimal inputs)

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
  make_option("--colors", type = "character", default = NULL, help = "Path to palette TSV (Only required if --show_ref is used)"),
  
  # Population 1 (Required)
  make_option("--umap_1", type = "character", help = "Path to UMAP TSV for Pop 1"),
  make_option("--master_1", type = "character", help = "Path to Master Annotations TSV for Pop 1"),
  make_option("--modules_1", type = "character", help = "Comma-separated Module IDs for Pop 1"),
  make_option("--name_1", type = "character", default = "Pop1", help = "Display name for Pop 1"),
  
  # Population 2 (Optional)
  make_option("--umap_2", type = "character", default = NULL, help = "Path to UMAP TSV for Pop 2"),
  make_option("--master_2", type = "character", default = NULL, help = "Path to Master Annotations TSV for Pop 2"),
  make_option("--modules_2", type = "character", default = NULL, help = "Comma-separated Module IDs for Pop 2"),
  make_option("--name_2", type = "character", default = "Pop2", help = "Display name for Pop 2"),
  
  # General Settings
  make_option("--join_col", type = "character", default = "cell", help = "Column name in UMAP for celltypes [default: %default]"),
  make_option("--anno_join_col", type = "character", default = NULL, help = "Column name in Master Annotations (if different from join_col)"),
  make_option("--gene_col", type = "character", default = "eGene_symbol", help = "Column name for gene symbol [default: %default]"),
  
  make_option("--out", type = "character", default = "umap_plot", help = "Output filename prefix"),
  make_option("--pt_size", type = "numeric", default = 0.25, help = "Point size"),
  make_option("--max_cells", type = "numeric", default = 250000, help = "Max cells to plot"),
  
  # Flipped Toggle Flags (PDF, Raster, No Ref, Label Betas are TRUE by default)
  make_option("--png", action="store_true", default=FALSE, help="Save as PNG instead of PDF (PDF is default)"),
  make_option("--no_raster", action="store_true", default=FALSE, help="Disable rasterizing points (Raster is default)"),
  make_option("--show_ref", action="store_true", default=FALSE, help="Show the reference UMAP plot (Hidden by default)"),
  make_option("--no_labels", action="store_true", default=FALSE, help="Hide active cell type labels (Shown by default)")
)

opt <- parse_args(OptionParser(option_list = option_list))
anno_col <- if(!is.null(opt$anno_join_col)) opt$anno_join_col else opt$join_col
dual_mode <- !is.null(opt$master_2)

# Inverse Logical Variables for readability
save_pdf <- !opt$png
use_raster <- !opt$no_raster
show_ref <- opt$show_ref
show_labels <- !opt$no_labels

# --- Validation & Parsing ---
required_args <- c("umap_1", "master_1", "modules_1")
missing_args <- setdiff(required_args, names(opt))
if (length(missing_args) > 0) stop(paste("Missing required arguments:", paste(missing_args, collapse = ", ")))

if (show_ref && is.null(opt$colors)) stop("You must provide the --colors file if you enable --show_ref.")

mods_1 <- trimws(unlist(strsplit(opt$modules_1, ",")))

if (dual_mode) {
  if (is.null(opt$umap_2) || is.null(opt$modules_2)) stop("If using Pop 2, --umap_2 and --modules_2 must be provided.")
  mods_2 <- trimws(unlist(strsplit(opt$modules_2, ",")))
  if (length(mods_1) != length(mods_2)) stop("For side-by-side mode, Pop 1 and Pop 2 module lists must have the same number of items.")
}

# --- Helper Functions ---
plot_ref <- function(df, title_text, join_col, pt_size, use_raster) {
  centroids <- df %>% group_by(.data[[join_col]]) %>% summarise(UMAP_1 = median(UMAP_1, na.rm = TRUE), UMAP_2 = median(UMAP_2, na.rm = TRUE), .groups = "drop")
  p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2))
  if (use_raster) p <- p + geom_scattermore(aes(color = hex_color), pointsize = pt_size + 1.5, pixels = c(2048, 2048)) else p <- p + geom_point(aes(color = hex_color), size = pt_size, stroke = 0)
  p + scale_color_identity() + geom_text_repel(data = centroids, aes(label = .data[[join_col]]), size = 3.5, fontface = "bold", bg.color = "white", bg.r = 0.1, min.segment.length = 0) + coord_equal() + labs(title = title_text, x = NULL, y = NULL) + theme_minimal() + theme(panel.grid = element_blank(), axis.text = element_blank(), plot.title = element_text(face="bold"))
}

plot_beta <- function(df, title, pt_size, join_col, show_legend = TRUE, use_raster = FALSE, show_labels = FALSE) {
  b_max <- max(abs(df$beta_val), na.rm = TRUE)
  if (!is.finite(b_max) || b_max <= 0) b_max <- 1e-9
  colors_beta <- c("#00008B", "#1E90FF", "#90CAF9", "#E8E8E8", "#FFEB3B", "#FF9800", "#E53935")
  
  p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = beta_val))
  if (use_raster) p <- p + geom_scattermore(pointsize = pt_size + 1.5, pixels = c(2048, 2048)) else p <- p + geom_point(size = pt_size, stroke = 0)
  
  p <- p + scale_color_gradientn(colors = colors_beta, limits = c(-b_max, b_max), na.value = "grey50", name = "Beta", values = c(0.0, 0.20, 0.45, 0.5, 0.55, 0.80, 1.0)) + 
    coord_equal() + labs(title = title, x = NULL, y = NULL) + theme_minimal() + 
    theme(panel.grid = element_blank(), axis.text = element_blank(), plot.title = element_text(size = 8.5, face="bold"), legend.title = element_text(size = 12, face="bold"), legend.text = element_text(size = 10))
    
  if (!show_legend) p <- p + theme(legend.position = "none")
  if (show_labels) {
    centroids <- df %>% group_by(.data[[join_col]]) %>% summarise(UMAP_1 = median(UMAP_1, na.rm = TRUE), UMAP_2 = median(UMAP_2, na.rm = TRUE), mod_id = if("mod_id_col" %in% names(.)) { v <- na.omit(mod_id_col); if(length(v) > 0) v[1] else NA_character_ } else { NA_character_ }, .groups = "drop") %>% filter(!is.na(mod_id) & mod_id != "") %>% mutate(label_text = .data[[join_col]])
    p <- p + geom_text_repel(data = centroids, aes(x = UMAP_1, y = UMAP_2, label = label_text), inherit.aes = FALSE, size = 4.5, fontface = "bold", bg.color = "white", bg.r = 0.15, color = "black", min.segment.length = 0)
  }
  return(p)
}

# Extracts all necessary data directly from Master Annotation based on module ID
get_module_betas <- function(df, mod_id, anno_col, join_col, gene_col) {
  sub_df <- df %>% filter(module == mod_id)
  
  if (nrow(sub_df) == 0) return(NULL) 
  if ("cs_max_pip" %in% names(sub_df)) sub_df <- sub_df %>% arrange(desc(as.numeric(cs_max_pip)))
  
  sub_df %>% 
    distinct(.data[[anno_col]], .keep_all = TRUE) %>%
    mutate(
      final_beta = if("cs_top_snp_beta" %in% names(.)) cs_top_snp_beta else if("most_likely_snp_beta" %in% names(.)) most_likely_snp_beta else beta, 
      gene_sym_col = if(gene_col %in% names(.)) .data[[gene_col]] else if("eGene" %in% names(.)) eGene else "Unknown_Gene", 
      mod_id_col = mod_id,
      snp_id = if("most_likely_snp_rsID" %in% names(.)) most_likely_snp_rsID else NA,
      pval_exact = if("most_likely_snp_chisq" %in% names(.)) pchisq(as.numeric(most_likely_snp_chisq), df=1, lower.tail=FALSE) else NA
    ) %>%
    select(all_of(anno_col), beta_val = final_beta, gene_sym_col, mod_id_col, snp_id, pval_exact) %>% 
    rename(!!join_col := all_of(anno_col))
}

# --- Main Logic ---
cat("Loading data...\n")

if (show_ref) {
  colors_df <- fread(opt$colors)
}

umap_1 <- fread(opt$umap_1)
if ("UMAP_0" %in% names(umap_1) && "UMAP_1" %in% names(umap_1)) umap_1 <- umap_1 %>% rename(UMAP_2 = UMAP_1, UMAP_1 = UMAP_0)
if (show_ref) {
  merged_1 <- umap_1 %>% left_join(colors_df %>% select(all_of(opt$join_col), hex_color = color_ct2), by = opt$join_col)
} else {
  merged_1 <- umap_1
}
if (nrow(merged_1) > opt$max_cells) { set.seed(42); merged_1 <- merged_1 %>% slice_sample(n = opt$max_cells) }
meta_1 <- fread(opt$master_1)

if (dual_mode) {
  umap_2 <- fread(opt$umap_2)
  if ("UMAP_0" %in% names(umap_2) && "UMAP_1" %in% names(umap_2)) umap_2 <- umap_2 %>% rename(UMAP_2 = UMAP_1, UMAP_1 = UMAP_0)
  if (show_ref) {
    merged_2 <- umap_2 %>% left_join(colors_df %>% select(all_of(opt$join_col), hex_color = color_ct2), by = opt$join_col)
  } else {
    merged_2 <- umap_2
  }
  if (nrow(merged_2) > opt$max_cells) { set.seed(42); merged_2 <- merged_2 %>% slice_sample(n = opt$max_cells) }
  meta_2 <- fread(opt$master_2)
}

cat("Generating plot grid...\n")
plot_list <- list()
n_items <- length(mods_1)

if (dual_mode && show_ref) {
  plot_list[[1]] <- plot_ref(merged_1, paste(opt$name_1, "Reference"), opt$join_col, opt$pt_size, use_raster)
  plot_list[[2]] <- plot_ref(merged_2, paste(opt$name_2, "Reference"), opt$join_col, opt$pt_size, use_raster)
}

for (i in 1:n_items) {
  # --- POP 1 LOGIC ---
  if (mods_1[i] != "NA") {
    beta_1 <- get_module_betas(meta_1, mods_1[i], anno_col, opt$join_col, opt$gene_col)
    
    if (!is.null(beta_1)) {
      plot_df_1 <- merged_1 %>% left_join(beta_1, by = opt$join_col) %>% mutate(beta_val = coalesce(as.numeric(beta_val), 0)) %>% arrange(abs(beta_val))
      
      gene_name_1 <- head(na.omit(plot_df_1$gene_sym_col), 1)
      if (length(gene_name_1) == 0) gene_name_1 <- "Unknown_Gene"
      
      snp_1 <- head(na.omit(plot_df_1$snp_id), 1)
      snp_str_1 <- if(length(snp_1) > 0 && !is.na(snp_1) && snp_1 != "") paste0("\n", snp_1) else "\nUnknown_SNP"
      
      pval_1 <- head(na.omit(plot_df_1$pval_exact), 1)
      pval_str_1 <- if(length(pval_1) > 0 && !is.na(pval_1)) sprintf(" | P = %.2e", as.numeric(pval_1)) else ""
      
      title_1 <- paste0(opt$name_1, " | ", gene_name_1, "\n[", mods_1[i], "]", snp_str_1, pval_str_1)
      plot_list[[length(plot_list) + 1]] <- plot_beta(plot_df_1, title_1, opt$pt_size, opt$join_col, show_legend = TRUE, use_raster = use_raster, show_labels = show_labels)
    } else {
      cat(paste("Warning: Module", mods_1[i], "not found in Pop 1 annotations - skipping panel.\n"))
    }
  } else {
    plot_list[[length(plot_list) + 1]] <- plot_spacer()
  }

  # --- POP 2 LOGIC ---
  if (dual_mode) {
    if (mods_2[i] != "NA") {
      beta_2 <- get_module_betas(meta_2, mods_2[i], anno_col, opt$join_col, opt$gene_col)
      
      if (!is.null(beta_2)) {
        plot_df_2 <- merged_2 %>% left_join(beta_2, by = opt$join_col) %>% mutate(beta_val = coalesce(as.numeric(beta_val), 0)) %>% arrange(abs(beta_val))
        
        gene_name_2 <- head(na.omit(plot_df_2$gene_sym_col), 1)
        if (length(gene_name_2) == 0) gene_name_2 <- "Unknown_Gene"
        
        snp_2 <- head(na.omit(plot_df_2$snp_id), 1)
        snp_str_2 <- if(length(snp_2) > 0 && !is.na(snp_2) && snp_2 != "") paste0("\n", snp_2) else "\nUnknown_SNP"
        
        pval_2 <- head(na.omit(plot_df_2$pval_exact), 1)
        pval_str_2 <- if(length(pval_2) > 0 && !is.na(pval_2)) sprintf(" | P = %.2e", as.numeric(pval_2)) else ""
        
        title_2 <- paste0(opt$name_2, " | ", gene_name_2, "\n[", mods_2[i], "]", snp_str_2, pval_str_2)
        plot_list[[length(plot_list) + 1]] <- plot_beta(plot_df_2, title_2, opt$pt_size, opt$join_col, show_legend = TRUE, use_raster = use_raster, show_labels = show_labels)
      } else {
        cat(paste("Warning: Module", mods_2[i], "not found in Pop 2 annotations - skipping panel.\n"))
      }
    } else {
      plot_list[[length(plot_list) + 1]] <- plot_spacer()
    }
  }
}

# --- Layout Logic ---
cols <- 2 

if (dual_mode) {
  n_rows <- length(plot_list) / cols
  final_plot <- wrap_plots(plot_list, ncol = cols)
  plot_width <- 14 
} else {
  n_rows <- ceiling(length(plot_list) / cols)
  grid_plot <- wrap_plots(plot_list, ncol = cols)
  
  if (show_ref) {
    p_ref <- plot_ref(merged_1, paste(opt$name_1, "Reference"), opt$join_col, opt$pt_size, use_raster)
    final_plot <- p_ref | grid_plot
    final_plot <- final_plot + plot_layout(widths = c(1, cols))
    plot_width <- (cols + 1) * 6 
  } else {
    final_plot <- grid_plot
    plot_width <- cols * 6
  }
}

base_height <- 4.5
plot_height <- max(base_height * n_rows, 5) 

# --- Save Logic ---
out_base <- sub("\\.png$|\\.pdf$", "", opt$out, ignore.case = TRUE)
if (save_pdf) {
  out_file <- paste0(out_base, ".pdf")
  cat(paste("Saving as vector PDF with rasterized points to:", out_file, "\n"))
  ggsave(out_file, plot = final_plot, width = plot_width, height = plot_height, device = cairo_pdf, limitsize = FALSE)
} else {
  out_file <- paste0(out_base, ".png")
  cat(paste("Saving as fast PNG to:", out_file, "\n"))
  ggsave(out_file, plot = final_plot, width = plot_width, height = plot_height, dpi = 300, bg = "white", device = ragg::agg_png, limitsize = FALSE)
}