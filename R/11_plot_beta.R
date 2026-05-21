# UMAP coloured by per-cell-type beta for the focal (module, cell, gene).

# scattermore_pointsize(): map ggplot-style --pt_size to geom_scattermore's
# pointsize (pixels on the internal raster); small values were nearly invisible.
#
# Example:
#   scattermore_pointsize(0.25)  # -> ~10 pixels
#   scattermore_pointsize(0.6)   # -> larger on-screen dots
scattermore_pointsize <- function(pt_size) {
  max(4, 6 + 22 * as.numeric(pt_size))
}

# beta_prgn_colors(): ColorBrewer "PrGn" (Purple–Green) diverging stops.
# Low / negative beta -> purple; high / positive beta -> green.
#
# Example:
#   beta_prgn_colors()
#   # -> character vector of 11 hex colours (purple -> green)
beta_prgn_colors <- function() {
  c("#40004B", "#762A83", "#9970AB", "#C2A5CF", "#E7D4E8", "#F7F7F7",
    "#D9F0D3", "#A6DBA0", "#5AAE61", "#1B7837", "#00441B")
}

# plot_beta(): UMAP coloured by beta_val with a symmetric PrGn diverging
# palette (purple = negative, green = positive, centred at 0). Non-focal
# cells use na.value (grey). Panel title is NULL; the orchestrator supplies
# the row-level figure title (module, cell, gene, SNP, beta).
#
# Arguments:
#   join_col  - column used for FILTERING (matches `label_cells_only`).
#   label_col - column used for centroid labels (defaults to join_col).
#
# Example:
#   plot_beta(plot_df, NULL, 0.25, "celltype_2", show_legend = TRUE,
#             label_cells_only = "T_CD4_CM")
#   # -> ggplot; focal cell on PrGn scale, others grey, Beta legend shown.
plot_beta <- function(df, title, pt_size, join_col, show_legend = TRUE,
                      use_raster = FALSE, show_labels = FALSE,
                      label_cells_only = NULL, label_col = NULL) {
  b_max <- max(abs(df$beta_val), na.rm = TRUE)
  if (!is.finite(b_max) || b_max <= 0) b_max <- 1e-9
  colors_beta <- beta_prgn_colors()
  n_beta_cols <- length(colors_beta)

  p <- ggplot(df, aes(x = UMAP_1, y = UMAP_2, color = beta_val))
  if (use_raster) {
    p <- p + geom_scattermore(pointsize = scattermore_pointsize(pt_size),
                              pixels = c(2048, 2048))
  } else {
    p <- p + geom_point(size = max(0.15, 1.2 + 5 * as.numeric(pt_size)), stroke = 0)
  }

  p <- p + scale_color_gradientn(colors = colors_beta,
                                 limits = c(-b_max, b_max), na.value = "#dadada",
                                 name = "Beta",
                                 values = seq(0, 1, length.out = n_beta_cols)) +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    coord_equal(clip = "off") +
    labs(title = NULL, x = NULL, y = NULL) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank(),
          legend.title = element_text(size = 12, face = "plain"),
          legend.text  = element_text(size = 10))

  if (!show_legend) p <- p + theme(legend.position = "none")

  if (show_labels) {
    group_col <- if (!is.null(label_col) && nzchar(label_col) &&
                     label_col %in% names(df)) label_col else join_col
    centroids <- df %>% group_by(.data[[group_col]]) %>%
      summarise(UMAP_1 = median(UMAP_1, na.rm = TRUE),
                UMAP_2 = median(UMAP_2, na.rm = TRUE),
                mod_id = if ("mod_id_col" %in% names(.)) {
                  v <- na.omit(mod_id_col)
                  if (length(v) > 0) v[1] else NA_character_
                } else {
                  NA_character_
                },
                .groups = "drop") %>%
      filter(!is.na(mod_id) & mod_id != "") %>%
      mutate(label_text = .data[[group_col]])
    if (!is.null(label_cells_only) && length(label_cells_only) > 0) {
      if (identical(group_col, join_col)) {
        keep_lab <- as.character(centroids[[group_col]]) %in% label_cells_only
      } else {
        keep_groups <- df %>%
          dplyr::filter(as.character(.data[[join_col]]) %in% label_cells_only) %>%
          dplyr::pull(.data[[group_col]]) %>% as.character() %>% unique()
        keep_lab <- as.character(centroids[[group_col]]) %in% keep_groups
      }
      centroids <- centroids[keep_lab, , drop = FALSE]
    }
    if (nrow(centroids) > 0) {
      p <- p + geom_text_repel(data = centroids,
                               aes(x = UMAP_1, y = UMAP_2, label = label_text),
                               inherit.aes = FALSE, size = 7.5,
                               family = "Helvetica", fontface = "plain",
                               bg.color = "white", bg.r = 0.15,
                               color = "black", min.segment.length = 0)
    }
  }
  return(p)
}
