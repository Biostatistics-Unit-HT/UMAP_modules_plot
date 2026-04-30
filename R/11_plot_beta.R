# UMAP coloured by per-cell-type beta for the focal (module, cell, gene).

# plot_beta(): UMAP coloured by a per-cell-type beta_val, with a symmetric
# blue-yellow-red diverging palette anchored at 0. Rows with NA beta_val use
# na.value (grey) so non-highlighted cell types can stay in the background.
#
# Example:
#   plot_beta(plot_df, "GH | RBFA\n[M_18448]", 0.25, "celltype_2",
#             show_legend=TRUE, use_raster=TRUE, show_labels=TRUE,
#             label_cells_only=NULL)
#   # -> ggplot; with label_cells_only=c("B_naive") only that type gets a label.
plot_beta <- function(df, title, pt_size, join_col, show_legend = TRUE,
                      use_raster = FALSE, show_labels = FALSE,
                      label_cells_only = NULL) {
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
    if (!is.null(label_cells_only) && length(label_cells_only) > 0) {
      keep_lab <- as.character(centroids[[join_col]]) %in% label_cells_only
      centroids <- centroids[keep_lab, , drop = FALSE]
    }
    p <- p + geom_text_repel(data = centroids,
                             aes(x = UMAP_1, y = UMAP_2, label = label_text),
                             inherit.aes = FALSE, size = 4.5,
                             fontface = "bold", bg.color = "white", bg.r = 0.15,
                             color = "black", min.segment.length = 0)
  }
  return(p)
}
