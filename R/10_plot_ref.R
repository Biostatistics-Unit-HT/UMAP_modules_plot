# Reference UMAP coloured by cell-type palette (enabled via --show_ref).

# plot_ref(): reference UMAP coloured by celltype palette; labels celltype
# centroids with ggrepel. Used when --show_ref is enabled. Raster point size
# uses scattermore_pointsize() from 11_plot_beta.R (sourced later in the same run).
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
  if (use_raster) {
    p <- p + geom_scattermore(aes(color = hex_color),
                              pointsize = scattermore_pointsize(pt_size),
                              pixels = c(2048, 2048))
  } else {
    p <- p + geom_point(aes(color = hex_color),
                        size = max(0.15, 1.2 + 5 * as.numeric(pt_size)), stroke = 0)
  }
  p + scale_color_identity() +
    geom_text_repel(data = centroids, aes(label = .data[[join_col]]),
                    family = "Helvetica", fontface = "plain",
                    size = 6, bg.color = "white",
                    bg.r = 0.1, min.segment.length = 0) +
    # UMAP dimensions carry no numeric meaning -- remove breaks at the scale
    # level so the global theme cascade can't reintroduce tick labels.
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    coord_equal() + labs(title = title_text, x = NULL, y = NULL) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text  = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(face = "bold"))
}
