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

# beta_brbg_colors(): ColorBrewer "BrBG" (Brown–Blue Green) diverging stops
# for the Beta UMAP. Low beta -> brown, 0 -> neutral, high beta -> teal/green.
#
# Example:
#   beta_brbg_colors()
#   # -> character vector of 11 hex colours (low -> high)
beta_brbg_colors <- function() {
  c("#543005", "#8C510A", "#BF812D", "#DFC27D", "#F6E8C3", "#F5F5F5",
    "#C7EAE5", "#80CDC1", "#35978F", "#01665E", "#003C30")
}

# plot_beta(): UMAP coloured by a per-cell-type beta_val, with a symmetric
# BrBG (brown–blue-green) diverging palette anchored at 0. Rows with NA
# beta_val use na.value (grey) so non-highlighted cell types stay in the
# background.
#
# Arguments:
#   join_col  - column used for FILTERING (matches `label_cells_only`, which
#               is what the orchestrator passes for the focal cell). Always
#               needed.
#   label_col - column used for the centroid LABELS drawn on top of the
#               panel. Defaults to `join_col` (i.e. labels and filter come
#               from the same column). When set to a different column (e.g.
#               "celltype_3"), centroids are grouped by `label_col` and only
#               those `label_col` groups containing at least one row whose
#               `join_col` is in `label_cells_only` are labelled.
#
# Example:
#   plot_beta(plot_df, "GH | RBFA\n[M_18448]", 0.25, "celltype_2",
#             show_legend=TRUE, use_raster=TRUE, show_labels=TRUE,
#             label_cells_only=NULL)
#   # -> ggplot; with label_cells_only=c("B_naive") only that type gets a label.
#   plot_beta(..., join_col="celltype_2", label_col="celltype_3",
#             label_cells_only=c("NK_CD16_adaptive"))
#   # -> labels are drawn at the celltype_3 granularity but restricted to
#   #    groups that contain at least one celltype_2 == "NK_CD16_adaptive".
plot_beta <- function(df, title, pt_size, join_col, show_legend = TRUE,
                      use_raster = FALSE, show_labels = FALSE,
                      label_cells_only = NULL, label_col = NULL) {
  b_max <- max(abs(df$beta_val), na.rm = TRUE)
  if (!is.finite(b_max) || b_max <= 0) b_max <- 1e-9
  colors_beta <- beta_brbg_colors()
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
    # UMAP_1 / UMAP_2 axes carry no numeric meaning, so remove their breaks
    # at the SCALE level. Doing it via theme(axis.text = element_blank())
    # alone is not enough because the outer `& big_helvetica_theme()`
    # cascade later overrides axis.text and resurrects the tick labels.
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    # `clip = "off"` lets the ggrepel cell-type centroid labels draw past
    # the panel edge -- at the bigger Helvetica size they sometimes start
    # near the panel boundary and get truncated (e.g. "T_CD4_C" instead of
    # "T_CD4_CM").
    coord_equal(clip = "off") + labs(title = title, x = NULL, y = NULL) + theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank(),
          plot.title  = element_text(size = 8.5, face = "bold"),
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 10))
  
  if (!show_legend) p <- p + theme(legend.position = "none")
  if (show_labels) {
    # `group_col` controls where labels are anchored. Defaults to join_col
    # (existing behaviour); when `label_col` is set AND present in the data
    # frame, labels use that finer/coarser annotation instead.
    group_col <- if (!is.null(label_col) && nzchar(label_col) &&
                     label_col %in% names(df)) label_col else join_col
    centroids <- df %>% group_by(.data[[group_col]]) %>%
      summarise(UMAP_1 = median(UMAP_1, na.rm = TRUE),
                UMAP_2 = median(UMAP_2, na.rm = TRUE),
                mod_id = if ("mod_id_col" %in% names(.)) { v <- na.omit(mod_id_col); if (length(v) > 0) v[1] else NA_character_ } else { NA_character_ },
                .groups = "drop") %>%
      filter(!is.na(mod_id) & mod_id != "") %>%
      mutate(label_text = .data[[group_col]])
    if (!is.null(label_cells_only) && length(label_cells_only) > 0) {
      if (identical(group_col, join_col)) {
        # Labels and filter share a column -> direct match.
        keep_lab <- as.character(centroids[[group_col]]) %in% label_cells_only
      } else {
        # Different granularity: keep label-column groups that contain at
        # least one row whose join_col is in label_cells_only.
        keep_groups <- df %>%
          dplyr::filter(as.character(.data[[join_col]]) %in% label_cells_only) %>%
          dplyr::pull(.data[[group_col]]) %>% as.character() %>% unique()
        keep_lab <- as.character(centroids[[group_col]]) %in% keep_groups
      }
      centroids <- centroids[keep_lab, , drop = FALSE]
    }
    p <- p + geom_text_repel(data = centroids,
                             aes(x = UMAP_1, y = UMAP_2, label = label_text),
                             inherit.aes = FALSE, size = 7.5,
                             family = "Helvetica", fontface = "plain",
                             bg.color = "white", bg.r = 0.15,
                             color = "black", min.segment.length = 0)
  }
  return(p)
}
