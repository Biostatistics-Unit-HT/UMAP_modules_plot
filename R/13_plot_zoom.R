# Zoom annotation panel (per-CS) and the module-level merged annotation box.

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
                                  track_order, full_tbl = NULL) {
  win_s <- lead_pos - win_half
  win_e <- lead_pos + win_half
  n_tracks <- max(length(track_order), 1L)
  pal <- c("#1565C0", "#2E7D32", "#EF6C00", "#6A1B9A", "#00838F",
           "#AD1457", "#558B2F", "#4E342E", "#283593", "#F57F17")
  track_colors <- setNames(rep(pal, length.out = n_tracks), track_order)
  
  # Pre-compute which lanes render as continuous profiles vs discrete blocks.
  continuous <- vapply(track_order, function(ft) is_continuous_track(ft, full_tbl),
                       logical(1))
  
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
    
    for (k in seq_along(track_order)) {
      ft <- track_order[k]
      rows <- anno_win[anno_win$feature_type == ft, , drop = FALSE]
      if (nrow(rows) == 0) next
      if (continuous[k]) {
        # Per-lane min-max normalization (inside the window). Bottom of the
        # lane is the baseline; bars grow upward up to 80% of lane height.
        max_score <- suppressWarnings(max(rows$score, na.rm = TRUE))
        if (!is.finite(max_score) || max_score <= 0) max_score <- 1
        rows$score_norm <- ifelse(is.na(rows$score), 0,
                                  pmin(rows$score / max_score, 1))
        p <- p + geom_rect(
          data = rows,
          aes(xmin = x_start, xmax = x_end,
              ymin = lane - 0.4,
              ymax = lane - 0.4 + 0.8 * score_norm,
              fill = feature_type),
          color = NA, alpha = 0.9)
      } else {
        p <- p + geom_rect(
          data = rows,
          aes(xmin = x_start, xmax = x_end,
              ymin = lane - 0.35, ymax = lane + 0.35,
              fill = feature_type),
          color = NA, alpha = 0.85)
      }
    }
  } else {
    p <- p + annotate("text", x = lead_pos, y = (n_tracks + 1) / 2,
                      label = "No annotations in window",
                      size = 3.2, color = "gray55", fontface = "italic")
  }
  p + geom_vline(xintercept = lead_pos, linetype = "dashed",
                 color = "#8E24AA", alpha = 0.8, linewidth = 0.5)
}

# build_merged_zoom_box(): module-level annotation box covering every
# credible-set lead SNP on a single continuous x-axis. Used when a module
# has >=1 credible set and the user supplied --annotations.
#
# Inputs:
#   chr_label        : chromosome (no "chr" prefix) for filtering annotations.
#   lead_positions   : numeric vector of CS lead-SNP positions (bp).
#   cs_labels        : optional character vector of CS identifiers; when
#                      supplied, each lead gets a short label "CS1 | L..."
#                      floated above the dashed vertical.
#   annotations_tbl  : full annotation data.table (chrom, start, end,
#                      feature_type, score).
#   annotation_tracks: character vector giving the shared lane order.
#   zoom_window      : half-width in bp around each lead; the final x-range
#                      is [min(lead)-zoom, max(lead)+zoom].
#   title            : optional panel title.
#
# Returns a ggplot. When no annotations overlap the merged window a small
# placeholder panel is still produced so the per-module layout is stable.
#
# Example:
#   build_merged_zoom_box("18", c(80005273, 80120000),
#                          cs_labels = c("...::L1", "...::L2"),
#                          annotations_tbl, annotation_tracks,
#                          zoom_window = 5000,
#                          title = "Merged annotations | M_18448")
build_merged_zoom_box <- function(chr_label, lead_positions,
                                  cs_labels = NULL,
                                  annotations_tbl, annotation_tracks,
                                  zoom_window,
                                  title = NULL,
                                  extra_lead_pos   = NA_real_,
                                  extra_lead_label = NA_character_) {
  lead_positions <- sort(unique(as.numeric(lead_positions)))
  lead_positions <- lead_positions[!is.na(lead_positions) &
                                    is.finite(lead_positions)]
  if (length(lead_positions) == 0) return(plot_spacer())
  
  # Extend the window so the module-level "most likely SNP" tick fits
  # even when it sits outside the union of CS-specific leads.
  all_positions <- lead_positions
  if (!is.na(extra_lead_pos) && is.finite(extra_lead_pos))
    all_positions <- c(all_positions, as.numeric(extra_lead_pos))
  win_s <- min(all_positions) - zoom_window
  win_e <- max(all_positions) + zoom_window
  n_tracks <- max(length(annotation_tracks), 1L)
  
  anno_win <- filter_annotations_window(annotations_tbl, chr_label, win_s, win_e)
  
  pal <- c("#1565C0", "#2E7D32", "#EF6C00", "#6A1B9A", "#00838F",
           "#AD1457", "#558B2F", "#4E342E", "#283593", "#F57F17")
  track_colors <- setNames(rep(pal, length.out = n_tracks), annotation_tracks)
  
  continuous <- vapply(annotation_tracks,
                       function(ft) is_continuous_track(ft, annotations_tbl),
                       logical(1))
  
  panel_title <- if (is.null(title)) {
    sprintf("Merged annotations (chr%s, %d credible set%s, %s bp span)",
            chr_label, length(lead_positions),
            if (length(lead_positions) > 1) "s" else "",
            format(win_e - win_s, big.mark = ",", scientific = FALSE))
  } else title
  
  p <- ggplot() +
    scale_y_continuous(breaks = seq_along(annotation_tracks),
                       labels = annotation_tracks,
                       limits = c(0.3, n_tracks + 1.6)) +
    scale_fill_manual(values = track_colors, guide = "none",
                      limits = annotation_tracks) +
    scale_x_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
    coord_cartesian(xlim = c(win_s, win_e), expand = FALSE) +
    labs(title = panel_title,
         x = paste0("Chromosome ", chr_label, " position (bp)"),
         y = NULL) +
    theme_minimal() +
    theme(plot.title         = element_text(size = 10, face = "bold"),
          axis.title.x       = element_text(size = 10, face = "bold"),
          axis.text          = element_text(size = 9),
          panel.grid.minor   = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.border       = element_rect(color = "gray80", fill = NA, linewidth = 0.8))
  
  if (!is.null(anno_win) && nrow(anno_win) > 0) {
    anno_win <- as.data.frame(anno_win)
    anno_win$lane    <- match(anno_win$feature_type, annotation_tracks)
    anno_win$x_start <- pmax(anno_win$start, win_s)
    anno_win$x_end   <- pmin(anno_win$end,   win_e)
    anno_win <- anno_win[!is.na(anno_win$lane), , drop = FALSE]
    
    for (k in seq_along(annotation_tracks)) {
      ft <- annotation_tracks[k]
      rows <- anno_win[anno_win$feature_type == ft, , drop = FALSE]
      if (nrow(rows) == 0) next
      if (continuous[k]) {
        max_score <- suppressWarnings(max(rows$score, na.rm = TRUE))
        if (!is.finite(max_score) || max_score <= 0) max_score <- 1
        rows$score_norm <- ifelse(is.na(rows$score), 0,
                                  pmin(rows$score / max_score, 1))
        p <- p + geom_rect(data = rows,
                           aes(xmin = x_start, xmax = x_end,
                               ymin = lane - 0.4,
                               ymax = lane - 0.4 + 0.8 * score_norm,
                               fill = feature_type),
                           color = NA, alpha = 0.9)
      } else {
        p <- p + geom_rect(data = rows,
                           aes(xmin = x_start, xmax = x_end,
                               ymin = lane - 0.35, ymax = lane + 0.35,
                               fill = feature_type),
                           color = NA, alpha = 0.85)
      }
    }
  } else {
    p <- p + annotate("text",
                      x = mean(lead_positions), y = (n_tracks + 1) / 2,
                      label = "No annotations in merged window",
                      size = 3.2, color = "gray55", fontface = "italic")
  }
  
  # Per-CS lead SNP dashed lines were removed at the user's request. Only
  # the module-level "most-likely SNP" line (below) is drawn on the merged
  # annotation panel. `lead_positions` / `cs_labels` are still used above
  # to set the window width so no data is cropped.
  
  # Marker: the module-level "most likely SNP" from the master annotation
  # (or the top-P SNP in the LZ file when the master was not provided).
  # Drawn as a dotted orange vertical line with a label above the tracks.
  if (!is.na(extra_lead_pos) && is.finite(extra_lead_pos)) {
    extra_df <- data.frame(
      lead  = as.numeric(extra_lead_pos),
      label = if (!is.na(extra_lead_label) && nzchar(extra_lead_label))
                as.character(extra_lead_label)
              else sprintf("most-likely SNP | %s",
                           format(extra_lead_pos, big.mark = ",", scientific = FALSE)),
      stringsAsFactors = FALSE
    )
    p <- p + geom_vline(data = extra_df, aes(xintercept = lead),
                        linetype = "dotted", color = "#EF6C00",
                        alpha = 0.95, linewidth = 0.7)
    p <- p + geom_text_repel(data = extra_df,
                             aes(x = lead, y = n_tracks + 1.3, label = label),
                             color = "#BF4F00", size = 2.7, fontface = "bold",
                             direction = "y", nudge_y = 0.2,
                             segment.size = 0.3, segment.alpha = 0.6,
                             box.padding = 0.25, point.padding = 0.1,
                             min.segment.length = 0, max.overlaps = Inf,
                             seed = 43)
  }
  p
}
