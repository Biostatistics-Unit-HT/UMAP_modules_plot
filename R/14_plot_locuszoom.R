# LocusZoom composite panel: -log10(P) scatter + gene track + optional
# per-CS zoom panel + optional connector strip. Uses helpers from
# 01_io_gtf.R, 02_io_annotations.R, 03_io_lz.R, 04_io_ld.R, 13_plot_zoom.R.

# clip_gene_segments_to_xlim(): clip gene arrow endpoints and label anchors to
# the LocusZoom Mb window so geom_segment / ggrepel stay inside the panel.
#
# Example:
#   df <- data.frame(seg_x = 79.99, seg_xend = 80.02, mid_mb = 80.005, lane = 1)
#   clip_gene_segments_to_xlim(df, c(80.0, 80.05))
clip_gene_segments_to_xlim <- function(track_df, x_lims_mb, eps_mb = 1e-6) {
  x_lo <- x_lims_mb[1]
  x_hi <- x_lims_mb[2]
  track_df$seg_x_clip <- pmax(pmin(track_df$seg_x, x_hi), x_lo)
  track_df$seg_xend_clip <- pmax(pmin(track_df$seg_xend, x_hi), x_lo)
  track_df <- track_df[abs(track_df$seg_x_clip - track_df$seg_xend_clip) > eps_mb, ,
                       drop = FALSE]
  # Inset label anchors slightly from the boundary so ggrepel connector
  # segments never originate exactly on (or outside) the panel border.
  label_inset <- (x_hi - x_lo) * 0.01
  track_df$mid_mb_clip <- pmax(pmin(track_df$mid_mb, x_hi - label_inset),
                                x_lo + label_inset)
  track_df
}

# plot_locuszoom(): composite LocusZoom panel.
#   - main: -log10(P) vs position (Mb) with lead SNP highlighted
#   - gene track: all protein-coding genes in the window (focal eGene red)
#   - optional zoom: regulatory annotations within +/- zoom_window bp
#   - xlim_mode "context" widens x-axis with eGene + GTF genes; "snp" uses
#     association POS range only (+ small pad).
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
                           zoom_window = 5000,
                           lz_filter_cell = NULL,
                           lz_filter_gene = NULL,
                           lz_filter_cs   = NULL,
                           include_zoom_panel = TRUE,
                           include_zoom_connector = TRUE,
                           ld_vec = NULL,
                           force_lead_pos = NULL,
                           merge_all_cs = FALSE,
                           xlim_mode = c("context", "snp")) {
  xlim_mode <- match.arg(xlim_mode)
  if (is.na(lz_file) || lz_file == "NA" || !file.exists(lz_file)) {
    if (!is.na(lz_file) && lz_file != "NA") cat(sprintf("Warning: LocusZoom file not found at %s\n", lz_file))
    return(plot_spacer())
  }
  df <- load_lz_file(lz_file)
  if (is.null(df)) return(plot_spacer())
  needed <- c("CHR", "POS", "P")
  if (!all(needed %in% names(df))) {
    cat(sprintf("Warning: %s is missing required columns (CHR, POS, P)\n", lz_file))
    return(plot_spacer())
  }
  # Optionally subset to a specific (cell, gene) context when the file
  # pools associations for multiple combinations.
  if (!is.null(lz_filter_cell) && "CELL" %in% names(df)) {
    df <- df[as.character(CELL) == lz_filter_cell]
  }
  if (!is.null(lz_filter_gene) && "GENE" %in% names(df)) {
    df <- df[sub("\\.\\d+$", "", as.character(GENE)) == sub("\\.\\d+$", "", lz_filter_gene)]
  }
  if (!is.null(lz_filter_cs) && "CS" %in% names(df) && !isTRUE(merge_all_cs)) {
    df <- df[as.character(CS) == lz_filter_cs]
  }
  df <- df[!is.na(P) & is.finite(as.numeric(P))]
  if (nrow(df) == 0) {
    cat(sprintf("Warning: %s has no rows for cell=%s gene=%s\n", lz_file,
                if (is.null(lz_filter_cell)) "*" else lz_filter_cell,
                if (is.null(lz_filter_gene)) "*" else lz_filter_gene))
    return(plot_spacer())
  }
  # `P` is the raw p-value; compute -log10 and floor to 0 to avoid -Inf.
  df[, logp := -log10(pmax(as.numeric(P), 1e-300))]
  
  if (isTRUE(merge_all_cs)) {
    if (!"CS" %in% names(df)) {
      cat(sprintf("Warning: merge_all_cs requires a CS column in %s; using standard single-region panel.\n", lz_file))
    } else {
      if (!is.null(ld_vec)) {
        cat("Note: LD colouring is disabled when merge_all_cs=TRUE (points coloured by credible set).\n")
        ld_vec <- NULL
      }
      chr_lab_m <- sub("^chr", "", as.character(df$CHR[1]), ignore.case = TRUE)
      return(plot_locuszoom_merged(
        df, locus_info, title, genes_df, annotations_tbl, annotation_tracks,
        zoom_window, include_zoom_panel, include_zoom_connector,
        chr_label = chr_lab_m, xlim_mode = xlim_mode))
    }
  }
  
  chr_label <- sub("^chr", "", as.character(df$CHR[1]), ignore.case = TRUE)
  
  # Preferred order for the lead SNP position:
  #   1. `force_lead_pos` when provided (used by the orchestrator to pin
  #      every LocusZoom in a module to the master's module-level
  #      most-likely SNP so the diamond / LD colouring are consistent),
  #   2. extracted from the CS column of this LZ file when present,
  #   3. locus_info$lead_pos from the master annotation,
  #   4. min-P fallback inside the file.
  cs_lead_pos <- NA_real_
  if ("CS" %in% names(df)) {
    cs_val <- as.character(df$CS[1])
    # Pull any chr<N>:<pos>:<ALLELE>:<ALLELE> token out of the CS string.
    m <- regmatches(cs_val,
                    regexpr("chr[0-9XYM]+:\\d+:[ACGTN]+:[ACGTN]+",
                            cs_val, ignore.case = TRUE))
    if (length(m) == 1 && nchar(m) > 0) {
      cs_lead_pos <- suppressWarnings(as.numeric(strsplit(m, ":")[[1]][2]))
    }
  }
  
  lead_pos <- if (!is.null(force_lead_pos) && !is.na(force_lead_pos) &&
                  is.finite(as.numeric(force_lead_pos)))
                as.numeric(force_lead_pos)
              else if (!is.na(cs_lead_pos)) cs_lead_pos
              else if (!is.null(locus_info) && !is.null(locus_info$lead_pos))
                locus_info$lead_pos
              else NA_real_
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
  # context: widen x-axis with focal eGene + all genes in genes_df (GTF).
  # snp: min/max association POS only (still adds small pad below).
  if (identical(xlim_mode, "context")) {
    if (has_focal) {
      x_min_bp <- min(x_min_bp, as.numeric(locus_info$gene_start))
      x_max_bp <- max(x_max_bp, as.numeric(locus_info$gene_end))
    }
    if (has_many) {
      x_min_bp <- min(x_min_bp, min(genes_df$start))
      x_max_bp <- max(x_max_bp, max(genes_df$end))
    }
  }
  pad <- max((x_max_bp - x_min_bp) * 0.02, 1)
  x_lims_mb <- c((x_min_bp - pad) / 1e6, (x_max_bp + pad) / 1e6)
  y_max <- max(df$logp, na.rm = TRUE) * 1.15
  
  # Attach LD r^2 from the lead SNP to every SNP in the panel (by chr:pos).
  if (!is.null(ld_vec)) {
    df[, ld_r2  := ld_vec[paste(sub("^chr","",CHR,ignore.case=TRUE),
                                 as.integer(POS), sep=":")]]
    df[, ld_bin := bin_ld_r2(ld_r2)]
  } else {
    df[, ld_bin := factor(NA_character_, levels = names(LD_COLORS))]
  }
  
  p_main <- ggplot(df[is_lead == FALSE], aes(x = POS / 1e6, y = logp))
  if (!is.null(ld_vec)) {
    # One point size for every SNP; only fill reflects LD bin. Draw
    # lower-LD points first so warmer colours sit on top.
    PT_LZ <- 2.5
    df_nonlead <- data.table::copy(df[is_lead == FALSE][order(as.integer(ld_bin),
                                                              na.last = FALSE)])
    # Invisible anchor layer: one off-screen point per bin so every LD
    # colour appears in the fill legend even when a bin has zero SNPs.
    anchor_df <- data.frame(
      POS    = rep(min(df$POS, na.rm = TRUE) - 1e9,
                   length(LD_COLORS)),
      logp   = rep(-1, length(LD_COLORS)),
      ld_bin = factor(names(LD_COLORS), levels = names(LD_COLORS)),
      stringsAsFactors = FALSE)
    p_main <- ggplot(df_nonlead, aes(x = POS / 1e6, y = logp)) +
      geom_point(data = anchor_df,
                 aes(x = POS / 1e6, y = logp, fill = ld_bin),
                 shape = 21, color = "white", stroke = 0.3,
                 size = PT_LZ, inherit.aes = FALSE, show.legend = TRUE) +
      geom_point(aes(fill = ld_bin), shape = 21,
                 color = "white", alpha = 0.9, stroke = 0.3, size = PT_LZ) +
      scale_fill_manual(name = expression(r^2), values = LD_COLORS,
                        drop = FALSE, limits = names(LD_COLORS),
                        breaks = names(LD_COLORS),
                        guide = guide_legend(
                          override.aes = list(
                            shape  = 21,
                            size   = PT_LZ,
                            stroke = 0.3,
                            color  = "white",
                            fill   = unname(LD_COLORS))))
  } else {
    p_main <- p_main +
      geom_point(shape = 21, fill = "#4F8EDC", color = "white",
                 size = 1.6, alpha = 0.85, stroke = 0.25)
  }
  p_main <- p_main +
    geom_point(data = df[is_lead == TRUE][1L], aes(x = POS / 1e6, y = logp),
               shape = 23, fill = "#8E24AA", color = "black",
               size = 3.5, stroke = 0.6, inherit.aes = FALSE) +
    geom_text_repel(data = df[is_lead == TRUE][1L],
                    aes(x = POS / 1e6, y = logp, label = lead_label),
                    size = 2.9, fontface = "bold", color = "#4A148C",
                    nudge_y = y_max * 0.08, segment.color = "#8E24AA",
                    min.segment.length = 0, box.padding = 0.3,
                    inherit.aes = FALSE) +
    coord_cartesian(xlim = x_lims_mb, ylim = c(0, y_max)) +
    labs(title = title, x = NULL, y = expression(-log[10](italic(P)))) +
    theme_minimal() +
    theme(plot.title = element_text(size = 10, face = "bold"),
          # Pull the "-log10(P)" title tight against the tick labels.
          axis.title.y = element_text(size = 10, face = "bold",
                                      margin = margin(r = 2)),
          axis.text = element_text(size = 9),
          axis.text.y = element_text(margin = margin(r = 1)),
          axis.text.x = element_blank(),
          axis.ticks.length = unit(0.12, "cm"),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.8),
          plot.margin = margin(5, 5, 0, 2),
          legend.position = if (!is.null(ld_vec)) "right" else "none",
          legend.title = element_text(size = 8, face = "bold"),
          legend.text  = element_text(size = 7),
          legend.key.size = unit(0.35, "cm"))
  
  build_zoom_panel <- function() {
    if (!include_zoom_panel) return(NULL)
    if (is.null(annotations_tbl) || is.null(annotation_tracks) ||
        length(annotation_tracks) == 0) return(NULL)
    anno_win <- filter_annotations_window(annotations_tbl, chr_label,
                                          lead_pos - zoom_window,
                                          lead_pos + zoom_window)
    plot_zoom_annotations(anno_win, chr_label, lead_pos, zoom_window,
                          annotation_tracks, full_tbl = annotations_tbl)
  }
  
  # Two diagonal dashed lines that visually link the zoom-window edges on
  # the main LocusZoom / gene-track x-range (in Mb) down to the full-width
  # left/right edges of the zoom annotation panel below. Because patchwork
  # aligns panel plot-areas, using `coord_cartesian(xlim = x_lims_mb)` lets
  # us draw the top of the lines at the exact Mb position of the zoom
  # window while the bottom sits at the panel's full width -- where the
  # zoom panel starts.
  build_zoom_connector <- function() {
    ggplot() +
      geom_segment(aes(x = (lead_pos - zoom_window) / 1e6,
                       xend = x_lims_mb[1],
                       y = 1, yend = 0),
                   linewidth = 0.55, color = "#8E24AA",
                   linetype = "dashed", alpha = 0.75) +
      geom_segment(aes(x = (lead_pos + zoom_window) / 1e6,
                       xend = x_lims_mb[2],
                       y = 1, yend = 0),
                   linewidth = 0.55, color = "#8E24AA",
                   linetype = "dashed", alpha = 0.75) +
      coord_cartesian(xlim = x_lims_mb, ylim = c(0, 1), expand = FALSE) +
      theme_void() +
      theme(plot.margin = margin(0, 5, 0, 5))
  }
  
  if (!has_focal && !has_many) {
    p_main <- p_main +
      labs(x = paste0("Chromosome ", chr_label, " position (Mb)")) +
      theme(axis.text.x = element_text(size = 9),
            axis.title.x = element_text(size = 10, face = "bold"))
    p_zoom <- build_zoom_panel()
    if (is.null(p_zoom)) return(p_main)
    p_link <- if (include_zoom_connector) build_zoom_connector() else NULL
    zoom_h <- max(1.2, 0.6 * length(annotation_tracks) + 0.6)
    if (is.null(p_link))
      return(p_main / p_zoom + plot_layout(heights = c(5, zoom_h)))
    return(p_main / p_link / p_zoom +
             plot_layout(heights = c(5, 0.35, zoom_h)))
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
  track_df$seg_x    <- ifelse(track_df$strand == "-", track_df$end,   track_df$start) / 1e6
  track_df$seg_xend <- ifelse(track_df$strand == "-", track_df$start, track_df$end)   / 1e6
  track_df$mid_mb   <- (track_df$start + track_df$end) / 2 / 1e6
  track_df$col      <- ifelse(track_df$is_focal, "#C62828", "#1A237E")
  track_df <- clip_gene_segments_to_xlim(track_df, x_lims_mb)
  if (nrow(track_df) == 0L) {
    n_lanes <- 1L
    p_gene <- ggplot() +
      coord_cartesian(xlim = x_lims_mb, ylim = c(0.3, 1.2),
                      clip = "on", expand = FALSE) +
      labs(x = paste0("Chromosome ", chr_label, " position (Mb)"), y = NULL) +
      theme_minimal() +
      theme(panel.grid = element_blank(),
            axis.text.y = element_blank(),
            axis.title.x = element_text(size = 10, face = "bold"),
            axis.text.x = element_text(size = 9),
            panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.8),
            plot.margin = margin(0, 5, 5, 5))
  } else {
    n_lanes <- max(track_df$lane, 1L)
    p_gene <- ggplot(track_df) +
      geom_segment(aes(x = seg_x_clip, xend = seg_xend_clip, y = lane, yend = lane,
                       color = I(col)),
                   linewidth = 1.1,
                   arrow = arrow(ends = "last", type = "closed",
                                 length = unit(0.07, "inches"))) +
      geom_text_repel(aes(x = mid_mb_clip, y = lane + 0.35, label = gene_name,
                          color = I(col),
                          fontface = ifelse(is_focal, "bold.italic", "italic")),
                      size = 2.8,
                      direction = "y",
                      nudge_y = 0.25,
                      box.padding = 0.2,
                      point.padding = 0.1,
                      segment.size = 0.45,
                      segment.alpha = 0.75,
                      min.segment.length = 0,
                      max.overlaps = Inf,
                      seed = 42,
                      xlim = c(x_lims_mb[1], x_lims_mb[2]),
                      ylim = c(0.35, n_lanes + 1.55)) +
      scale_y_continuous(limits = c(0.1, n_lanes + 1.7), breaks = NULL) +
      coord_cartesian(xlim = x_lims_mb, ylim = c(0.1, n_lanes + 1.7),
                      clip = "on", expand = FALSE) +
      labs(x = paste0("Chromosome ", chr_label, " position (Mb)"), y = NULL) +
      theme_minimal() +
      theme(panel.grid = element_blank(),
            axis.text.y = element_blank(),
            axis.title.x = element_text(size = 10, face = "bold"),
            axis.text.x = element_text(size = 9),
            panel.border = element_rect(color = "gray80", fill = NA, linewidth = 0.8),
            plot.margin = margin(0, 5, 5, 5))
  }
  
  # Row heights scale with the number of stacked gene lanes and annotation
  # tracks so labels stay readable when the locus is gene-dense.
  gene_weight <- max(1, min(n_lanes, 5))
  p_zoom <- build_zoom_panel()
  # When annotations are supplied we also draw a small zoom connector strip
  # right under the gene track -- two diagonal dashed lines from this CS's
  # lead_pos +/- zoom_window on the main LocusZoom's Mb axis to the outer
  # panel edges. Gives a visible "zoom-out" cue that points toward the
  # merged annotation box lower down in the module composite.
  has_link <- include_zoom_connector &&
              !is.null(annotations_tbl) && !is.null(annotation_tracks) &&
              length(annotation_tracks) > 0
  p_link <- if (has_link) build_zoom_connector() else NULL
  
  if (is.null(p_zoom)) {
    if (is.null(p_link)) {
      return(p_main / p_gene + plot_layout(heights = c(5, gene_weight)))
    }
    return(p_main / p_gene / p_link +
             plot_layout(heights = c(5, gene_weight, 0.35)))
  }
  zoom_h <- max(1.2, 0.6 * length(annotation_tracks) + 0.6)
  p_main / p_gene / p_link / p_zoom +
    plot_layout(heights = c(5, gene_weight, 0.35, zoom_h))
}

# plot_locuszoom_merged(): one LocusZoom panel with all credible sets overlaid,
# points coloured by CS, plus a bottom strip showing each CS genomic span on
# the same Mb axis (not a zoom window).
#
# Example:
#   plot_locuszoom_merged(df, NULL, "LZ merged", NULL, NULL, NULL, 5000,
#                         FALSE, FALSE, "18")
#   # -> patchwork of main + optional gene + CS span strip.
plot_locuszoom_merged <- function(df, locus_info, title, genes_df,
                                  annotations_tbl, annotation_tracks,
                                  zoom_window, include_zoom_panel,
                                  include_zoom_connector, chr_label,
                                  xlim_mode = c("context", "snp")) {
  xlim_mode <- match.arg(xlim_mode)
  data.table::setDT(df)
  data.table::setorderv(df, c("CS", "logp"), c(1, -1))
  df[, is_lead := (seq_len(.N) == 1L), by = CS]
  df[, cs_id := factor(as.character(CS), levels = sort(unique(as.character(CS))))]
  ucs <- levels(df$cs_id)
  nlev <- length(ucs)
  pal_named <- stats::setNames(
    grDevices::rainbow(nlev, s = 0.58, v = 0.88, start = 0.05, end = 0.92),
    ucs
  )
  has_focal <- !is.null(locus_info) &&
    !is.null(locus_info$gene_start) && !is.na(locus_info$gene_start) &&
    !is.null(locus_info$gene_end) && !is.na(locus_info$gene_end)
  has_many <- !is.null(genes_df) && nrow(genes_df) > 0
  x_min_bp <- min(df$POS, na.rm = TRUE)
  x_max_bp <- max(df$POS, na.rm = TRUE)
  if (identical(xlim_mode, "context")) {
    if (has_focal) {
      x_min_bp <- min(x_min_bp, as.numeric(locus_info$gene_start))
      x_max_bp <- max(x_max_bp, as.numeric(locus_info$gene_end))
    }
    if (has_many) {
      x_min_bp <- min(x_min_bp, min(genes_df$start, na.rm = TRUE))
      x_max_bp <- max(x_max_bp, max(genes_df$end, na.rm = TRUE))
    }
  }
  pad <- max((x_max_bp - x_min_bp) * 0.02, 1)
  x_lims_mb <- c((x_min_bp - pad) / 1e6, (x_max_bp + pad) / 1e6)
  y_max <- max(df$logp, na.rm = TRUE) * 1.12

  p_main <- ggplot2::ggplot(df[is_lead == FALSE],
                            ggplot2::aes(x = POS / 1e6, y = logp, color = cs_id)) +
    ggplot2::geom_point(size = 2.15, alpha = 0.86) +
    ggplot2::geom_point(
      data = df[is_lead == TRUE],
      mapping = ggplot2::aes(x = POS / 1e6, y = logp, fill = cs_id),
      inherit.aes = FALSE,
      shape = 23, color = "black", size = 3.3, stroke = 0.45) +
    ggplot2::scale_color_manual(values = pal_named, name = "Credible set", drop = FALSE) +
    ggplot2::scale_fill_manual(values = pal_named, guide = "none") +
    ggplot2::coord_cartesian(xlim = x_lims_mb, ylim = c(0, y_max), expand = FALSE) +
    ggplot2::labs(title = title, x = NULL,
                  y = expression(-log[10](italic(P)))) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(size = 10, face = "bold"),
      axis.title.y = ggplot2::element_text(size = 10, face = "bold",
                                           margin = ggplot2::margin(r = 2)),
      axis.text = ggplot2::element_text(size = 9),
      axis.text.x = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "gray80", fill = NA, linewidth = 0.8),
      plot.margin = ggplot2::margin(5, 5, 0, 2),
      legend.position = "right",
      legend.title = ggplot2::element_text(size = 8, face = "bold"),
      legend.text = ggplot2::element_text(size = 7))

  rng <- df[, .(xmin_mb = min(POS, na.rm = TRUE) / 1e6,
                xmax_mb = max(POS, na.rm = TRUE) / 1e6), by = .(CS = as.character(CS))]
  rng[, lab := vapply(CS, function(s) {
    if (nchar(s) > 44) paste0(substr(s, 1, 41), "...") else s
  }, character(1))]
  rng[, cs_f := factor(CS, levels = ucs)]

  p_strip <- ggplot2::ggplot(rng,
                             ggplot2::aes(xmin = xmin_mb, xmax = xmax_mb,
                                          ymin = 0, ymax = 1, fill = cs_f)) +
    ggplot2::geom_rect(color = "gray35", linewidth = 0.25) +
    ggplot2::geom_text(ggplot2::aes(x = (xmin_mb + xmax_mb) / 2, y = 0.5, label = lab),
                       size = 2.35, color = "black") +
    ggplot2::scale_fill_manual(values = pal_named, drop = FALSE, guide = "none") +
    ggplot2::coord_cartesian(xlim = x_lims_mb, ylim = c(0, 1), expand = FALSE) +
    ggplot2::labs(
      x = paste0("Chromosome ", chr_label, " position (Mb)"),
      y = NULL,
      subtitle = "Credible sets (genomic span; colours match points above)") +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.text.y = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(size = 9, face = "bold"),
      axis.text.x = ggplot2::element_text(size = 8),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "gray70", fill = NA, linewidth = 0.4),
      plot.subtitle = ggplot2::element_text(size = 7.5, face = "italic", hjust = 0),
      plot.margin = ggplot2::margin(0, 5, 4, 5))

  if (!has_focal && !has_many) {
    return(p_main / p_strip + patchwork::plot_layout(heights = c(5.2, 1.05)))
  }

  focal_sym <- if (!is.null(locus_info) && !is.null(locus_info$gene_symbol) &&
                    !is.na(locus_info$gene_symbol)) as.character(locus_info$gene_symbol) else NA_character_
  if (has_many) {
    track_df <- data.frame(
      start = as.numeric(genes_df$start), end = as.numeric(genes_df$end),
      strand = as.character(genes_df$strand), gene_name = as.character(genes_df$gene_name),
      stringsAsFactors = FALSE)
  } else {
    track_df <- data.frame(
      start = as.numeric(locus_info$gene_start), end = as.numeric(locus_info$gene_end),
      strand = if (!is.null(locus_info$gene_strand)) as.character(locus_info$gene_strand) else "+",
      gene_name = if (!is.na(focal_sym)) focal_sym else "gene",
      stringsAsFactors = FALSE)
  }
  track_df$is_focal <- !is.na(focal_sym) & track_df$gene_name == focal_sym
  track_df <- assign_gene_lanes(track_df)
  track_df$seg_x <- ifelse(track_df$strand == "-", track_df$end, track_df$start) / 1e6
  track_df$seg_xend <- ifelse(track_df$strand == "-", track_df$start, track_df$end) / 1e6
  track_df$mid_mb <- (track_df$start + track_df$end) / 2 / 1e6
  track_df$col <- ifelse(track_df$is_focal, "#C62828", "#1A237E")
  track_df <- clip_gene_segments_to_xlim(track_df, x_lims_mb)

  if (nrow(track_df) == 0L) {
    n_lanes <- 1L
    p_gene <- ggplot2::ggplot() +
      ggplot2::coord_cartesian(xlim = x_lims_mb, ylim = c(0.3, 1.2),
                               clip = "on", expand = FALSE) +
      ggplot2::labs(x = NULL, y = NULL) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(), axis.title.x = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(color = "gray80", fill = NA, linewidth = 0.8),
        plot.margin = ggplot2::margin(0, 5, 0, 5))
  } else {
    n_lanes <- max(track_df$lane, 1L)
    p_gene <- ggplot2::ggplot(track_df) +
      ggplot2::geom_segment(ggplot2::aes(x = seg_x_clip, xend = seg_xend_clip, y = lane, yend = lane,
                                         color = I(col)),
                            linewidth = 1.1,
                            arrow = grid::arrow(ends = "last", type = "closed",
                                                length = grid::unit(0.07, "inches"))) +
      ggrepel::geom_text_repel(ggplot2::aes(x = mid_mb_clip, y = lane + 0.35, label = gene_name,
                                             color = I(col),
                                             fontface = ifelse(is_focal, "bold.italic", "italic")),
                               size = 2.8, direction = "y", nudge_y = 0.25,
                               box.padding = 0.2, point.padding = 0.1,
                               segment.size = 0.45, segment.alpha = 0.75,
                               min.segment.length = 0, max.overlaps = Inf, seed = 42,
                               xlim = c(x_lims_mb[1], x_lims_mb[2]),
                               ylim = c(0.35, n_lanes + 1.55)) +
      ggplot2::scale_y_continuous(limits = c(0.1, n_lanes + 1.7), breaks = NULL) +
      ggplot2::coord_cartesian(xlim = x_lims_mb, ylim = c(0.1, n_lanes + 1.7),
                               clip = "on", expand = FALSE) +
      ggplot2::labs(x = NULL, y = NULL) +
      ggplot2::theme_minimal() +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(), axis.text.y = ggplot2::element_blank(),
        axis.text.x = ggplot2::element_blank(), axis.title.x = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(color = "gray80", fill = NA, linewidth = 0.8),
        plot.margin = ggplot2::margin(0, 5, 0, 5))
  }

  gene_weight <- max(1, min(n_lanes, 5))
  p_main / p_gene / p_strip +
    patchwork::plot_layout(heights = c(5.2, gene_weight, 1.05))
}
