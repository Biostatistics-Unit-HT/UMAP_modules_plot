# QTL beta vs disease beta scatter for a single (module, CS).

# resolve_beta_columns(): map data.frame column names to QTL-beta and disease-beta
# fields; avoids treating cs_qtl as a QTL beta column. Legacy z_qtl / z_disease
# columns are still accepted so old CSVs keep working.
#
# Example:
#   resolve_beta_columns(c("snp", "beta_disease", "cs_qtl", "beta_qtl"))
#   # -> list(qtl_col = "beta_qtl", dis_col = "beta_disease")
resolve_beta_columns <- function(nms) {
  qtl_col <- if ("beta_qtl" %in% nms) {
    "beta_qtl"
  } else if ("z_qtl" %in% nms) {
    "z_qtl"
  } else {
    hit <- nms[grepl("beta_qtl|z_qtl", nms, ignore.case = TRUE)][1]
    if (!is.na(hit) && nzchar(hit)) hit else {
      qhit <- nms[grepl("qtl", nms, ignore.case = TRUE) &
                  !nms %in% c("cs_qtl") & !grepl("^cs_", nms, ignore.case = TRUE)][1]
      if (!is.na(qhit) && nzchar(qhit)) qhit else NA_character_
    }
  }
  dis_col <- if ("beta_disease" %in% nms) {
    "beta_disease"
  } else if ("z_disease" %in% nms) {
    "z_disease"
  } else {
    hit <- nms[grepl("beta_dis|z_icd10|z_dis", nms, ignore.case = TRUE)][1]
    if (!is.na(hit) && nzchar(hit)) hit else {
      skip <- unique(c("snp", qtl_col, "cs_qtl", nms[grepl("^cs_", nms, ignore.case = TRUE)]))
      skip <- skip[!is.na(skip)]
      cand <- setdiff(nms, skip)
      if (length(cand)) cand[[1]] else NA_character_
    }
  }
  list(qtl_col = qtl_col, dis_col = dis_col)
}

# shorten_cs_label(): compact label for a cs_qtl token (last :: segment or trim).
#
# Example:
#   shorten_cs_label("chr18::GH::T_CD4_CM::ENSG1::chr18:80005273:C:T::L1")
#   # -> "L1" or similar tail token.
shorten_cs_label <- function(s, max_len = 56L) {
  s <- as.character(s)
  if (!nzchar(s)) return("")
  parts <- strsplit(s, "::", fixed = TRUE)[[1]]
  tail <- if (length(parts)) parts[[length(parts)]] else s
  if (nchar(tail) > max_len) paste0(substr(tail, 1, max_len - 3L), "...") else tail
}

# beta_axis_limits(): symmetric limits for QTL-beta (x) and disease-beta (y).
# When linked_axes is TRUE both axes share max(abs(all beta)); when FALSE each
# axis uses only its own column (independent scale).
#
# Example:
#   beta_axis_limits(c(-3, 5), c(-1, 2), linked_axes = TRUE)
#   # -> list(xlim = c(-5.5, 5.5), ylim = c(-5.5, 5.5))
#   beta_axis_limits(c(-3, 5), c(-1, 2), linked_axes = FALSE)
#   # -> list(xlim = c(-5.5, 5.5), ylim = c(-2.2, 2.2))
beta_axis_limits <- function(x_vals, y_vals, linked_axes = TRUE, pad = 1.1) {
  lim_from <- function(v) {
    m <- max(abs(v), na.rm = TRUE)
    if (!is.finite(m) || m <= 0) m <- 1
    m * pad
  }
  if (isTRUE(linked_axes)) {
    L <- lim_from(c(x_vals, y_vals))
    list(xlim = c(-L, L), ylim = c(-L, L))
  } else {
    list(xlim = c(-lim_from(x_vals), lim_from(x_vals)),
         ylim = c(-lim_from(y_vals), lim_from(y_vals)))
  }
}

# plot_beta_scatter_df(): scatter QTL beta vs disease beta from an in-memory table
# (uniform blue points, no colour scale). Optional panel_cs_label becomes subtitle.
#
# Example:
#   plot_beta_scatter_df(data.frame(beta_qtl = 1:3, beta_disease = 3:1), "Beta | cell | g\n[M]")
#   # -> ggplot with diagonals and dashed zero lines (linked axes).
#   plot_beta_scatter_df(..., linked_axes = FALSE)
#   # -> x and y limits from each column separately.
plot_beta_scatter_df <- function(df, title = "Coloc Betas", panel_cs_label = NULL,
                                 linked_axes = TRUE) {
  if (is.null(df) || nrow(df) == 0L) return(plot_spacer())
  nms <- names(df)
  cols <- resolve_beta_columns(nms)
  qtl_col <- cols$qtl_col
  dis_col <- cols$dis_col
  if (is.na(qtl_col) || is.na(dis_col) || !nzchar(qtl_col) || !nzchar(dis_col)) {
    cat("Warning: Could not identify beta columns in supplied data.frame\n")
    return(plot_spacer())
  }
  lims <- beta_axis_limits(df[[qtl_col]], df[[dis_col]], linked_axes = linked_axes)
  p <- ggplot(df, aes(x = .data[[qtl_col]], y = .data[[dis_col]])) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = 0, color = "black", linetype = "dashed", alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, color = "gray80", linetype = "dotted") +
    geom_abline(slope = -1, intercept = 0, color = "gray80", linetype = "dotted") +
    geom_point(fill = "#1E90FF", color = "white", shape = 21,
               size = 3.5, alpha = 0.8, stroke = 0.3)
  p <- if (isTRUE(linked_axes)) {
    p + coord_fixed(xlim = lims$xlim, ylim = lims$ylim)
  } else {
    p + coord_cartesian(xlim = lims$xlim, ylim = lims$ylim)
  }
  p <- p +
    theme_minimal() +
    labs(title = title, x = "QTL Beta", y = "Disease Beta") +
    theme(plot.title  = element_text(size = 10, face = "bold"),
          plot.subtitle = element_text(size = 7, hjust = 0),
          axis.title  = element_text(size = 10, face = "bold"),
          axis.text   = element_text(size = 9),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "gray90", fill = NA, linewidth = 1))
  if (!is.null(panel_cs_label) && nzchar(as.character(panel_cs_label))) {
    lab <- shorten_cs_label(panel_cs_label, max_len = 72L)
    p <- p + labs(subtitle = lab)
  } else if ("cs_qtl" %in% nms) {
    cs_u <- unique(stats::na.omit(as.character(df[["cs_qtl"]])))
    cs_u <- cs_u[nzchar(cs_u)]
    if (length(cs_u) == 1L) {
      cap <- shorten_cs_label(cs_u[[1]], max_len = 120L)
      p <- p + labs(caption = cap) +
        theme(plot.caption = element_text(size = 6, hjust = 0))
    }
  }
  p
}

# plot_beta_scatter(): read a beta CSV and draw one scatter (see plot_beta_scatter_df()).
#
# Example:
#   plot_beta_scatter("icd10_j15_beta_merge_notnull.csv", "Coloc Betas\n[M_18448]")
#   # -> ggplot with diagonals and dashed zero lines.
plot_beta_scatter <- function(beta_file, title = "Coloc Betas", linked_axes = TRUE) {
  if (is.na(beta_file) || beta_file == "NA" || !file.exists(beta_file)) {
    if (!is.na(beta_file) && beta_file != "NA") cat(sprintf("Warning: beta file not found at %s\n", beta_file))
    return(plot_spacer())
  }
  df <- fread(beta_file)
  plot_beta_scatter_df(df, title, panel_cs_label = NULL, linked_axes = linked_axes)
}

# build_beta_column(): one or more beta panels from the module beta table — split by
# cs_qtl when multiple values exist; if this_cs matches a cs_qtl row, only that
# subset is drawn for this grid row.
#
# Example:
#   z <- data.table::fread(text = "snp,beta_disease,cs_qtl,beta_qtl\na,1,c1,2\nb,-1,c2,-2")
#   build_beta_column(z, "Coloc Beta | cell | g\n[M]", this_cs = NA_character_)
#   # -> patchwork column with two stacked scatters (c1 and c2).
build_beta_column <- function(z_tbl, title, this_cs = NA_character_,
                              linked_axes = TRUE) {
  if (is.null(z_tbl) || nrow(z_tbl) == 0L) return(plot_spacer())
  if (!"cs_qtl" %in% names(z_tbl)) return(plot_beta_scatter_df(z_tbl, title, linked_axes = linked_axes))
  keys <- unique(stats::na.omit(as.character(z_tbl[["cs_qtl"]])))
  keys <- keys[nzchar(keys)]
  if (length(keys) == 0L) return(plot_beta_scatter_df(z_tbl, title, linked_axes = linked_axes))
  this_cs <- as.character(this_cs)
  matched <- (!is.na(this_cs) && length(this_cs) == 1L && nzchar(this_cs) &&
              this_cs %in% keys)
  if (matched) {
    sub <- dplyr::filter(z_tbl, as.character(.data[["cs_qtl"]]) == this_cs)
    return(plot_beta_scatter_df(sub, title, panel_cs_label = this_cs,
                                linked_axes = linked_axes))
  }
  if (length(keys) == 1L) {
    k1 <- keys[[1]]
    sub <- dplyr::filter(z_tbl, as.character(.data[["cs_qtl"]]) == k1)
    return(plot_beta_scatter_df(sub, title, panel_cs_label = k1,
                                linked_axes = linked_axes))
  }
  pieces <- lapply(keys, function(k) {
    sub <- dplyr::filter(z_tbl, as.character(.data[["cs_qtl"]]) == k)
    plot_beta_scatter_df(sub, title, panel_cs_label = k, linked_axes = linked_axes)
  })
  wrap_plots(pieces, ncol = 1L) +
    patchwork::plot_layout(heights = rep(1, length(pieces)))
}
