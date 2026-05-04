# QTL Z vs disease Z scatter for a single (module, CS).

# resolve_zscore_columns(): map data.frame column names to QTL-Z and disease-Z
# fields; avoids treating cs_qtl as a QTL Z column.
#
# Example:
#   resolve_zscore_columns(c("snp", "z_disease", "cs_qtl", "z_qtl"))
#   # -> list(qtl_col = "z_qtl", dis_col = "z_disease")
resolve_zscore_columns <- function(nms) {
  qtl_col <- if ("z_qtl" %in% nms) {
    "z_qtl"
  } else {
    hit <- nms[grepl("z_qtl", nms, ignore.case = TRUE)][1]
    if (!is.na(hit) && nzchar(hit)) hit else {
      qhit <- nms[grepl("qtl", nms, ignore.case = TRUE) &
                  !nms %in% c("cs_qtl") & !grepl("^cs_", nms, ignore.case = TRUE)][1]
      if (!is.na(qhit) && nzchar(qhit)) qhit else NA_character_
    }
  }
  dis_col <- if ("z_disease" %in% nms) {
    "z_disease"
  } else {
    hit <- nms[grepl("z_icd10|z_dis", nms, ignore.case = TRUE)][1]
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

# plot_zscore_df(): scatter QTL Z vs disease Z from an in-memory table (uniform
# blue points, no Z-based colour scale). Optional panel_cs_label becomes subtitle.
#
# Example:
#   plot_zscore_df(data.frame(z_qtl = 1:3, z_disease = 3:1), "Z | cell | g\n[M]")
#   # -> ggplot with diagonals and dashed zero lines.
plot_zscore_df <- function(df, title = "Coloc Z-Scores", panel_cs_label = NULL) {
  if (is.null(df) || nrow(df) == 0L) return(plot_spacer())
  nms <- names(df)
  cols <- resolve_zscore_columns(nms)
  qtl_col <- cols$qtl_col
  dis_col <- cols$dis_col
  if (is.na(qtl_col) || is.na(dis_col) || !nzchar(qtl_col) || !nzchar(dis_col)) {
    cat("Warning: Could not identify Z-score columns in supplied data.frame\n")
    return(plot_spacer())
  }
  max_val   <- max(abs(c(df[[qtl_col]], df[[dis_col]])), na.rm = TRUE)
  limit_val <- max_val * 1.1
  if (!is.finite(limit_val) || limit_val <= 0) limit_val <- 1
  p <- ggplot(df, aes(x = .data[[qtl_col]], y = .data[[dis_col]])) +
    geom_hline(yintercept = 0, color = "black", linetype = "dashed", alpha = 0.5) +
    geom_vline(xintercept = 0, color = "black", linetype = "dashed", alpha = 0.5) +
    geom_abline(slope = 1, intercept = 0, color = "gray80", linetype = "dotted") +
    geom_abline(slope = -1, intercept = 0, color = "gray80", linetype = "dotted") +
    geom_point(fill = "#1E90FF", color = "white", shape = 21,
               size = 2.5, alpha = 0.8, stroke = 0.3) +
    coord_fixed(xlim = c(-limit_val, limit_val), ylim = c(-limit_val, limit_val)) +
    theme_minimal() +
    labs(title = title, x = "QTL Z-Score", y = "Disease Z-Score") +
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

# plot_zscore(): read a Z-score CSV and draw one scatter (see plot_zscore_df()).
#
# Example:
#   plot_zscore("icd10_j15_zscore_merge_notnull.csv", "Coloc Z-Scores\n[M_18448]")
#   # -> ggplot with diagonals and dashed zero lines.
plot_zscore <- function(z_file, title = "Coloc Z-Scores") {
  if (is.na(z_file) || z_file == "NA" || !file.exists(z_file)) {
    if (!is.na(z_file) && z_file != "NA") cat(sprintf("Warning: Z-score file not found at %s\n", z_file))
    return(plot_spacer())
  }
  df <- fread(z_file)
  plot_zscore_df(df, title, panel_cs_label = NULL)
}

# build_zscore_column(): one or more Z panels from the module Z table — split by
# cs_qtl when multiple values exist; if this_cs matches a cs_qtl row, only that
# subset is drawn for this grid row.
#
# Example:
#   z <- data.table::fread(text = "snp,z_disease,cs_qtl,z_qtl\na,1,c1,2\nb,-1,c2,-2")
#   build_zscore_column(z, "Coloc Z | cell | g\n[M]", this_cs = NA_character_)
#   # -> patchwork column with two stacked scatters (c1 and c2).
build_zscore_column <- function(z_tbl, title, this_cs = NA_character_) {
  if (is.null(z_tbl) || nrow(z_tbl) == 0L) return(plot_spacer())
  if (!"cs_qtl" %in% names(z_tbl)) return(plot_zscore_df(z_tbl, title))
  keys <- unique(stats::na.omit(as.character(z_tbl[["cs_qtl"]])))
  keys <- keys[nzchar(keys)]
  if (length(keys) == 0L) return(plot_zscore_df(z_tbl, title))
  this_cs <- as.character(this_cs)
  matched <- (!is.na(this_cs) && length(this_cs) == 1L && nzchar(this_cs) &&
              this_cs %in% keys)
  if (matched) {
    sub <- dplyr::filter(z_tbl, as.character(.data[["cs_qtl"]]) == this_cs)
    return(plot_zscore_df(sub, title, panel_cs_label = this_cs))
  }
  if (length(keys) == 1L) {
    k1 <- keys[[1]]
    sub <- dplyr::filter(z_tbl, as.character(.data[["cs_qtl"]]) == k1)
    return(plot_zscore_df(sub, title, panel_cs_label = k1))
  }
  pieces <- lapply(keys, function(k) {
    sub <- dplyr::filter(z_tbl, as.character(.data[["cs_qtl"]]) == k)
    plot_zscore_df(sub, title, panel_cs_label = k)
  })
  wrap_plots(pieces, ncol = 1L) +
    patchwork::plot_layout(heights = rep(1, length(pieces)))
}
