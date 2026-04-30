# QTL Z vs disease Z scatter for a single (module, CS).

# plot_zscore(): QTL vs disease Z-score scatter for a single module. Accepts the
# classic `snp, z_disease, z_qtl` table or the extended form with a `cs_qtl`
# credible-set id column (must not be mistaken for the QTL Z column).
#
# Example:
#   plot_zscore("icd10_j15_zscore_merge_notnull.csv", "Coloc Z-Scores\n[M_18448]")
#   # -> ggplot with diagonals and dashed zero lines.
#   plot_zscore("z_ext.csv", "Coloc Z-Scores\n[M_1]")  # z_ext.csv has snp,z_disease,cs_qtl,z_qtl
#   # -> same scatter; optional caption shows one unique cs_qtl when present.
plot_zscore <- function(z_file, title = "Coloc Z-Scores") {
  if (is.na(z_file) || z_file == "NA" || !file.exists(z_file)) {
    if (!is.na(z_file) && z_file != "NA") cat(sprintf("Warning: Z-score file not found at %s\n", z_file))
    return(plot_spacer())
  }
  df <- fread(z_file)
  nms <- names(df)
  # Prefer exact names from the extract_z_lz / coloc pipeline (avoid matching cs_qtl for "qtl").
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
  if (is.na(qtl_col) || is.na(dis_col) || !nzchar(qtl_col) || !nzchar(dis_col)) {
    cat(sprintf("Warning: Could not identify Z-score columns in %s\n", z_file))
    return(plot_spacer())
  }
  max_val   <- max(abs(c(df[[qtl_col]], df[[dis_col]])), na.rm = TRUE)
  limit_val <- max_val * 1.1
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
          axis.title  = element_text(size = 10, face = "bold"),
          axis.text   = element_text(size = 9),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "gray90", fill = NA, linewidth = 1))
  if ("cs_qtl" %in% nms) {
    cs_u <- unique(stats::na.omit(as.character(df[["cs_qtl"]])))
    cs_u <- cs_u[nzchar(cs_u)]
    if (length(cs_u) >= 1L) {
      cap <- if (length(cs_u) == 1L) cs_u[[1]] else paste(length(cs_u), "credible sets")
      if (nchar(cap) > 120) cap <- paste0(substr(cap, 1, 117), "...")
      p <- p + labs(caption = cap) +
        theme(plot.caption = element_text(size = 6, hjust = 0))
    }
  }
  p
}
