# QTL Z vs disease Z scatter for a single (module, CS).

# plot_zscore(): QTL vs disease Z-score scatter for a single module.
#
# Example:
#   plot_zscore("icd10_j15_zscore_merge_notnull.csv", "Coloc Z-Scores\n[M_18448]")
#   # -> ggplot showing points with diagonals and dashed zero lines.
plot_zscore <- function(z_file, title = "Coloc Z-Scores") {
  if (is.na(z_file) || z_file == "NA" || !file.exists(z_file)) {
    if (!is.na(z_file) && z_file != "NA") cat(sprintf("Warning: Z-score file not found at %s\n", z_file))
    return(plot_spacer())
  }
  df <- fread(z_file)
  qtl_col <- names(df)[grepl("z_qtl", names(df), ignore.case = TRUE)][1]
  if (is.na(qtl_col)) qtl_col <- names(df)[grepl("qtl", names(df), ignore.case = TRUE)][1]
  dis_col <- names(df)[grepl("z_icd10|z_dis", names(df), ignore.case = TRUE)][1]
  if (is.na(dis_col)) dis_col <- names(df)[!names(df) %in% c("snp", qtl_col)][1]
  if (is.na(qtl_col) || is.na(dis_col)) {
    cat(sprintf("Warning: Could not identify Z-score columns in %s\n", z_file))
    return(plot_spacer())
  }
  max_val   <- max(abs(c(df[[qtl_col]], df[[dis_col]])), na.rm = TRUE)
  limit_val <- max_val * 1.1
  ggplot(df, aes(x = .data[[qtl_col]], y = .data[[dis_col]])) +
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
}
