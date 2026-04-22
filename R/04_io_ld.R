# LD matrix loader + SNP matching helpers + bin palette for the LocusZoom.

# load_ld_matrix(): read a plink2 "--export A" .raw genotype file and return
# the r^2 correlation matrix among its variants. Adapted from make_ld.R:
#   * starts straight from the .raw file (no plink2 / dentist / pfile step),
#   * removes constant-genotype columns,
#   * parses SNP ids of the form "<chr>:<pos>:<ref>:<alt>_<ea>(/<oa>)",
#   * drops duplicate SNP ids,
#   * flips dosages so the effect allele is the one that sorts first
#     alphabetically (consistent orientation across SNPs),
#   * mean-imputes missing dosages,
#   * computes crossprod(scale(geno)) / (n-1) (== cor(geno) but faster),
#   * squares the result to get r^2.
# Row/column names are the short SNP ids (e.g. "18:79869986:A:G").
#
# Returns NULL if the file has fewer than 2 variable SNPs.
#
# Example:
#   ld <- load_ld_matrix("UMAP_lz_values.raw")
#   ld["18:80005273:C:T", "18:79869986:A:G"]  # r^2 of the two variants
load_ld_matrix <- function(raw_path) {
  if (is.na(raw_path) || raw_path == "NA" || !file.exists(raw_path)) {
    if (!is.na(raw_path) && raw_path != "NA")
      cat(sprintf("Warning: LD .raw file not found at %s\n", raw_path))
    return(NULL)
  }
  geno <- fread(raw_path)[, -c(1:6)]
  not_same <- which(vapply(geno, function(x) length(unique(x)) > 1, logical(1)))
  if (length(not_same) < 2) {
    cat(sprintf("Warning: %s has fewer than 2 variable SNPs; skipping LD.\n", raw_path))
    return(NULL)
  }
  geno <- geno[, ..not_same]
  
  snp_info <- strsplit(colnames(geno), "_|\\(/|\\)")
  snp_info <- Reduce(rbind, snp_info) |> as.data.frame()
  if (ncol(snp_info) == 1) snp_info <- as.data.frame(t(snp_info))
  colnames(snp_info) <- c("SNP", "ea", "oa")[seq_len(ncol(snp_info))]
  rownames(snp_info) <- NULL
  colnames(geno) <- snp_info$SNP
  
  dup <- which(duplicated(snp_info$SNP) |
               duplicated(snp_info$SNP, fromLast = TRUE))
  if (length(dup) > 0) {
    snp_info <- snp_info[-dup, , drop = FALSE]
    geno <- geno[, -..dup]
  }
  
  if (all(c("ea", "oa") %in% colnames(snp_info))) {
    idx_flip <- which(!(snp_info$ea <= snp_info$oa))
    if (length(idx_flip) > 0) {
      switch_0_2 <- function(x) (x * -1) + 2
      geno[, (idx_flip) := lapply(.SD, switch_0_2), .SDcols = idx_flip]
    }
  }
  
  geno <- apply(geno, 2, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE); x
  })
  
  X_scaled <- scale(geno)
  ld <- crossprod(X_scaled) / (nrow(geno) - 1)
  ld^2
}

# ld_pos_key(): turn any SNP identifier into a "<chr_no_chr_prefix>:<pos>"
# key for matching. Works on both LD matrix names like "18:79869986:A:G"
# and CS strings like "chr18:80005273:C:T".
#
# Example:
#   ld_pos_key("chr18:80005273:C:T")  # -> "18:80005273"
#   ld_pos_key("18:79869986:A:G")     # -> "18:79869986"
ld_pos_key <- function(snp) {
  if (is.na(snp)) return(NA_character_)
  s <- sub("^chr", "", as.character(snp), ignore.case = TRUE)
  parts <- strsplit(s, ":", fixed = TRUE)[[1]]
  if (length(parts) < 2) return(NA_character_)
  paste(parts[1], parts[2], sep = ":")
}

# ld_vec_for_lead(): extract the r^2 vector of a given lead SNP vs every
# other SNP in the matrix, keyed by "<chr>:<pos>" so it can be joined onto
# LocusZoom data rows. Returns NULL when the lead is absent.
#
# Example:
#   ld_vec_for_lead(ld, "chr18:80005273:C:T")
#   # -> named numeric vector c("18:79869986" = 0.023, "18:79870170" = 0.001, ...)
ld_vec_for_lead <- function(ld_mat, lead_id) {
  if (is.null(ld_mat) || is.null(lead_id) || is.na(lead_id)) return(NULL)
  lead_key <- ld_pos_key(lead_id)
  mat_keys <- vapply(rownames(ld_mat), ld_pos_key, character(1))
  hit <- which(mat_keys == lead_key)
  if (length(hit) == 0) {
    cat(sprintf("Warning: lead SNP %s (key %s) not in LD matrix; skipping LD colouring for this CS.\n",
                lead_id, lead_key))
    return(NULL)
  }
  v <- ld_mat[hit[1], ]
  setNames(as.numeric(v), mat_keys)
}

# bin_ld_r2(): convert a numeric r^2 vector (or NAs) to a factor with the
# classic 5-level LocusZoom palette bins. The factor levels are fixed so
# scale_fill_manual() doesn't drop an unused colour.
#
# Example:
#   bin_ld_r2(c(0.15, 0.5, 0.92, NA))
#   # -> factor with levels "<0.2" "0.2-0.4" "0.4-0.6" "0.6-0.8" ">=0.8" "NA"
bin_ld_r2 <- function(x) {
  lev <- c("<0.2", "0.2-0.4", "0.4-0.6", "0.6-0.8", ">=0.8", "NA")
  lab <- character(length(x))
  for (i in seq_along(x)) {
    v <- x[i]
    lab[i] <- if (is.na(v)) "NA"
              else if (v < 0.2) "<0.2"
              else if (v < 0.4) "0.2-0.4"
              else if (v < 0.6) "0.4-0.6"
              else if (v < 0.8) "0.6-0.8"
              else ">=0.8"
  }
  factor(lab, levels = lev)
}

LD_COLORS <- c("<0.2"    = "#1F3B8B",  # navy
               "0.2-0.4" = "#71acf5",  # light blue
               "0.4-0.6" = "#38d940",  # green
               "0.6-0.8" = "#f5e618",  # orange
               ">=0.8"   = "#f53333",  # red
               "NA"      = "#B0B0B0")  # gray for SNPs absent from matrix
