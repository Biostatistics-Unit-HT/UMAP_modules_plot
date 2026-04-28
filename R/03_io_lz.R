# LocusZoom CSV loader + CS-string parser.

# load_lz_file(): read a LocusZoom CSV and normalise its columns. Supports:
#   * 5 columns with header -- CHR, CELL, GENE, POS, P
#   * 6 columns without header -- CS, CHR, CELL, GENE, POS, P
#     (CS encodes the lead SNP, e.g.
#      "chr18::GH_meta_cardinal::T_CD4_CM:ENSG00000101546::chr18:80005273:C:T::L1")
# Returns a data.table with those canonical column names, or NULL when the
# file is missing or has an unexpected shape.
#
# Example:
#   load_lz_file("UMAP_lz.csv")
#   # -> data.table with columns CS, CHR, CELL, GENE, POS, P (1,055 rows)
load_lz_file <- function(path) {
  if (is.na(path) || path == "NA" || !file.exists(path)) {
    if (!is.na(path) && path != "NA") cat(sprintf("Warning: LocusZoom file not found at %s\n", path))
    return(NULL)
  }
  first_line <- readLines(path, n = 1, warn = FALSE)
  has_header <- length(first_line) > 0 &&
                grepl("\\bCHR\\b", first_line) &&
                grepl("\\bPOS\\b", first_line) &&
                grepl("\\bP\\b",   first_line)
  df <- if (has_header) fread(path) else fread(path, header = FALSE)
  if (!has_header) {
    ncol_df <- ncol(df)
    if (ncol_df == 6)      setnames(df, c("CS", "CHR", "CELL", "GENE", "POS", "P"))
    else if (ncol_df == 5) setnames(df, c("CHR", "CELL", "GENE", "POS", "P"))
    else {
      cat(sprintf("Warning: %s has an unexpected number of columns (%d)\n",
                  path, ncol_df))
      return(NULL)
    }
  }
  df
}

# parse_cs(): split a credible-set identifier like
# "chr18::GH_meta_cardinal::T_CD4_CM:ENSG00000101546::chr18:80005273:C:T::L1"
# into its components. Returns a named list with cs, cell, eGene (no version),
# chrom (no "chr"), lead_pos (numeric, bp). Any field may be NA when the
# string doesn't match the expected pattern.
#
# Example:
#   parse_cs("chr18::GH_meta_cardinal::T_CD4_CM:ENSG00000101546::chr18:80005273:C:T::L1")
#   # -> list(cs=..., cell="T_CD4_CM", eGene="ENSG00000101546",
#   #         chrom="18", lead_pos=80005273)
parse_cs <- function(cs) {
  out <- list(cs = cs, cell = NA_character_, eGene = NA_character_,
              chrom = NA_character_, lead_pos = NA_real_)
  if (is.na(cs) || !nzchar(cs)) return(out)
  parts <- strsplit(as.character(cs), "::", fixed = TRUE)[[1]]
  # Cell / gene come from e.g. "T_CD4_CM:ENSG00000101546" (3rd :: segment).
  if (length(parts) >= 3) {
    cg <- strsplit(parts[3], ":", fixed = TRUE)[[1]]
    if (length(cg) >= 2) {
      out$cell  <- cg[1]
      out$eGene <- sub("\\.\\d+$", "", cg[2])
    }
  }
  # Lead SNP is any chr<N>:<pos>:<ALLELE>:<ALLELE> token.
  tok <- regmatches(cs, regexpr("chr[0-9XYM]+:\\d+:[ACGTN]+:[ACGTN]+",
                                cs, ignore.case = TRUE))
  if (length(tok) == 1 && nzchar(tok)) {
    sn <- strsplit(tok, ":", fixed = TRUE)[[1]]
    out$chrom    <- sub("^chr", "", sn[1], ignore.case = TRUE)
    out$lead_pos <- suppressWarnings(as.numeric(sn[2]))
  }
  out
}
