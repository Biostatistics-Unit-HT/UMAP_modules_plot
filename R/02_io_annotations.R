# Regulatory annotation (BED / GFF) loader + window filter + track-kind probe.

# load_annotation_file(): read a regulatory annotation file into a tidy
# data.table with columns chrom, start, end, feature_type, score. Accepts
# either:
#   * BED-style (no header):
#       4 cols: chrom, start, end, feature_type
#       5+ cols: chrom, start, end, feature_type, score [, ...ignored]
#     A numeric 5th column marks the feature_type as continuous and is later
#     rendered as a per-position profile in the zoom panel. Empty/blank
#     score values are treated as NA (and the feature_type as discrete when
#     every row has NA).
#   * Full 9-column GFF/GTF (comment lines starting with # are skipped).
# Chromosome names are normalised by stripping any leading "chr".
#
# Example:
#   load_annotation_file("promoter.bed")  # 5 cols with numeric 5th
#   # -> data.table chrom=18, start=..., end=..., feature_type="promoter",
#   #    score=0.73
load_annotation_file <- function(path) {
  if (!file.exists(path)) stop(sprintf("Annotation file not found: %s", path))
  con <- file(path, "r"); on.exit(close(con), add = TRUE)
  first_data_line <- NULL
  while (length(line <- readLines(con, n = 1)) > 0) {
    if (!grepl("^#", line) && nchar(line) > 0) { first_data_line <- line; break }
  }
  if (is.null(first_data_line)) return(NULL)
  n_cols <- length(strsplit(first_data_line, "\t", fixed = TRUE)[[1]])
  
  if (n_cols >= 9) {
    raw <- fread(path, sep = "\t", header = FALSE, skip = "#",
                 col.names = c("chrom", "source", "feature_type", "start", "end",
                               "score", "strand", "frame", "attrs"),
                 select = 1:9, fill = TRUE)
    out <- raw[, .(chrom = as.character(chrom),
                   start = as.integer(start),
                   end   = as.integer(end),
                   feature_type = as.character(feature_type),
                   score = suppressWarnings(as.numeric(score)))]
  } else if (n_cols >= 4) {
    # BED: 4 required columns plus an optional numeric 5th `score` column.
    # We peek at the file's actual column count via a dry read so trailing
    # tabs or variable-length rows don't trip fread.
    peek <- fread(path, sep = "\t", header = FALSE, fill = TRUE,
                  nrows = 5)
    actual_cols <- ncol(peek)
    if (actual_cols >= 5) {
      raw <- fread(path, sep = "\t", header = FALSE, fill = TRUE,
                   col.names = c("chrom", "start", "end", "feature_type", "score"),
                   select = 1:5)
      out <- raw[, .(chrom = as.character(chrom),
                     start = as.integer(start),
                     end   = as.integer(end),
                     feature_type = as.character(feature_type),
                     score = suppressWarnings(as.numeric(score)))]
    } else {
      raw <- fread(path, sep = "\t", header = FALSE, fill = TRUE,
                   col.names = c("chrom", "start", "end", "feature_type"),
                   select = 1:4)
      out <- raw[, .(chrom = as.character(chrom),
                     start = as.integer(start),
                     end   = as.integer(end),
                     feature_type = as.character(feature_type),
                     score = NA_real_)]
    }
  } else {
    stop(sprintf("Unsupported annotation format (%d columns) in %s", n_cols, path))
  }
  out[, chrom := sub("^chr", "", chrom, ignore.case = TRUE)]
  out
}

# is_continuous_track(): TRUE when any row of `tbl` with this feature_type
# has a non-NA score, i.e. the lane should be rendered as a continuous
# profile rather than a block.
#
# Example:
#   is_continuous_track("atac_signal", annotations_tbl) # -> TRUE when scored
is_continuous_track <- function(ft, tbl) {
  if (is.null(tbl) || nrow(tbl) == 0) return(FALSE)
  if (!"score" %in% names(tbl)) return(FALSE)
  any(tbl$feature_type == ft & !is.na(tbl$score))
}

# filter_annotations_window(): subset an annotation data.table to features
# overlapping [win_start, win_end] on chromosome chr.
#
# Example:
#   filter_annotations_window(anno_tbl, "chr18", 80000273, 80010273)
#   # -> data.table containing the promoter_flanking_region at
#   #    80009802-80011199.
filter_annotations_window <- function(tbl, chr, win_start, win_end) {
  if (is.null(tbl) || nrow(tbl) == 0) return(tbl[0])
  chr_no <- sub("^chr", "", chr, ignore.case = TRUE)
  tbl[chrom == chr_no & start <= win_end & end >= win_start]
}
