# GTF / gene-cache I/O + region queries for the LocusZoom gene track.

# load_gene_table(): return a data.table of protein-coding gene bodies ready
# to use for the LocusZoom gene track. Accepts either:
#   (a) a pre-built cache TSV with columns chrom, start, end, strand,
#       gene_name, gene_id -- loaded directly; or
#   (b) a raw GENCODE GTF -- in which case a cache TSV is built next to the
#       GTF (re-built only when the GTF is newer than the cache).
#
# Cache header is detected by reading the first non-empty line of the file.
#
# Example:
#   load_gene_table("gencode.v49.annotation.gtf")
#   # -> builds "gencode.v49.annotation.gtf.coding_genes.tsv" (~20k rows)
#   load_gene_table("gencode.v49.annotation.gtf.coding_genes.tsv")
#   # -> loads the cached TSV directly, no awk step.
load_gene_table <- function(gtf_path) {
  if (!file.exists(gtf_path)) stop(sprintf("GTF / gene cache not found: %s", gtf_path))
  expected_header <- c("chrom", "start", "end", "strand", "gene_name", "gene_id")
  first_line <- readLines(gtf_path, n = 1, warn = FALSE)
  if (length(first_line) > 0) {
    tokens <- strsplit(first_line, "\t", fixed = TRUE)[[1]]
    if (identical(tokens, expected_header)) {
      cat(sprintf("Loading pre-built gene cache %s ...\n", gtf_path))
      return(fread(gtf_path))
    }
  }
  cache_path <- paste0(gtf_path, ".coding_genes.tsv")
  if (file.exists(cache_path) &&
      file.info(cache_path)$mtime >= file.info(gtf_path)$mtime) {
    cat(sprintf("Loading existing gene cache %s ...\n", cache_path))
    return(fread(cache_path))
  }
  cat(sprintf("Building protein-coding gene cache from %s (one-off, ~1 min) ...\n", gtf_path))
  tmp <- tempfile(fileext = ".tsv")
  on.exit(unlink(tmp), add = TRUE)
  cmd <- sprintf("awk -F'\\t' '$0 !~ /^#/ && $3 == \"gene\"' %s > %s",
                 shQuote(gtf_path), shQuote(tmp))
  if (system(cmd) != 0) stop("Failed to pre-filter GTF with awk.")
  raw <- fread(tmp, sep = "\t", header = FALSE,
               col.names = c("chrom", "source", "feature", "start", "end",
                             "score", "strand", "frame", "attrs"))
  raw[, gene_type := sub('.*gene_type "([^"]+)".*', "\\1", attrs)]
  raw <- raw[gene_type == "protein_coding"]
  raw[, gene_name := sub('.*gene_name "([^"]+)".*', "\\1", attrs)]
  raw[, gene_id   := sub('.*gene_id "([^"]+)".*',   "\\1", attrs)]
  out <- raw[, .(chrom, start, end, strand, gene_name, gene_id)]
  fwrite(out, cache_path, sep = "\t")
  cat(sprintf("Cached %d protein-coding genes to %s\n", nrow(out), cache_path))
  out
}

# get_genes_in_region(): protein-coding genes overlapping [win_start,win_end]
# on chromosome chr. When `whitelist_ids` is non-NULL, rows are further
# restricted to those whose `gene_id` (version stripped, e.g.
# "ENSG00000186092.7" -> "ENSG00000186092") is in the whitelist.
#
# Example:
#   get_genes_in_region(gtf_tbl, "chr18", 79800000, 80200000,
#                       whitelist_ids = c("ENSG00000101546"))
#   # -> data.table with only RBFA.
get_genes_in_region <- function(gtf_tbl, chr, win_start, win_end,
                                whitelist_ids = NULL) {
  if (is.null(gtf_tbl) || nrow(gtf_tbl) == 0) return(NULL)
  chr_no <- sub("^chr", "", chr, ignore.case = TRUE)
  targets <- unique(c(chr, chr_no, paste0("chr", chr_no)))
  hit <- gtf_tbl[chrom %in% targets & start <= win_end & end >= win_start]
  if (nrow(hit) == 0) return(NULL)
  if (!is.null(whitelist_ids) && length(whitelist_ids) > 0) {
    bare <- sub("\\.\\d+$", "", hit$gene_id)
    hit <- hit[bare %in% whitelist_ids]
    if (nrow(hit) == 0) return(NULL)
  }
  hit
}

# assign_gene_lanes(): greedy lane packing so overlapping gene arrows stack
# in separate lanes. Adds a `lane` column (1 = bottom).
#
# Example:
#   assign_gene_lanes(data.frame(start=c(1,5,20), end=c(10,15,30)))
#   # -> lanes = c(1, 2, 1)
assign_gene_lanes <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  df <- df[order(df$start), , drop = FALSE]
  lanes <- integer(nrow(df))
  lane_ends <- numeric(0)
  for (i in seq_len(nrow(df))) {
    placed <- FALSE
    for (k in seq_along(lane_ends)) {
      if (df$start[i] > lane_ends[k]) {
        lanes[i] <- k; lane_ends[k] <- df$end[i]; placed <- TRUE; break
      }
    }
    if (!placed) { lanes[i] <- length(lane_ends) + 1L; lane_ends <- c(lane_ends, df$end[i]) }
  }
  df$lane <- lanes
  df
}
