# Module / credible-set enumeration and per-CS master-annotation lookups.

# module_cs_list(): return one row per credible set for a module. Prefers
# the CS column of the LocusZoom CSV (when `has_lz`); falls back to unique
# (cell, gene) pairs discovered in either the LZ file or the master
# annotation. Columns: cs, cell, eGene, eGene_symbol, chrom, lead_pos.
# `cell_filter` / `gene_filter` vectors (NULL = no filter) are applied on
# cell and on (eGene | eGene_symbol) respectively.
#
# Example:
#   module_cs_list("UMAP_lz.csv", meta, "M_18448", "cell", NULL, NULL)
#   # -> data.table with one row:
#   #    cs=..., cell="T_CD4_CM", eGene="ENSG00000101546",
#   #    eGene_symbol="RBFA", chrom="18", lead_pos=80005273
module_cs_list <- function(lz_path, meta, mod_id, anno_col,
                           cell_filter = NULL, gene_filter = NULL) {
  rows <- NULL
  lz_tbl <- if (!is.null(lz_path) && !is.na(lz_path) && lz_path != "NA")
              load_lz_file(lz_path) else NULL
  
  if (!is.null(lz_tbl) && "CS" %in% names(lz_tbl)) {
    uniq_cs <- unique(as.character(lz_tbl$CS))
    parsed  <- lapply(uniq_cs, parse_cs)
    rows <- data.table(
      cs       = sapply(parsed, `[[`, "cs"),
      cell     = sapply(parsed, `[[`, "cell"),
      eGene    = sapply(parsed, `[[`, "eGene"),
      chrom    = sapply(parsed, `[[`, "chrom"),
      lead_pos = sapply(parsed, `[[`, "lead_pos")
    )
  } else if (!is.null(lz_tbl) && all(c("CELL", "GENE") %in% names(lz_tbl))) {
    rows <- unique(lz_tbl[, .(cell = as.character(CELL),
                               eGene = sub("\\.\\d+$","", as.character(GENE)))])
    rows[, `:=`(cs = NA_character_, chrom = NA_character_,
                lead_pos = NA_real_)]
  } else if (!is.null(lz_tbl)) {
    # LZ file without any CS/CELL/GENE context: one implicit CS.
    rows <- data.table(cs = NA_character_, cell = NA_character_,
                       eGene = NA_character_, chrom = NA_character_,
                       lead_pos = NA_real_)
  } else {
    # No LZ file: fall back to the master annotation (requires `meta`).
    if (is.null(meta)) return(data.table())
    pairs <- module_pairs(meta, mod_id, anno_col, NULL, NULL)
    if (nrow(pairs) == 0) return(data.table())
    rows <- data.table(
      cs = NA_character_,
      cell = pairs$cell,
      eGene = pairs$eGene,
      chrom = NA_character_,
      lead_pos = NA_real_
    )
  }
  
  # Attach eGene_symbol from master annotation when possible.
  if (!is.null(meta) && all(c("module","cell","eGene","eGene_symbol") %in% names(meta))) {
    sym <- unique(meta[module == mod_id,
                       .(cell = as.character(cell),
                         eGene = sub("\\.\\d+$","", as.character(eGene)),
                         eGene_symbol = as.character(eGene_symbol))])
    rows <- merge(rows, sym, by = c("cell", "eGene"), all.x = TRUE)
  } else {
    rows[, eGene_symbol := NA_character_]
  }
  # If no gene symbol found, fall back to the ENSG id.
  rows[is.na(eGene_symbol) | eGene_symbol == "", eGene_symbol := eGene]
  
  if (!is.null(cell_filter) && length(cell_filter) > 0)
    rows <- rows[cell %in% cell_filter]
  if (!is.null(gene_filter) && length(gene_filter) > 0) {
    keep <- rows$eGene %in% gene_filter |
            (!is.na(rows$eGene_symbol) & rows$eGene_symbol %in% gene_filter)
    rows <- rows[keep]
  }
  rows
}

# module_pairs(): distinct (cell, eGene, eGene_symbol) triples for a module,
# optionally filtered by `cell_filter` (character vector of cells to keep)
# and `gene_filter` (character vector matched against eGene_symbol or bare
# eGene ENSG). Returns a data.table with columns cell, eGene, eGene_symbol;
# empty when nothing survives.
#
# Example:
#   module_pairs(meta, "M_18448", cell_filter=NULL, gene_filter=NULL)
#   # -> data.table(cell="T_CD4_CM", eGene="ENSG00000101546",
#   #               eGene_symbol="RBFA")
module_pairs <- function(df, mod_id, anno_col,
                         cell_filter = NULL, gene_filter = NULL) {
  sub_df <- df %>% filter(module == mod_id)
  if (nrow(sub_df) == 0) return(data.table(cell = character(),
                                           eGene = character(),
                                           eGene_symbol = character()))
  pairs <- sub_df %>%
    mutate(
      cell_val    = if (anno_col %in% names(.)) as.character(.data[[anno_col]]) else NA_character_,
      eGene_val   = if ("eGene" %in% names(.)) as.character(eGene) else NA_character_,
      eSymbol_val = if ("eGene_symbol" %in% names(.)) as.character(eGene_symbol) else NA_character_
    ) %>%
    filter(!is.na(cell_val) & !is.na(eGene_val)) %>%
    distinct(cell_val, eGene_val, .keep_all = FALSE) %>%
    transmute(cell = cell_val, eGene = eGene_val,
              eGene_symbol = NA_character_)
  # Re-attach a gene symbol for each (cell, eGene) from the first matching row.
  if ("eGene_symbol" %in% names(sub_df)) {
    sym_lookup <- sub_df %>%
      mutate(cell_val = if (anno_col %in% names(.)) as.character(.data[[anno_col]]) else NA_character_,
             eGene_val = if ("eGene" %in% names(.)) as.character(eGene) else NA_character_) %>%
      distinct(cell_val, eGene_val, .keep_all = TRUE) %>%
      select(cell_val, eGene_val, eGene_symbol = eGene_symbol)
    pairs <- pairs %>% left_join(sym_lookup,
                                 by = c("cell" = "cell_val", "eGene" = "eGene_val")) %>%
      mutate(eGene_symbol = coalesce(eGene_symbol.y, eGene_symbol.x)) %>%
      select(cell, eGene, eGene_symbol)
  }
  out <- as.data.table(pairs)
  if (!is.null(cell_filter) && length(cell_filter) > 0) {
    out <- out[cell %in% cell_filter]
  }
  if (!is.null(gene_filter) && length(gene_filter) > 0) {
    keep <- out$eGene %in% gene_filter |
            (!is.na(out$eGene_symbol) & out$eGene_symbol %in% gene_filter)
    out <- out[keep]
  }
  out
}

# get_module_locus_info(): locus coordinates + lead SNP for a (module, cell,
# gene) triple. Returns a list with chrom, gene_start, gene_end, gene_strand,
# gene_symbol, lead_pos (any field may be NA).
#
# Example:
#   get_module_locus_info(meta, "M_18448", "T_CD4_CM", "ENSG00000101546",
#                         anno_col="cell", gene_col="eGene_symbol")
#   # -> list(chrom="chr18", gene_start=79850000, gene_end=79960000,
#   #         gene_strand="+", gene_symbol="RBFA", lead_pos=80005273)
get_module_locus_info <- function(df, mod_id, focal_cell = NULL,
                                  focal_gene = NULL, anno_col = "cell",
                                  gene_col = "eGene_symbol") {
  sub_df <- df %>% filter(module == mod_id)
  if (!is.null(focal_cell) && anno_col %in% names(sub_df))
    sub_df <- sub_df %>% filter(.data[[anno_col]] == focal_cell)
  if (!is.null(focal_gene) && "eGene" %in% names(sub_df))
    sub_df <- sub_df %>% filter(eGene == focal_gene)
  if (nrow(sub_df) == 0) return(NULL)
  if ("cs_max_pip" %in% names(sub_df))
    sub_df <- sub_df %>% arrange(desc(as.numeric(cs_max_pip)))
  row1 <- sub_df[1, ]
  pick <- function(col) if (col %in% names(row1)) row1[[col]][1] else NA
  list(
    chrom       = pick("chrom"),
    gene_start  = suppressWarnings(as.numeric(pick("eGene_start"))),
    gene_end    = suppressWarnings(as.numeric(pick("eGene_end"))),
    gene_strand = pick("eGene_strand"),
    gene_symbol = if (gene_col %in% names(row1)) row1[[gene_col]][1] else pick("eGene_symbol"),
    lead_pos    = suppressWarnings(as.numeric(pick("most_likely_snp_pos")))
  )
}

# get_module_betas(): per-celltype beta table for a (module, gene). Filters
# master rows to module==mod_id AND eGene==focal_gene before collapsing to
# one row per cell, so the Beta UMAP shows the effect of this specific gene
# across all cell types.
#
# Example:
#   get_module_betas(meta, "M_18448", focal_gene="ENSG00000101546",
#                    anno_col="cell", join_col="celltype_2",
#                    gene_col="eGene_symbol")
#   # -> tibble with one row per celltype: beta_val, gene_sym_col, snp_id, ...
get_module_betas <- function(df, mod_id, focal_gene, anno_col, join_col, gene_col) {
  sub_df <- df %>% filter(module == mod_id)
  if (!is.null(focal_gene) && "eGene" %in% names(sub_df))
    sub_df <- sub_df %>% filter(eGene == focal_gene)
  if (nrow(sub_df) == 0) return(NULL)
  if ("cs_max_pip" %in% names(sub_df)) sub_df <- sub_df %>% arrange(desc(as.numeric(cs_max_pip)))
  sub_df %>%
    distinct(.data[[anno_col]], .keep_all = TRUE) %>%
    mutate(
      final_beta   = if ("cs_top_snp_beta" %in% names(.)) cs_top_snp_beta else if ("most_likely_snp_beta" %in% names(.)) most_likely_snp_beta else beta,
      gene_sym_col = if (gene_col %in% names(.)) .data[[gene_col]] else if ("eGene" %in% names(.)) eGene else "Unknown_Gene",
      mod_id_col   = mod_id,
      snp_id       = if ("most_likely_snp_rsID" %in% names(.)) most_likely_snp_rsID else NA,
      pval_exact   = if ("most_likely_snp_chisq" %in% names(.)) pchisq(as.numeric(most_likely_snp_chisq), df = 1, lower.tail = FALSE) else NA
    ) %>%
    select(all_of(anno_col), beta_val = final_beta, gene_sym_col, mod_id_col, snp_id, pval_exact) %>%
    rename(!!join_col := all_of(anno_col))
}
