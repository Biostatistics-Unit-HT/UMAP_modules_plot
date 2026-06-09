# UMAP coloured by per-cell gene expression (yellow-red ramp).

# get_palette_colors(): return a colour ramp for the requested palette name.
# Falls back to viridis if the palette is unknown.
#
# Example:
#   get_palette_colors("magma", 256)
#   # -> character vector of 256 hex colours from magma.
get_palette_colors <- function(name = "yellowred", n = 256L) {
  name <- tolower(as.character(name))
  if (name %in% c("yellowred", "yelloworangered", "ylorrd", "greyred")) {
    return(grDevices::colorRampPalette(c("#D9D9D9", "#F0E1A8", "#FED976",
                                         "#FEB24C", "#FD8D3C", "#F03B20",
                                         "#BD0026"))(n))
  }
  if (name %in% c("brightyellowred", "ylorrd_bright")) {
    return(grDevices::colorRampPalette(c("#FFF5B1", "#FED98E", "#FE9929",
                                         "#D95F0E", "#B30000"))(n))
  }
  if (!requireNamespace("viridisLite", quietly = TRUE)) {
    return(grDevices::colorRampPalette(c("#D9D9D9", "#FED976", "#FEB24C",
                                         "#FD8D3C", "#F03B20", "#BD0026"))(n))
  }
  opt_map <- c(viridis = "viridis", magma = "magma", inferno = "inferno",
               plasma = "plasma", turbo = "turbo", cividis = "cividis")
  opt <- if (name %in% names(opt_map)) opt_map[[name]] else "viridis"
  viridisLite::viridis(n, option = opt)
}

# clip_expression_values(): clip the upper tail of expression at a quantile so
# outliers do not squash the colour scale. Returns clipped values and limits.
#
# Example:
#   clip_expression_values(c(0, 0, 0, 5, 100), 0.95)
#   # -> list(values = ..., vmin = 0, vmax = quantile of non-zero)
clip_expression_values <- function(values, clip_q = 0.99) {
  vals <- as.numeric(values)
  finite <- vals[is.finite(vals)]
  if (length(finite) == 0L) return(list(values = vals, vmin = 0, vmax = 1))
  vmin <- min(finite, na.rm = TRUE)
  if (is.numeric(clip_q) && length(clip_q) == 1L &&
      is.finite(clip_q) && clip_q > 0 && clip_q < 1) {
    vmax <- as.numeric(stats::quantile(finite, probs = clip_q, na.rm = TRUE))
  } else {
    vmax <- max(finite, na.rm = TRUE)
  }
  if (!is.finite(vmax) || vmax <= vmin) vmax <- vmin + 1e-9
  vals[is.finite(vals) & vals > vmax] <- vmax
  list(values = vals, vmin = vmin, vmax = vmax)
}

# expr_meta_columns(): default metadata column names to exclude when resolving
# gene-expression columns in a wide UMAP table.
#
# Example:
#   expr_meta_columns("celltype_2", "celltype_3")
#   # -> c("cell_barcode", "cell_id", "UMAP_1", "UMAP_2", "celltype_2", "celltype_3")
expr_meta_columns <- function(join_col, label_col = NULL) {
  unique(c("cell_barcode", "cell_id", "UMAP_1", "UMAP_2", join_col, label_col))
}

# resolve_expr_gene_column(): pick the expression column from a wide UMAP table.
# Priority: match eGene (with/without Ensembl version), then gene_symbol, then
# the sole non-metadata column (single-gene files).
#
# Example:
#   resolve_expr_gene_column(c("cell_id","UMAP_1","UMAP_2","celltype_2",
#                              "ENSG00000101546"),
#                            "ENSG00000101546",
#                            expr_meta_columns("celltype_2"))
#   # -> "ENSG00000101546"
resolve_expr_gene_column <- function(col_names, eGene, meta_cols,
                                     gene_symbol = NULL) {
  gene_cols <- setdiff(col_names, meta_cols)
  if (length(gene_cols) == 0L) {
    stop("No gene-expression column found in --umap (only metadata columns present).",
         call. = FALSE)
  }
  bare_egene <- if (!is.null(eGene) && !is.na(eGene) && nzchar(as.character(eGene)))
    sub("\\.\\d+$", "", as.character(eGene)) else NA_character_

  match_gene <- function(cols, target_bare, target_full = NULL) {
    hits <- cols[
      (!is.na(target_full) && cols == target_full) |
      (!is.na(target_bare) && sub("\\.\\d+$", "", cols) == target_bare)
    ]
    unique(hits)
  }

  if (!is.na(bare_egene)) {
    hits <- match_gene(gene_cols, bare_egene, eGene)
    if (length(hits) == 1L) return(hits)
    if (length(hits) > 1L) {
      stop(sprintf("Multiple expression columns match eGene '%s': %s",
                   eGene, paste(hits, collapse = ", ")), call. = FALSE)
    }
  }
  if (!is.null(gene_symbol) && nzchar(as.character(gene_symbol)) &&
      as.character(gene_symbol) %in% gene_cols) {
    return(as.character(gene_symbol))
  }
  if (length(gene_cols) == 1L) {
    if (!is.na(bare_egene)) {
      col_bare <- sub("\\.\\d+$", "", gene_cols)
      if (col_bare != bare_egene) {
        warning(sprintf(
          "Using sole expression column '%s' (does not match eGene '%s').",
          gene_cols, eGene), call. = FALSE)
      }
    }
    return(gene_cols)
  }
  stop(sprintf(
    "Cannot resolve expression column for eGene '%s'. Candidates: %s",
    ifelse(is.na(bare_egene), "NA", bare_egene),
    paste(gene_cols, collapse = ", ")), call. = FALSE)
}

# resolve_gene_columns(): list all non-metadata gene columns (standalone batch).
#
# Example:
#   resolve_gene_columns(c("UMAP_0","celltype_2","ENSG1","ENSG2"),
#                        meta_cols = c("UMAP_0","UMAP_1","celltype_2"),
#                        requested = NULL)
#   # -> c("ENSG1", "ENSG2")
resolve_gene_columns <- function(all_cols, meta_cols, requested = NULL) {
  if (!is.null(requested) && nzchar(requested)) {
    want <- trimws(strsplit(requested, ",", fixed = TRUE)[[1]])
    want <- want[nzchar(want)]
    missing <- setdiff(want, all_cols)
    if (length(missing)) {
      stop(sprintf("These --genes are not columns of the input: %s",
                   paste(missing, collapse = ", ")), call. = FALSE)
    }
    return(want)
  }
  setdiff(all_cols, meta_cols)
}

# plot_expression(): UMAP coloured by per-cell expression from gene_col.
# Cells with expression == 0 or NA are drawn as grey background; positive
# cells use a yellow-red ramp. Panel title is NULL (orchestrator supplies
# the row-level figure title).
#
# Example:
#   plot_expression(df, "ENSG00000101546", 0.25, "celltype_2",
#                   show_legend = TRUE, use_raster = TRUE, show_labels = TRUE)
#   # -> ggplot with expression-coloured UMAP
plot_expression <- function(df, gene_col, pt_size, join_col,
                            show_legend = TRUE, use_raster = FALSE,
                            show_labels = FALSE, label_cells_only = NULL,
                            label_col = NULL, log1p = FALSE,
                            clip_q = 0.99, palette = "yellowred",
                            zero_color = "grey85") {
  if (!gene_col %in% names(df)) {
    stop(sprintf("Expression column '%s' not found in UMAP data.", gene_col),
         call. = FALSE)
  }
  plot_df <- data.frame(
    UMAP_1 = as.numeric(df$UMAP_1),
    UMAP_2 = as.numeric(df$UMAP_2),
    expr   = as.numeric(df[[gene_col]]),
    stringsAsFactors = FALSE
  )
  group_col <- if (!is.null(label_col) && nzchar(label_col) &&
                   label_col %in% names(df)) label_col else join_col
  if (group_col %in% names(df)) plot_df[[group_col]] <- df[[group_col]]
  if (join_col %in% names(df))   plot_df[[join_col]]   <- df[[join_col]]

  plot_df <- plot_df[is.finite(plot_df$UMAP_1) & is.finite(plot_df$UMAP_2), ,
                     drop = FALSE]
  if (isTRUE(log1p)) plot_df$expr <- log1p(pmax(plot_df$expr, 0, na.rm = FALSE))

  is_pos  <- is.finite(plot_df$expr) & plot_df$expr > 0
  zero_df <- plot_df[!is_pos, , drop = FALSE]
  pos_df  <- plot_df[ is_pos, , drop = FALSE]

  if (nrow(pos_df) > 0L) {
    clipped <- clip_expression_values(pos_df$expr, clip_q = clip_q)
    pos_df$expr <- clipped$values
    vmin <- clipped$vmin
    vmax <- clipped$vmax
    pos_df <- pos_df[order(pos_df$expr), , drop = FALSE]
  } else {
    vmin <- 0
    vmax <- 1
  }

  pal_colors <- get_palette_colors(palette, n = 256L)
  raster_ok  <- isTRUE(use_raster) && requireNamespace("scattermore", quietly = TRUE)
  pt_pixels  <- scattermore_pointsize(pt_size)
  pt_geom    <- max(0.15, 1.2 + 5 * as.numeric(pt_size))

  p <- ggplot()
  if (nrow(zero_df) > 0L) {
    if (raster_ok) {
      p <- p + scattermore::geom_scattermore(
        data = zero_df, mapping = aes(x = UMAP_1, y = UMAP_2),
        color = zero_color, pointsize = pt_pixels, pixels = c(2048, 2048))
    } else {
      p <- p + geom_point(data = zero_df, mapping = aes(x = UMAP_1, y = UMAP_2),
                          color = zero_color, size = pt_geom, stroke = 0,
                          alpha = 0.85)
    }
  }
  if (nrow(pos_df) > 0L) {
    if (raster_ok) {
      p <- p + scattermore::geom_scattermore(
        data = pos_df, mapping = aes(x = UMAP_1, y = UMAP_2, color = expr),
        pointsize = pt_pixels, pixels = c(2048, 2048))
    } else {
      p <- p + geom_point(data = pos_df, mapping = aes(x = UMAP_1, y = UMAP_2,
                                                       color = expr),
                          size = pt_geom, stroke = 0, alpha = 0.95)
    }
  }

  legend_title <- if (isTRUE(log1p)) "log1p(expr)" else "Expression"
  p <- p + scale_color_gradientn(colors = pal_colors,
                                 limits = c(vmin, vmax),
                                 oob    = scales::squish,
                                 na.value = zero_color,
                                 name = legend_title) +
    scale_x_continuous(breaks = NULL) +
    scale_y_continuous(breaks = NULL) +
    coord_equal(clip = "off") +
    labs(title = NULL, x = NULL, y = NULL) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          axis.text = element_blank(), axis.ticks = element_blank(),
          legend.title = element_text(size = 12, face = "plain"),
          legend.text  = element_text(size = 10))

  if (!show_legend) p <- p + theme(legend.position = "none")

  if (show_labels && group_col %in% names(plot_df)) {
    centroids <- plot_df %>%
      group_by(.data[[group_col]]) %>%
      summarise(UMAP_1 = median(UMAP_1, na.rm = TRUE),
                UMAP_2 = median(UMAP_2, na.rm = TRUE),
                .groups = "drop") %>%
      filter(!is.na(.data[[group_col]]) &
             nzchar(as.character(.data[[group_col]]))) %>%
      mutate(label_text = .data[[group_col]])
    if (!is.null(label_cells_only) && length(label_cells_only) > 0) {
      if (identical(group_col, join_col)) {
        keep_lab <- as.character(centroids[[group_col]]) %in% label_cells_only
      } else if (join_col %in% names(plot_df)) {
        keep_groups <- plot_df %>%
          dplyr::filter(as.character(.data[[join_col]]) %in% label_cells_only) %>%
          dplyr::pull(.data[[group_col]]) %>% as.character() %>% unique()
        keep_lab <- as.character(centroids[[group_col]]) %in% keep_groups
      } else {
        keep_lab <- rep(TRUE, nrow(centroids))
      }
      centroids <- centroids[keep_lab, , drop = FALSE]
    }
    if (nrow(centroids) > 0) {
      p <- p + geom_text_repel(
        data = centroids,
        aes(x = UMAP_1, y = UMAP_2, label = label_text),
        inherit.aes = FALSE, size = 7.5,
        family = "Helvetica", fontface = "plain",
        bg.color = "white", bg.r = 0.15,
        color = "black", min.segment.length = 0)
    }
  }
  p
}

# build_gene_umap(): standalone helper for plot_gene_umaps.R — one gene column
# with configurable coordinate and cell-type column names.
#
# Example:
#   build_gene_umap(df, "ENSG00000101546", "UMAP_0", "UMAP_1", "celltype_2",
#                   pt_size = 0.25, use_raster = TRUE, log1p = FALSE,
#                   clip_q = 0.99, palette = "yellowred")
#   # -> ggplot titled with the gene name
build_gene_umap <- function(df, gene, x_col, y_col, ct_col,
                            pt_size, use_raster, log1p, clip_q,
                            palette = "yellowred", zero_color = "grey85",
                            label_celltypes = FALSE,
                            font_family = "Helvetica",
                            base_text_size = 18) {
  if (!gene %in% names(df)) stop(sprintf("Column '%s' not found.", gene))
  tmp <- df
  tmp$UMAP_1 <- as.numeric(df[[x_col]])
  tmp$UMAP_2 <- as.numeric(df[[y_col]])
  tmp[[gene]] <- as.numeric(df[[gene]])
  if (ct_col %in% names(df) && ct_col != "celltype_2") {
    tmp$celltype_2 <- df[[ct_col]]
    label_join <- "celltype_2"
  } else if (ct_col %in% names(df)) {
    label_join <- ct_col
  } else {
    label_join <- "celltype_2"
  }
  p <- plot_expression(
    tmp, gene, pt_size, label_join,
    show_legend = TRUE, use_raster = use_raster,
    show_labels = isTRUE(label_celltypes),
    label_col = label_join, log1p = log1p, clip_q = clip_q,
    palette = palette, zero_color = zero_color)
  bs <- as.numeric(base_text_size)
  p + labs(title = gene) +
    theme(text = element_text(family = font_family, face = "plain", size = bs),
          plot.title = element_text(family = font_family, face = "plain",
                                    size = bs + 4, hjust = 0),
          plot.title.position = "plot",
          legend.title = element_text(family = font_family, face = "plain",
                                      size = bs),
          legend.text  = element_text(family = font_family, face = "plain",
                                      size = bs - 2),
          legend.key.height = grid::unit(1.2, "cm"),
          legend.key.width  = grid::unit(0.7, "cm"))
}
