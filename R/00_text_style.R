# Global text styling shared across every plot in this script. Applied to the
# final patchwork via `& big_helvetica_theme()` so it cascades to every
# sub-panel without rewriting each individual theme() call inside the
# per-plot helpers.

# big_helvetica_theme(): theme() partial that forces every text element to
# Helvetica, plain face (no bold) and a much larger size. Override individual
# pieces with the `base_size` argument.
#
# Example:
#   p <- p_lz / p_beta
#   p & big_helvetica_theme()           # default sizes
#   p & big_helvetica_theme(base_size = 22)  # even larger
big_helvetica_theme <- function(base_size = 18, family = "Helvetica") {
  fam <- family
  bs  <- as.numeric(base_size)
  ggplot2::theme(
    text          = ggplot2::element_text(family = fam, face = "plain", size = bs),
    # Titles are LEFT-aligned (hjust = 0) so a wide title can't overflow into
    # the adjacent column when the LocusZoom column is narrower than the
    # others. The original ggplot default `hjust = 0.5` re-centres the title
    # over its panel and a long title spills sideways into the next panel.
    plot.title    = ggplot2::element_text(family = fam, face = "plain",
                                          size = bs + 2, hjust = 0,
                                          margin = ggplot2::margin(b = bs * 0.25)),
    plot.subtitle = ggplot2::element_text(family = fam, face = "plain",
                                          size = bs - 2, hjust = 0),
    # Anchor titles to the whole plot/slot (not just the inner panel) so a
    # multi-line title in a narrower column doesn't lean into the neighbour.
    plot.title.position    = "plot",
    plot.caption.position  = "plot",
    plot.caption  = ggplot2::element_text(family = fam, face = "plain", size = bs - 4),
    # Axis labels (titles) and tick text are bumped one step above the base
    # size so the chromosome / Mb / -log10(P) labels read clearly even when
    # the LocusZoom column shares width with the beta / Z panels.
    axis.title    = ggplot2::element_text(family = fam, face = "plain", size = bs + 3),
    axis.title.x  = ggplot2::element_text(family = fam, face = "plain", size = bs + 3),
    axis.title.y  = ggplot2::element_text(family = fam, face = "plain", size = bs + 3),
    axis.text     = ggplot2::element_text(family = fam, face = "plain", size = bs + 2),
    axis.text.x   = ggplot2::element_text(family = fam, face = "plain", size = bs + 2),
    axis.text.y   = ggplot2::element_text(family = fam, face = "plain", size = bs + 2),
    legend.title  = ggplot2::element_text(family = fam, face = "plain", size = bs - 2),
    legend.text   = ggplot2::element_text(family = fam, face = "plain", size = bs - 4),
    strip.text    = ggplot2::element_text(family = fam, face = "plain", size = bs)
  )
}
