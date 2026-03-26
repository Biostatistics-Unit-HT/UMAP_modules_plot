# `plot_umap_simplified.R`

R script that builds **UMAP figures** combining a **reference** plot (cells colored by a fixed palette per cell type) with one or two **beta-colored** panels for chosen module–eGene pairs. Output is **PNG** (default) or **PDF**, with optional **rasterization** of points for large cell counts.

## Requirements

Install these R packages before running:

| Package      | Role                                      |
|-------------|-------------------------------------------|
| `ggplot2`   | Plotting                                  |
| `dplyr`     | Data manipulation                         |
| `data.table`| Fast TSV reading (`fread`)                |
| `optparse`  | Command-line arguments                    |
| `patchwork` | Combining panels            |
| `ggrepel`   | Non-overlapping labels on reference/beta  |
| `ragg`      | PNG device (`agg_png`)                    |
| `scattermore` | Optional fast rasterized scatter layers |

Run the script with:

```bash
Rscript plot_umap_simplified.R [options]
```

---

## Required arguments

| Option | Description |
|--------|-------------|
| `--umap` | Path to a **tab-separated** file with UMAP coordinates. Expects columns **`UMAP_1`** and **`UMAP_2`**. If your file uses **`UMAP_0`** and **`UMAP_1`** instead, the script renames them to `UMAP_1` / `UMAP_2`. |
| `--colors` | Path to a **TSV** with the join column (see `--join_col`).
| `--annotations` | Path to a **TSV** with module/eGene annotations and beta values (see [Annotation file](#annotation-file)). |
| `--module_a` | Primary module ID (e.g. `M_13916`). |
| `--ensg_a` | Primary eGene Ensembl ID (e.g. `ENSG0000...`). |

---

## Optional arguments

| Option | Default | Description |
|--------|---------|-------------|
| `--module_b` | *(none)* | Secondary module ID; use together with `--ensg_b` for a second beta panel. |
| `--ensg_b` | *(none)* | Secondary eGene ID. |
| `--join_col` | `celltype_2` | Column present in **all three** inputs used to align cells / cell types. |
| `--out` | `umap_module_betas.png` | Output filename. For PDF, `.png` is replaced or `.pdf` is appended. |
| `--pt_size` | `0.25` | Point size (scatter layers use a slightly larger size when raster is on). |
| `--max_cells` | `250000` | If the merged table has more rows, a **random subsample** (seed `42`) is taken for speed. |
| `--pdf` | off | Save as **PDF** (`cairo_pdf`) instead of PNG. |
| `--raster` | off | Draw points with **`scattermore`** (recommended with `--pdf` for smaller files and faster drawing on huge point clouds). |
| `--no_ref` | off | Omit the **left** reference UMAP; layout becomes only the beta panel(s). |
| `--label_betas` | off | Add **cell-type labels** (medians per `join_col`) on the beta panel(s), similar to the reference plot. |

---

## Input file expectations

### UMAP file (`--umap`)

- TSV with at least **`UMAP_1`**, **`UMAP_2`**, and the column named by **`--join_col`** (default `celltype_2`).
- Alternate naming **`UMAP_0`**, **`UMAP_1`** is supported and mapped to `UMAP_1`, `UMAP_2`.

### Colors file (`--colors`)

- Must include **`--join_col`** and **`color_ct2`** (hex color per group). These drive the **“UMAP Reference”** panel.

### Annotation file (`--annotations`)

Used to attach **beta values** per `join_col` for a given `module` + `eGene`:

- Required conceptually: columns **`module`**, **`eGene`**, **`--join_col`**, and a beta source.
- Beta column: the script prefers **`most_likely_snp_beta`** if present; otherwise **`beta`**.
- If **`cs_max_pip`** exists, rows are ordered by it (descending) before deduplication so one row per `join_col` is kept.
- Missing combinations cause **`stop`** with a clear error.

---

## Layout and output sizes

| Mode | Layout | Width × height (inches) |
|------|--------|-------------------------|
| Reference + one beta | `reference \| beta A` | 11 × 5 |
| Reference + two betas | `reference \| (beta A / beta B)` | 11 × 8 |
| `--no_ref`, one beta | single panel | 5 × 5 |
| `--no_ref`, two betas | `beta A / beta B` | 5 × 8 |

- **PNG**: 300 DPI, white background, via **`ragg::agg_png`**.
- **PDF**: **`cairo_pdf`**; the script message mentions rasterized points when saving PDF—use **`--raster`** for actual rasterized scatter layers.

---

## Example commands

**Minimal (reference + one beta, PNG):**

```bash
Rscript plot_umap_simplified.R \
  --umap path/to/umap.tsv \
  --colors path/to/palette.tsv \
  --annotations path/to/annotations.tsv \
  --module_a M_13916 \
  --ensg_a ENSG00000123456 \
  --out my_umap.png
```

**Two modules + PDF + raster (smaller PDF):**

```bash
Rscript plot_umap_simplified.R \
  --umap umap.tsv \
  --colors colors.tsv \
  --annotations annotations.tsv \
  --module_a M_13916 \
  --ensg_a ENSG00000123456 \
  --module_b M_20000 \
  --ensg_b ENSG00000987654 \
  --pdf --raster \
  --out figure.png
```

**Beta panels only, with labels:**

```bash
Rscript plot_umap_simplified.R \
  --umap umap.tsv \
  --colors colors.tsv \
  --annotations annotations.tsv \
  --module_a M_13916 \
  --ensg_a ENSG00000123456 \
  --no_ref --label_betas \
  --out betas_only.png
```

**Custom join column:**

```bash
Rscript plot_umap_simplified.R \
  --umap umap.tsv \
  --colors colors.tsv \
  --annotations annotations.tsv \
  --join_col my_celltype_column \
  --module_a M_13916 \
  --ensg_a ENSG00000123456
```

---

## Beta color scale

The beta panel uses a diverging **blue → white → yellow/orange/red** gradient symmetric around zero (`± max(|beta|)` in the plotted data). Missing betas after join are coerced to **0** for plotting; `NA` in the scale is shown as **grey**.

---

## Troubleshooting

- **“Missing required arguments”** — Pass all of `--umap`, `--colors`, `--annotations`, `--module_a`, `--ensg_a`.
- **“No rows found in annotations for module … and eGene …”** — Check `module` / `eGene` strings and that the annotation table covers those pairs.
- **Join failures / empty plots** — Ensure **`--join_col`** values match across UMAP, colors, and annotations (same spelling and level set).
