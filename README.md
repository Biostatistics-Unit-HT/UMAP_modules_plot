# `plot_umap_simplified_multimodules.R`

R script that builds **UMAP figures** with beta-colored panels for one or two populations side-by-side, supporting multiple modules in a single run. Gene names, SNP IDs, and p-values are extracted automatically from master annotation files. Output is **PDF** (default) or **PNG**, with **rasterization** of points enabled by default for large cell counts.

## Requirements

Install these R packages before running:

| Package      | Role                                      |
|-------------|-------------------------------------------|
| `ggplot2`   | Plotting                                  |
| `dplyr`     | Data manipulation                         |
| `data.table`| Fast TSV reading (`fread`)                |
| `optparse`  | Command-line arguments                    |
| `patchwork` | Combining panels                          |
| `ggrepel`   | Non-overlapping labels on reference/beta  |
| `ragg`      | PNG device (`agg_png`)                    |
| `scattermore` | Fast rasterized scatter layers          |

Run the script with:

```bash
Rscript plot_umap_simplified_multimodules.R [options]
```

---

## Required arguments (Pop 1)

| Option | Description |
|--------|-------------|
| `--umap_1` | Path to a **tab-separated** UMAP file for Pop 1. Expects columns **`UMAP_1`** and **`UMAP_2`**. If your file uses **`UMAP_0`** and **`UMAP_1`** instead, the script renames them automatically. |
| `--master_1` | Path to the **Master Annotations TSV** for Pop 1 (see [Master Annotation file](#master-annotation-file)). |
| `--modules_1` | **Comma-separated** list of module IDs for Pop 1 (e.g. `M_13916,M_20000`). Use `NA` as a placeholder to skip a slot in dual mode. |

---

## Optional arguments (Pop 2 — dual/side-by-side mode)

| Option | Default | Description |
|--------|---------|-------------|
| `--umap_2` | *(none)* | Path to UMAP TSV for Pop 2. |
| `--master_2` | *(none)* | Path to Master Annotations TSV for Pop 2. |
| `--modules_2` | *(none)* | Comma-separated module IDs for Pop 2. Must contain the **same number of items** as `--modules_1`. |
| `--name_2` | `Pop2` | Display name for Pop 2 shown in panel titles. |

---

## General settings

| Option | Default | Description |
|--------|---------|-------------|
| `--colors` | *(none)* | Path to palette TSV. **Required only when `--show_ref` is used.** |
| `--name_1` | `Pop1` | Display name for Pop 1 shown in panel titles. |
| `--join_col` | `cell` | Column name in the UMAP file used to identify cell types. |
| `--anno_join_col` | *(same as `--join_col`)* | Column name in the Master Annotations file if it differs from `--join_col`. |
| `--gene_col` | `eGene_symbol` | Column name for the gene symbol in the Master Annotations file. |
| `--out` | `umap_plot` | Output filename **prefix** (extension is added automatically). |
| `--pt_size` | `0.25` | Point size for scatter layers. |
| `--max_cells` | `250000` | If the data has more rows, a random subsample (seed `42`) is taken. |

---

## Toggle flags

Note the defaults carefully — PDF output and rasterization are **on** by default, while the reference panel and label hiding are **off** by default:

| Flag | Default | Effect when passed |
|------|---------|--------------------|
| `--png` | off (PDF is default) | Save output as **PNG** (300 DPI, white background) instead of PDF. |
| `--no_raster` | off (raster is default) | **Disable** rasterized points; draw with standard `geom_point` instead. |
| `--show_ref` | off (reference hidden) | Add a **reference UMAP** panel colored by cell type (requires `--colors`). |
| `--no_labels` | off (labels shown) | **Hide** cell-type labels on beta panels. |

---

## Input file expectations

### UMAP file (`--umap_1`, `--umap_2`)

- TSV with at least **`UMAP_1`**, **`UMAP_2`**, and the column named by **`--join_col`** (default `cell`).
- Alternate naming **`UMAP_0`**, **`UMAP_1`** is supported and mapped to `UMAP_1` / `UMAP_2`.

### Colors file (`--colors`)

- Only required when **`--show_ref`** is used.
- Must include the **`--join_col`** column and **`color_ct2`** (hex color per group).

### Master Annotation file (`--master_1`, `--master_2`)

Used to extract **beta values**, gene names, SNP IDs, and p-values per module:

- Must include a **`module`** column and the annotation join column (see `--anno_join_col`).
- **Beta source priority**: `cs_top_snp_beta` → `most_likely_snp_beta` → `beta`.
- If **`cs_max_pip`** exists, rows are sorted descending by it before one row per cell type is kept.
- Optional columns: `cs_max_pip`, `most_likely_snp_rsID`, `most_likely_snp_chisq` (used to compute p-value), and the gene symbol column (see `--gene_col`; falls back to `eGene`).

---

## Beta extraction logic

For each module ID the `get_module_betas` function:

1. Filters the master annotation to rows where `module == mod_id`.
2. If `cs_max_pip` is present, sorts rows descending by it.
3. Deduplicates to **one row per `anno_join_col`** value.
4. Selects beta using the priority chain above.
5. Extracts gene symbol, SNP rsID, and p-value (computed from `most_likely_snp_chisq` via chi-squared with df = 1).

---

## Panel title format

Each beta panel title is composed as:

```
{PopName} | {GeneName}
[{ModuleID}]
{SNP_rsID} | P = {p-value}
```

---

## Layout and output sizes

| Mode | Layout | Width (inches) |
|------|--------|----------------|
| Dual mode (Pop 1 + Pop 2) | 2-column grid, Pop 1 left / Pop 2 right per row | 14 |
| Single mode with `--show_ref` | Reference panel left, beta grid right | `(cols + 1) × 6` |
| Single mode without reference | Beta grid only | `cols × 6` |

Height is `max(4.5 × n_rows, 5)` inches, where `n_rows` is the number of module rows in the grid (i.e. `ceil(n_panels / cols)`).

- **PDF** (default): **`cairo_pdf`** with vector output and rasterized scatter points.
- **PNG** (with `--png`): 300 DPI, white background, via **`ragg::agg_png`**.

---

## Example commands

**Minimal single-population (PDF, rasterized, no reference):**

```bash
Rscript plot_umap_simplified_multimodules.R \
  --umap_1 path/to/umap.tsv \
  --master_1 path/to/master_annotations.tsv \
  --modules_1 M_13916,M_20000 \
  --name_1 "European"
```

**Dual-population side-by-side:**

```bash
Rscript plot_umap_simplified_multimodules.R \
  --umap_1 pop1_umap.tsv \
  --master_1 pop1_master.tsv \
  --modules_1 M_13916,M_20000 \
  --name_1 "European" \
  --umap_2 pop2_umap.tsv \
  --master_2 pop2_master.tsv \
  --modules_2 M_13916,M_20000 \
  --name_2 "African"
```

**With reference panel and PNG output:**

```bash
Rscript plot_umap_simplified_multimodules.R \
  --umap_1 umap.tsv \
  --master_1 master.tsv \
  --modules_1 M_13916 \
  --colors palette.tsv \
  --show_ref --png \
  --out my_figure
```

**Skip a module in one population (use `NA`):**

```bash
Rscript plot_umap_simplified_multimodules.R \
  --umap_1 pop1_umap.tsv \
  --master_1 pop1_master.tsv \
  --modules_1 M_13916,NA,M_30000 \
  --name_1 "Pop1" \
  --umap_2 pop2_umap.tsv \
  --master_2 pop2_master.tsv \
  --modules_2 NA,M_20000,M_30000 \
  --name_2 "Pop2"
```

---

## Beta color scale

The beta panel uses a diverging **blue → white → yellow/orange/red** gradient symmetric around zero (`± max(|beta|)` in the plotted data). Missing betas after the join are coerced to **0** for plotting; `NA` in the scale is shown as **grey**.

---

## Troubleshooting

- **"Missing required arguments"** — Must provide all of `--umap_1`, `--master_1`, and `--modules_1`.
- **"--colors required with --show_ref"** — Provide `--colors` when using `--show_ref`.
- **"Pop 1 and Pop 2 module lists must have the same number of items"** — Ensure the comma-separated lists in `--modules_1` and `--modules_2` have equal length.
- **Module not found warnings** — Check that the module ID exists in the `module` column of the master annotation file.
- **Join failures / empty plots** — Ensure `--join_col` (and `--anno_join_col` if used) values match across the UMAP and master annotation files (same spelling and level set).
