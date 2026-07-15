# `plot_umap_simplified_multimodules.R`

R script that builds **UMAP figures** with beta-colored panels for one or two populations side-by-side, supporting multiple modules in a single run. Gene names, SNP IDs, and p-values are extracted automatically from master annotation files. Output is **PDF** (default) or **PNG**, with **rasterization** of points enabled by default for large cell counts.

```
[ LocusZoom ]   →   [ UMAP (beta or expression) ]   →   [ Coloc betas ]
```

- **LocusZoom** – `-log10(P)` vs genomic position, with
  - the **module-level most-likely SNP** (from the master table, or the top `-log10(P)` SNP when no module is provided) highlighted as a purple diamond; the same SNP anchors the LD r² colouring;
  - LD r² bins (`<0.2` navy, `0.2-0.4` light blue, `0.4-0.6` green, `0.6-0.8` orange, `≥0.8` red, `NA` gray) drawn as classic LocusZoom colours with bin-dependent point sizes so high-LD points stand out;
  - the bottom-most LZ row of a module emits a zoom‑out connector strip that fans toward the merged annotation box below.
- **UMAP middle panel** – `--umap_color_mode beta` (default): per-cell-type QTL beta from `--master`. `--umap_color_mode expression`: per-cell expression from a gene column in `--umap` (shared code with `plot_gene_umaps.R`).
- **Coloc betas** – QTL‑beta vs disease‑beta scatter for the colocalised credible sets.
- **Merged annotation box** (once per module, spans the full grid width)
  – regulatory tracks overlapping the module’s window, with a dotted
  orange line at the module-level most-likely SNP.

Every panel / file is optional. The **only hard requirement** is `--lz_files` OR the classic `--umap` + `--master` + `--modules` triple. Columns disappear from the layout entirely when their file(s) are omitted (no empty placeholder slots).

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

## Quick start
Get the pvalue for the locus if you don't have already together with the LD (not mandatory this last one)

```bash
Rscript plot_umap_simplified_multimodules.R [options]
```

This generate the files UMAP_lz_values.csv and UMAP_lz_values.raw

## Required arguments (Pop 1)

| Option | Description |
|--------|-------------|
| `--umap_1` | Path to a **tab-separated** UMAP file for Pop 1. Expects columns **`UMAP_1`** and **`UMAP_2`**. If your file uses **`UMAP_0`** and **`UMAP_1`** instead, the script renames them automatically. |
| `--master_1` | Path to the **Master Table file** for Pop 1 . |
| `--modules_1` | **Comma-separated** list of module IDs for Pop 1 (e.g. `M_13916,M_20000`). Use `NA` as a placeholder to skip a slot in dual mode. |

---

## Optional arguments (Pop 2 — dual/side-by-side mode)

| Option | Default | Description |
|--------|---------|-------------|
| `--umap_2` | *(none)* | Path to UMAP TSV for Pop 2. |
| `--master_2` | *(none)* | Path to **Master Table file** Pop 2. |
| `--modules_2` | *(none)* | Comma-separated module IDs for Pop 2. Must contain the **same number of items** as `--modules_1`. |
| `--name_2` | `Pop2` | Display name for Pop 2 shown in panel titles. |

---

## General settings

| Option | Default | Description |
|--------|---------|-------------|
| `--colors` | *(none)* | Path to palette TSV. **Required only when `--show_ref` is used.** |
| `--name_1` | `Pop1` | Display name for Pop 1 shown in panel titles. |
| `--join_col` | `cell` | Column name in the UMAP file used to identify cell types. |
| `--umap_color_mode` | `beta` | Middle UMAP colouring: `beta` or `expression`. |
| `--log1p` | off | log1p transform before expression colouring. |
| `--clip_quantile` | `0.99` | Upper quantile clip for non-zero expression cells. |
| `--expr_palette` | `yellowred` | Expression colour ramp. |
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

## Input files

Everything is optional — the minimal run is just `--lz_files` (and `--out`). When a file is omitted, the corresponding panel is dropped from the grid (no empty placeholder). The table below lists the three "classic" inputs that were previously required; they now behave as described in the "Minimal vs full runs" section further below.

### UMAP file (`--umap_1`, `--umap_2`)

- TSV with at least **`UMAP_1`**, **`UMAP_2`**, and the column named by **`--join_col`** (default `cell`).
- Alternate naming **`UMAP_0`**, **`UMAP_1`** is supported and mapped to `UMAP_1` / `UMAP_2`.

| column      | description                                       |
|-------------|---------------------------------------------------|
| `UMAP_1`    | first UMAP coordinate                             |
| `UMAP_2`    | second UMAP coordinate                            |
| `<celltype>`| a cell‑type label column (name via `--join_col`)  |

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

## Optional side panels

### 4. `--beta_files` – Coloc beta CSVs  *(enables the coloc panel)*

Fully optional. If you do not pass `--beta_files`, the Coloc betas column is not drawn at all (there is no empty slot in the layout). `--z_files` is accepted as a deprecated alias.

When supplied, pass a comma‑separated list, **same length and order as `--modules`**, one CSV per module. Use the literal string `NA` for a module that has no beta file.

Expected columns (auto‑detected, case‑insensitive; legacy `z_qtl` / `z_disease` / `z_icd10` still accepted):

| column         | description                     |
|----------------|---------------------------------|
| `snp`          | SNP identifier                  |
| `beta_qtl`     | QTL beta (effect size)          |
| `beta_disease` | disease beta (effect size)      |

Example:

```csv
snp,beta_disease,beta_qtl
chr18:79925373:C:T,0.12,-0.08
chr18:79925702:A:G,-0.11, 0.09
```

### 5. `--lz_files` – LocusZoom P‑value CSVs  *(enables the LocusZoom panel)*

Comma‑separated list, **same length and order as `--modules`** (one CSV per module). Use `NA` to skip a module.

Two layouts are auto‑detected:

1. **5 columns with header** — `CHR, CELL, GENE, POS, P`.
2. **6 columns, no header** — `CS, CHR, CELL, GENE, POS, P`. The `CS`
   string encodes the credible‑set lead SNP as `chr:pos:ref:alt` somewhere in it (e.g `chr18::GH_meta_cardinal::T_CD4_CM:ENSG00000101546::chr18:80005273:C:T::L1`).
   When present, the embedded position is used as the lead SNP on the LocusZoom (overriding the master table's `most_likely_snp_pos`).

| column | description                                                    |
|--------|----------------------------------------------------------------|
| `CS`   | (layout 2 only) credible‑set id carrying the lead SNP          |
| `CHR`  | chromosome (plain `18` or `chr18`)                             |
| `CELL` | cell type the association was computed in                      |
| `GENE` | Ensembl gene ID (e.g. `ENSG00000101546`, version suffix OK)    |
| `POS`  | genomic position (bp, 1‑based)                                 |
| `P`    | **raw P‑value**; the script computes `-log10(P)` on the fly    |

A single file may pool SNPs for **multiple `(CELL, GENE)` combinations** — the script automatically filters each file to the current pair it is plotting (matching `CELL` exactly and `GENE` with the ENSG version suffix stripped). If a file has only one cell/gene, that works too.

Example (layout 2, no header):

```csv
chr18::GH_meta_cardinal::T_CD4_CM:ENSG00000101546::chr18:80005273:C:T::L1,chr18,T_CD4_CM,ENSG00000101546,79869986,0.254
chr18::GH_meta_cardinal::T_CD4_CM:ENSG00000101546::chr18:80005273:C:T::L1,chr18,T_CD4_CM,ENSG00000101546,79870170,0.775
...
```

### 6. `--summary_table` – per‑module disease summary (CSV)

Optional. Used only to enrich the Beta panel title with the colocalised disease and its β / P. One row per module with any of:

| column (case‑insensitive)                        | description                         |
|--------------------------------------------------|-------------------------------------|
| `module`                                         | module ID                           |
| `coloc_trait`                                    | disease / trait label               |
| `most_likely_snp`                                | disease lead SNP                    |
| `most_likely_beta_disease` / `beta.*disease`     | disease β                           |
| `most_likely_se_disease` / `se.*disease`         | disease SE (used to compute P)      |

Modules whose β and SE are both exactly 0 in this table are skipped.

### 7. `--gtf` – gene track for the LocusZoom panel *(optional but recommended)*

Accepts **either**:

1. A raw **GENCODE GTF** (e.g. `gencode.v49.annotation.gtf`).
   The first call builds a small cache next to the GTF named
   `<gtf>.coding_genes.tsv` (~20k rows, <1 min, runs once).
   Subsequent calls load the cache instantly.
2. A pre‑built cache TSV with header
   `chrom\tstart\tend\tstrand\tgene_name\tgene_id`.
   Loaded as‑is; the awk pre‑filter is skipped.

When supplied, the LocusZoom gene track shows **every protein‑coding gene** overlapping the plotted window (stacked into lanes to avoid overlap), with the module's focal eGene coloured red. Without `--gtf`, only the focal
eGene from the master table is drawn.

### 8. `--annotations` – regulatory tracks for the zoom panel *(optional)*

Comma‑separated list of annotation files. Each becomes extra lanes in the zoom panel underneath the gene track.

Two formats are auto‑detected:

| Mode | Layout | Width (inches) |
|------|--------|----------------|
| Dual mode (Pop 1 + Pop 2) | 2-column grid, Pop 1 left / Pop 2 right per row | 14 |
| Single mode with `--show_ref` | Reference panel left, beta grid right | `(cols + 1) × 6` |
| Single mode without reference | Beta grid only | `cols × 6` |

Height is `max(4.5 × n_rows, 5)` inches, where `n_rows` is the number of module rows in the grid (i.e. `ceil(n_panels / cols)`).

- **PDF** (default): **`cairo_pdf`** with vector output and rasterized scatter points.
- **PNG** (with `--png`): 300 DPI, white background, via **`ragg::agg_png`**.

---

## Example runs

### Minimal – LocusZoom only, no UMAP / master / modules

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
  --name_1 "UKB" \
  --umap_2 pop2_umap.tsv \
  --master_2 pop2_master.tsv \
  --modules_2 M_13916,M_20000 \
  --name_2 "GH"
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

## Output layout

- One row of the grid per **(module, cell, gene)** triple (after any `--cell` / `--gene` filter).
- Number of columns per row is `1 + has(beta_files) + has(lz_files)` (1–3 columns). Columns you don't supply are not drawn at all (no empty slots left behind).
- Each LocusZoom cell is itself a stack of sub‑panels:
  1. `-log10(P)` scatter with the highlighted lead SNP,
  2. zoom annotation panel (only if `--annotations` is supplied). Discrete tracks draw solid blocks; continuous tracks (BED with a numeric 5th column) draw score‑scaled bars.

The beta panel uses a diverging **blue → white → yellow/orange/red** gradient symmetric around zero (`± max(|beta|)` in the plotted data). Missing betas after the join are coerced to **0** for plotting; `NA` in the scale is shown as **grey**.

---

## Tips

- If the run is slow because of very large UMAPs, lower `--max_cells` (default `250 000`) or keep rasterisation on (default).
- Passing an already‑built `*.coding_genes.tsv` to `--gtf` is auto‑detected by its header (`chrom start end strand gene_name gene_id`) — no awk pre‑filter runs.
- To change the order of annotation lanes, concatenate the files in the order you want in `--annotations`; the first‑seen feature types end up at the bottom of the zoom panel.
- Continuous annotation scores are min‑max normalised **within each feature type in the current window**, so profiles fill their lane. If you need absolute values across modules, pre‑normalise the scores in your BED file.

---

## Generating the inputs with `extract_z_lz.py` (optional)

The R script consumes ready-made `--lz_files`, `--ld_files`, and (optionally) `--beta_files` CSVs. The companion Python script `extract_z_lz.py` pulls those tables straight out of an anndata / TileDB QTL store so you don't have to wire them up by hand. It is entirely optional — you can feed the R script anything with the expected columns.

### What it produces

Given one QTL module, the script writes:

| file                   | purpose                                                          | where it plugs in     |
|------------------------|------------------------------------------------------------------|-----------------------|
| `<out>.csv`            | headerless `CS, CHR, CELL, GENE, POS, P` rows (the LZ format)    | `--lz_files`          |
| `<out>.raw`            | plink2 `--export A include-alt` genotype dosages for every SNP  | `--ld_files`          |
| `<out>_betas.csv`      | `snp, beta_disease, beta_qtl` coloc table (only if `--dis_adata`/`--dis_gz` set) | `--beta_files`        |

One call emits one set of files per module; the R script then accepts a comma-separated list of these (one per module) in the corresponding CLI flags.

### Inputs

| arg                  | what                                                               |
|----------------------|--------------------------------------------------------------------|
| `--qtl_module_adata` | H5AD with QTL modules (`.obs` has `phenotype_id` and `list_of_cs`) |
| `--qtl_cs_adata`     | H5AD with per-CS SNP tables (`.var` has `chr`, `pos`; `.obs` has `chr`, `start`, `end`) |
| `--tiledb`           | TileDB array with raw QTL P-values (queried by `(chr, cell, gene, start:end)`) |
| `--qtl_module`       | Module id to extract (e.g. `M_18448`, `M_3934`)                   |
| `--dis_adata`        | *(optional)* H5AD with disease fine-mapping, triggers the beta CSV |
| `--dis_cs`           | *(optional)* disease CS name to align against the QTL            |
| `--safeld`           | Path where a pgen file for ld are stored |
| `--out`              | output prefix (no extension)                                      |

### Example
python extract_z_lz.py --qtl_module_adata /ssu/bsssu/anndata_finemapping_repository/GH_meta_cardinal_F3_qtl_cis_07_01_hypr_modules_modules_JF_PIP_isCS99.h5ad --qtl_module M_3934 --out UMAP_lz_values --tiledb /project/cardinal/QTLs/TileDB_GH_f3_05_1_26/TileDB_tiledb_gh_meta_celltype2_f3_05_01_26 --qtl_cs_adata /ssu/bsssu/anndata_finemapping_repository/GH_meta_cardinal_F3_qtl_cis_07_01_anndata_PIP_isCS99.h5ad --safeld /project/cardinal/safeld_storage/GH_5k/final_output
```bash
python extract_z_lz.py \
  --qtl_module_adata /ssu/bsssu/anndata_finemapping_repository/GH_meta_cardinal_F3_qtl_cis_07_01_hypr_modules_modules_JF_PIP_isCS99.h5ad \
  --qtl_cs_adata     /ssu/bsssu/anndata_finemapping_repository/GH_meta_cardinal_F3_qtl_cis_07_01_anndata_PIP_isCS99.h5ad \
  --tiledb           /project/cardinal/QTLs/TileDB_GH_f3_05_1_26/TileDB_tiledb_gh_meta_celltype2_f3_05_01_26 \
  --qtl_module       M_18448 \
  --dis_adata        /ssu/.../GH:gwas:ICD10codes_anndata.h5ad \
  --dis_cs           "chr18::GH:gwas:ICD10codes::J15::chr18:80004668:A:T::L1" \
  --out              UMAP_lz_values
```

This writes `UMAP_lz_values.csv`, `UMAP_lz_values.raw`, and
`UMAP_lz_values_betas.csv` that you can then feed straight into:

```bash
Rscript plot_umap_simplified_multimodules_one_side.R \
  --lz_files    UMAP_lz_values.csv \
  --ld_files    UMAP_lz_values.raw \
  --beta_files  UMAP_lz_values_betas.csv \
  --out         umap_M18448.pdf --raster
```

### Python dependencies
You can use on the hpc the tdbsumstat conda environment or you need the following packages
```bash
pip install scanpy pandas tiledb scipy numpy
```

- **"Missing required arguments"** — Must provide all of `--umap_1`, `--master_1`, and `--modules_1`.
- **"--colors required with --show_ref"** — Provide `--colors` when using `--show_ref`.
- **"Pop 1 and Pop 2 module lists must have the same number of items"** — Ensure the comma-separated lists in `--modules_1` and `--modules_2` have equal length.
- **Module not found warnings** — Check that the module ID exists in the `module` column of the master annotation file.
- **Join failures / empty plots** — Ensure `--join_col` (and `--anno_join_col` if used) values match across the UMAP and master annotation files (same spelling and level set).
