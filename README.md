# UMAP + LocusZoom + Coloc plots

`plot_umap_simplified_multimodules_one_side.R` produces a composite figure for
one or more fine‑mapped credible sets in a single population.
Each CS of a module produces **one row** of the final grid with up to three panels:

```
[ LocusZoom ]   →   [ Beta UMAP ]   →   [ Coloc Z‑scores ]
```

- **LocusZoom** – `-log10(P)` vs genomic position, with
  - the **module-level most-likely SNP** (from the master table, or the top `-log10(P)` SNP when no module is provided) highlighted as a purple diamond; the same SNP anchors the LD r² colouring;
  - LD r² bins (`<0.2` navy, `0.2-0.4` light blue, `0.4-0.6` green, `0.6-0.8` orange, `≥0.8` red, `NA` gray) drawn as classic LocusZoom colours with bin-dependent point sizes so high-LD points stand out;
  - a gene-body track of every *protein-coding* gene in the window (focal eGene in red, others in navy; optionally restricted by `--gene_list`);
  - the bottom-most LZ row of a module emits a zoom‑out connector strip that fans toward the merged annotation box below.
- **Beta UMAP** – the population's UMAP coloured by the per-cell-type effect size for the current (cell, gene). Only drawn when **both** `--umap` and `--master` are provided.
- **Coloc Z-scores** – QTL‑Z vs disease‑Z scatter for the colocalised credible sets.
- **Merged annotation box** (once per module, spans the full grid width)
  – regulatory tracks overlapping the module’s window, with a dotted
  orange line at the module-level most-likely SNP.

Every panel / file is optional. The **only hard requirement** is `--lz_files` OR the classic `--umap` + `--master` + `--modules` triple. Columns disappear from the layout entirely when their file(s) are omitted (no empty placeholder slots).

---

## Quick start

```bash
Rscript plot_umap_simplified_multimodules_one_side.R \
  --umap          UMAP_GH_QCed_f3_ct2.tsv \
  --master        freeze3_gh_cs_full_annotations.20260306.tsv \
  --modules       M_18448 \
  --name          GH \
  --anno_join_col cell \
  --z_files       icd10_j15_zscore_merge_notnull.csv \
  --lz_files      j15_locuszoom_ENSG00000101546.csv \
  --gtf           gencode.v49.annotation.gtf \
  --gene_list     gene_list_j15.csv \
  --annotations   homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329.promoter.gff \
  --zoom_window   5000 \
  --out           umap_GH_J15.pdf \
  --raster
```

The output is a multi‑page vector PDF by default (`--png` switches to PNG). Remove `--z_files` and the right‑most Z-plot for coloc disappears from the layout. Add `--cell "T_CD4_CM"` and/or `--gene "ENSG00.."` to restrict which `(cell, gene)` rows of a module are plotted.

---

## R dependencies

```r
install.packages(c(
  "ggplot2", "dplyr", "data.table", "optparse",
  "patchwork", "ggrepel", "ragg", "scattermore"
))
```

Tested with R ≥ 4.3.

---

## Input files

Everything is optional — the minimal run is just `--lz_files` (and `--out`). When a file is omitted, the corresponding panel is dropped from the grid (no empty placeholder). The table below lists the three "classic" inputs that were previously required; they now behave as described in the "Minimal vs full runs" section further below.

### 1. `--umap` – UMAP coordinates table (TSV, *optional*)

One row per cell, tab separated. Expected columns:

| column      | description                                       |
|-------------|---------------------------------------------------|
| `UMAP_1`    | first UMAP coordinate                             |
| `UMAP_2`    | second UMAP coordinate                            |
| `<celltype>`| a cell‑type label column (name via `--join_col`)  |

If the file ships with `UMAP_0` / `UMAP_1` instead, they are renamed internally to `UMAP_1` / `UMAP_2`.

The default cell‑type column is `cell`; change it with
`--join_col <column_name>`.

### 2. `--master` – master credible-set annotations (TSV, *optional*)

Usually something like `freeze3_gh_cs_full_annotations.*.tsv`. One row per (module, cell type) with fine‑mapping info. When provided, this file drives the Beta UMAP (which needs `--umap` too), the gene-track eGene symbol, and the module-level most-likely SNP. When omitted, the script falls back to the top `-log10(P)` SNP in each LZ file as the module anchor. Columns consumed when present:

| column                        | used for                                           |
|-------------------------------|----------------------------------------------------|
| `module`                      | key to match `--modules`                           |
| `<join_col>` or `--anno_join_col` | cell‑type label to join against the UMAP       |
| `cs_top_snp_beta` / `most_likely_snp_beta` / `beta` | per cell‑type β shown on the Beta UMAP |
| `eGene_symbol` (or `--gene_col`) / `eGene` | gene name in titles                       |
| `eGene_start`, `eGene_end`, `eGene_strand`, `chrom` | focal eGene coords (LZ gene track) |
| `most_likely_snp_rsID`, `most_likely_snp_pos`      | lead SNP used in LocusZoom         |
| `most_likely_snp_chisq`       | exact P‑value in the Beta panel title              |
| `cs_max_pip`                  | used to pick the "best" row when multiple exist    |

### 3. `--modules` – module IDs (*optional*)

Comma-separated list: `--modules "M_18448,M_11607,M_10171"`.

Each module can have multiple `(cell, gene)` pairs in the master table (e.g. the same credible set mapping to several eGenes or cell types). The script enumerates those pairs and produces **one grid row per CS**. Use `--cell` and/or `--gene` to restrict them (see below). Without any filter every CS the module has is plotted.

**When `--modules` is omitted**, the script treats each `--lz_files`
entry as one implicit module (ID auto-derived from the filename, e.g.
`auto_UMAP_lz_values`) and uses the **top `-log10(P)` SNP in the file**
as the module anchor (diamond + LD source).

---

## Optional side panels

### 4. `--z_files` – Coloc Z‑score CSVs  *(enables the coloc panel)*

Fully optional. If you do not pass `--z_files`, the Coloc Z‑scores column is not drawn at all (there is no empty slot in the layout).

When supplied, pass a comma‑separated list, **same length and order as `--modules`**, one CSV per module. Use the literal string `NA` for a module that has no Z‑score file.

Expected columns (auto‑detected, case‑insensitive):

| column   | description                     |
|----------|---------------------------------|
| `snp`    | SNP identifier                  |
| `z_qtl`  | QTL Z‑score                     |
| `z_icd10` *or* `z_dis` | disease Z‑score    |

Example:

```csv
snp,z_icd10,z_qtl
chr18:79925373:C:T,4.26,-3.22
chr18:79925702:A:G,-4.19, 3.30
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

- **BED‑style** (no header) — 4 required columns, with an optional numeric 5th `score` column:

  ```
  chrom  start     end       feature_type              [score]
  18     80009802  80011199  promoter_flanking_region  0.73
  ```

  - 4 columns → a **discrete** track (each feature drawn as a solid block).
  - 5 columns with any numeric `score` → a **continuous** track: features
    are drawn as bars whose height is proportional to `score`, normalised
    within the track and window, so the lane looks like a distribution /
    profile around the lead SNP.

- **Standard 9‑column GFF/GTF** — `#` comment lines are skipped; column 1 is chromosome, column 3 is feature type, columns 4/5 are start/end, column 6 is interpreted as an optional score.

Every distinct `feature_type` across all supplied files becomes its own lane (row), with a consistent lane order across modules/rows.

`--zoom_window N` (default `5000`) controls the half‑width in bp of the zoom panel around the lead SNP.

### 9. `--gene_list` – gene‑track whitelist (CSV) *(optional)*

CSV with two columns, `CELL,GENE` (header required). `GENE` is a bare Ensembl ID (without version, e.g. `ENSG00000101546`). The script strips any `.N` suffix from the GTF's `gene_id` before matching, so GENCODE cache entries like `ENSG00000186092.7` still work.

```csv
CELL,GENE
T_CD4_CM,ENSG00000060069
T_CD4_CM,ENSG00000101544
T_CD4_CM,ENSG00000101546
...
```

When plotting a row for cell `X`, only genes whose `CELL == X` row exists in this file (by `GENE`) are drawn as arrows in the LocusZoom gene track. If the file has no rows for the current cell, the whitelist is disabled for that row (every gene in the window is drawn) and a notice is logged.

### 10. `--cell` / `--gene` – per‑module pair filters *(optional)*

Comma‑separated lists applied to the distinct `(cell, gene)` pairs discovered in the master table for each module.

- `--cell "T_CD4_CM,T_CD8_CM"` — keep only pairs whose cell is in the list.
- `--gene "RBFA,ATP9B,ENSG00000101546"` — keep pairs whose `eGene_symbol` **or** bare `eGene` (ENSG) is in the list, so you can mix symbols and IDs in a single argument.

Both filters are optional; omit them to plot every `(cell, gene)` pair the module has.

### 11. `--colors` – cell‑type colour palette (TSV) *(only with `--show_ref`)*

Required only when you pass `--show_ref`. Tab separated with at least:

| column          | description                                 |
|-----------------|---------------------------------------------|
| `<join_col>`    | cell‑type label (same column as the UMAP)   |
| `color_ct2`     | hex colour used in the reference UMAP       |

---

## All CLI options

| Option              | Required | Default         | What it does                                                              |
|---------------------|----------|-----------------|---------------------------------------------------------------------------|
| `--umap`            | ✔        | —               | UMAP coordinates TSV                                                      |
| `--master`          | ✔        | —               | Master table from joint finemappping TSV                                       |
| `--modules`         | ✔        | —               | Comma‑separated module IDs                                                |
| `--name`            |          | `Pop`           | Display name used in panel titles                                         |
| `--join_col`        |          | `celltype_2`    | Cell‑type column name in the UMAP                                         |
| `--anno_join_col`   |          | `--join_col`    | Cell‑type column name in the master table (if different)           |
| `--gene_col`        |          | `eGene_symbol`  | Gene‑symbol column in the master table                              |
| `--z_files`         |          | —               | Coloc Z‑score CSVs (comma‑separated, one per module, `NA` to skip). Omit to drop the panel entirely. |
| `--lz_files`        |          | —               | LocusZoom P‑value CSVs (comma‑separated, one per module, `NA` to skip). Omit to drop the panel entirely. |
| `--summary_table`   |          | —               | Per‑module disease summary CSV                                            |
| `--gtf`             |          | —               | Raw GTF or pre‑built `.coding_genes.tsv`                                  |
| `--gene_list`       |          | —               | CSV `CELL,GENE` whitelist filtering the LocusZoom gene track per cell     |
| `--cell`            |          | —               | Comma‑separated list of cell types to keep                                |
| `--gene`            |          | —               | Comma‑separated list of gene symbols or bare ENSGs to keep                |
| `--annotations`     |          | —               | Comma‑separated annotation files (BED 4/5‑col or 9‑col GFF)               |
| `--zoom_window`     |          | `5000`          | Half‑width in bp for the zoom annotation panel                            |
| `--colors`          | (with `--show_ref`) | —    | Cell‑type colour palette TSV                                              |
| `--show_ref`        |          | off             | Also render a reference UMAP on the left                                   |
| `--no_labels`       |          | off             | Hide the cell‑type text labels on the Beta UMAP                            |
| `--raster` (default)/`--no_raster` |  | raster on    | Rasterise the UMAP points for fast PDF rendering                          |
| `--png`             |          | off             | Save PNG instead of PDF                                                    |
| `--pt_size`         |          | `0.25`          | UMAP point size                                                            |
| `--max_cells`       |          | `250000`        | Downsample UMAP to this many cells                                         |
| `--out`             |          | `umap_plot`     | Output filename (the extension is set automatically)                       |

---

## Example runs

### Minimal – LocusZoom only, no UMAP / master / modules

Only the LocusZoom CSV is mandatory. The script auto-derives the module name from the filename and uses the top-P SNP as the anchor.

```bash
Rscript plot_umap_simplified_multimodules_one_side.R \
  --lz_files     UMAP_lz_values.csv \
  --ld_files     UMAP_lz_values.raw \
  --gtf          gencode.v49.annotation.gtf \
  --annotations  homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329.promoter.gff \
  --zoom_window  5000 \
  --name         GH \
  --out          umap_minimal.pdf --raster
```

### Beta UMAP only (no LocusZoom)

```bash
Rscript plot_umap_simplified_multimodules_one_side.R \
  --umap    UMAP_GH_QCed_f3_ct2.tsv \
  --master  freeze3_gh_cs_full_annotations.20260306.tsv \
  --modules M_18448 \
  --name    GH \
  --anno_join_col cell \
  --out     umap_GH_min.pdf \
  --raster
```

### Multiple modules – one row each

```bash
Rscript plot_umap_simplified_multimodules_one_side.R \
  --umap    UMAP_GH_QCed_f3_ct2.tsv \
  --master  freeze3_gh_cs_full_annotations.20260306.tsv \
  --modules "M_18448,M_11607,M_10171" \
  --name    GH --anno_join_col cell \
  --z_files "icd10_j15_zscore.csv,NA,bmi_zscore.csv" \
  --lz_files "j15_lz_ENSG00000101546.csv,NA,NA" \
  --gtf     gencode.v49.annotation.gtf \
  --out     umap_GH_multi.pdf --raster
```

Each comma‑separated list must be the same length as `--modules`; use `NA`
for any module that has no Z‑score or LocusZoom file.

### Restrict to one (cell, gene) pair of a multi‑gene module

```bash
Rscript plot_umap_simplified_multimodules_one_side.R \
  --umap    UMAP_GH_QCed_f3_ct2.tsv \
  --master  freeze3_gh_cs_full_annotations.20260306.tsv \
  --modules M_18448 \
  --name    GH --anno_join_col cell \
  --cell    T_CD4_CM \
  --gene    "RBFA,ENSG00000101544" \
  --lz_files j15_locuszoom_ENSG00000101546.csv \
  --gtf     gencode.v49.annotation.gtf \
  --out     umap_GH_one_pair.pdf --raster
```

### Full pipeline (the canonical example)

```bash
Rscript plot_umap_simplified_multimodules_one_side.R \
  --umap          UMAP_GH_QCed_f3_ct2.tsv \
  --master        freeze3_gh_cs_full_annotations.20260306.tsv \
  --modules       M_18448 \
  --name          GH \
  --anno_join_col cell \
  --z_files       icd10_j15_zscore_merge_notnull.csv \
  --lz_files      j15_locuszoom_ENSG00000101546.csv \
  --gtf           gencode.v49.annotation.gtf \
  --gene_list     gene_list_j15.csv \
  --annotations   homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329.promoter.gff \
  --zoom_window   5000 \
  --out           umap_GH_J15.pdf \
  --raster
```

---

## Output layout

- One row of the grid per **(module, cell, gene)** triple (after any `--cell` / `--gene` filter).
- Number of columns per row is `1 + has(z_files) + has(lz_files)` (1–3 columns). Columns you don't supply are not drawn at all (no empty slots left behind).
- Each LocusZoom cell is itself a stack of sub‑panels:
  1. `-log10(P)` scatter with the highlighted lead SNP,
  2. gene‑body track (genes stacked in lanes, focal eGene in red; filtered by `--gene_list` when supplied),
  3. zoom annotation panel (only if `--annotations` is supplied). Discrete tracks draw solid blocks; continuous tracks (BED with a numeric 5th column) draw score‑scaled bars.

Sizing scales automatically: page width = `n_cols × 6"`, page height = `n_rows × 4.5"` (minimum `5"`), where `n_rows` = number of surviving `(module, cell, gene)` triples.

---

## Tips

- If the run is slow because of very large UMAPs, lower `--max_cells` (default `250 000`) or keep rasterisation on (default).
- Passing an already‑built `*.coding_genes.tsv` to `--gtf` is auto‑detected by its header (`chrom start end strand gene_name gene_id`) — no awk pre‑filter runs.
- To change the order of annotation lanes, concatenate the files in the order you want in `--annotations`; the first‑seen feature types end up at the bottom of the zoom panel.
- Continuous annotation scores are min‑max normalised **within each feature type in the current window**, so profiles fill their lane. If you need absolute values across modules, pre‑normalise the scores in your BED file.

---

## Generating the inputs with `extract_z_lz.py` (optional)

The R script consumes ready-made `--lz_files`, `--ld_files`, and (optionally) `--z_files` CSVs. The companion Python script `extract_z_lz.py` pulls those tables straight out of an anndata / TileDB QTL store so you don't have to wire them up by hand. It is entirely optional — you can feed the R script anything with the expected columns.

### What it produces

Given one QTL module, the script writes:

| file                   | purpose                                                          | where it plugs in     |
|------------------------|------------------------------------------------------------------|-----------------------|
| `<out>.csv`            | headerless `CS, CHR, CELL, GENE, POS, P` rows (the LZ format)    | `--lz_files`          |
| `<out>.raw`            | plink2 `--export A include-alt` genotype dosages for every SNP  | `--ld_files`          |
| `<out>_zscores.csv`    | `snp, z_disease, z_qtl` coloc table (only if `--dis_adata` set) | `--z_files`           |

One call emits one set of files per module; the R script then accepts a comma-separated list of these (one per module) in the corresponding CLI flags.

### Inputs

| arg                  | what                                                               |
|----------------------|--------------------------------------------------------------------|
| `--qtl_module_adata` | H5AD with QTL modules (`.obs` has `phenotype_id` and `list_of_cs`) |
| `--qtl_cs_adata`     | H5AD with per-CS SNP tables (`.var` has `chr`, `pos`; `.obs` has `chr`, `start`, `end`) |
| `--tiledb`           | TileDB array with raw QTL P-values (queried by `(chr, cell, gene, start:end)`) |
| `--qtl_module`       | Module id to extract (e.g. `M_18448`, `M_3934`)                   |
| `--dis_adata`        | *(optional)* H5AD with disease fine-mapping, triggers the Z-score CSV |
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
`UMAP_lz_values_zscores.csv` that you can then feed straight into:

```bash
Rscript plot_umap_simplified_multimodules_one_side.R \
  --lz_files UMAP_lz_values.csv \
  --ld_files UMAP_lz_values.raw \
  --z_files  UMAP_lz_values_zscores.csv \
  --out      umap_M18448.pdf --raster
```

### Python dependencies
You can use on the hpc the tdbsumstat conda environment or you need the following packages
```bash
pip install scanpy pandas tiledb scipy numpy
```

`extract_z_lz.py` also shells out to `plink2` (must be on `$PATH`) to produce the `.raw` genotype matrix from the SNP list.
