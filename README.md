# UMAP + LocusZoom + Coloc plots

`plot_umap_simplified_multimodules_one_side.R` produces a composite figure for
one or more fine‑mapped credible sets ("modules") in a single population.
Each module contributes **one row** of the final grid with up to three panels:

```
[ LocusZoom ]   →   [ Beta UMAP ]   →   [ Coloc Z‑scores ]
```

- **LocusZoom** – `-log10(P)` vs genomic position, with
  - the fine‑mapped lead SNP highlighted as a purple diamond;
  - a gene‑body track of every *protein‑coding* gene in the window
    (the module's focal eGene in red, others in navy);
  - an optional zoomed‑in annotation panel around the lead SNP.
- **Beta UMAP** – the population's UMAP coloured by the per‑cell‑type effect
  size for the module.
- **Coloc Z‑scores** – QTL‑Z vs disease‑Z scatter for the colocalised SNPs.

The LocusZoom and Coloc panels are optional (enabled by passing the
corresponding files).

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
  --summary_table Nicole_summary_table.csv \
  --gtf           gencode.v49.annotation.gtf \
  --annotations   homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329.promoter.gff \
  --zoom_window   5000 \
  --out           umap_GH_J15.pdf \
  --raster
```

The output is a multi‑page vector PDF by default (`--png` switches to PNG).

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

## Required input files

### 1. `--umap` – the UMAP coordinates table (TSV)

One row per cell, tab separated. Expected columns:

| column      | description                                       |
|-------------|---------------------------------------------------|
| `UMAP_1`    | first UMAP coordinate                             |
| `UMAP_2`    | second UMAP coordinate                            |
| `<celltype>`| a cell‑type label column (name via `--join_col`)  |

If the file ships with `UMAP_0` / `UMAP_1` instead, they are renamed
internally to `UMAP_1` / `UMAP_2`.

The default cell‑type column is `celltype_2`; change it with
`--join_col <column_name>`.

### 2. `--master` – master credible‑set annotations (TSV)

Usually something like `freeze3_gh_cs_full_annotations.*.tsv`. One row per
(module, cell type) with fine‑mapping info. The script reads any of the
following columns that are present:

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

### 3. `--modules` – module IDs

Comma‑separated list, one entry per module you want plotted:
`--modules "M_18448,M_11607,M_10171"`. Each module produces one grid row.

---

## Optional side panels

### 4. `--z_files` – Coloc Z‑score CSVs  *(enables the coloc panel)*

Comma‑separated list, **same length and order as `--modules`**, one CSV per
module. Use the literal string `NA` for a module that has no Z‑score file.

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

Comma‑separated list, **same length and order as `--modules`**. Use `NA`
to skip a module. Expected columns (header **must** match):

| column | description                                    |
|--------|------------------------------------------------|
| `CHR`  | chromosome (plain `18` or `chr18`)             |
| `CELL` | cell type the association was computed in      |
| `GENE` | Ensembl gene ID (e.g. `ENSG00000101546`)       |
| `POS`  | genomic position (bp, 1‑based)                 |
| `P`    | association P‑value                            |

Example:

```csv
CHR,CELL,GENE,POS,P
18,T_CD4_CM,ENSG00000101546,79869986,0.254
18,T_CD4_CM,ENSG00000101546,79870170,0.775
...
```

### 6. `--summary_table` – per‑module disease summary (CSV)

Optional. Used only to enrich the Beta panel title with the colocalised
disease and its β / P. One row per module with any of:

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

When supplied, the LocusZoom gene track shows **every protein‑coding gene**
overlapping the plotted window (stacked into lanes to avoid overlap), with
the module's focal eGene coloured red. Without `--gtf`, only the focal
eGene from the master annotation is drawn.

### 8. `--annotations` – regulatory tracks for the zoom panel *(optional)*

Comma‑separated list of annotation files. Each becomes extra lanes in the
zoom panel underneath the gene track.

Two formats are auto‑detected:

- **Simple 4‑column TSV** (no header) — the format used by
  `homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.*.promoter.gff`:

  ```
  chrom  feature_type              start     end
  18     promoter_flanking_region  80009802  80011199
  ```

- **Standard 9‑column GFF/GTF** — `#` comment lines are skipped, column 1 is
  chromosome, column 3 is feature type, columns 4/5 are start/end.

Every distinct `feature_type` across all the supplied files becomes its own
lane (row) in the zoom panel, with a consistent lane order across modules.

`--zoom_window N` (default `5000`) controls the half‑width in bp of the
zoom panel around the lead SNP.

### 9. `--colors` – cell‑type colour palette (TSV) *(only with `--show_ref`)*

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
| `--master`          | ✔        | —               | Master credible‑set annotations TSV                                       |
| `--modules`         | ✔        | —               | Comma‑separated module IDs                                                |
| `--name`            |          | `Pop`           | Display name used in panel titles                                         |
| `--join_col`        |          | `celltype_2`    | Cell‑type column name in the UMAP                                         |
| `--anno_join_col`   |          | `--join_col`    | Cell‑type column name in the master annotations (if different)           |
| `--gene_col`        |          | `eGene_symbol`  | Gene‑symbol column in the master annotations                              |
| `--z_files`         |          | —               | Coloc Z‑score CSVs (comma‑separated, one per module, `NA` to skip)        |
| `--lz_files`        |          | —               | LocusZoom P‑value CSVs (comma‑separated, one per module, `NA` to skip)    |
| `--summary_table`   |          | —               | Per‑module disease summary CSV                                            |
| `--gtf`             |          | —               | Raw GTF or pre‑built `.coding_genes.tsv`                                  |
| `--annotations`     |          | —               | Comma‑separated annotation files (4‑col or 9‑col GFF)                     |
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

### Minimal – beta UMAP only

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
  --summary_table Nicole_summary_table.csv \
  --gtf     gencode.v49.annotation.gtf \
  --out     umap_GH_multi.pdf --raster
```

Each comma‑separated list must be the same length as `--modules`; use `NA`
for any module that has no Z‑score or LocusZoom file.

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
  --summary_table Nicole_summary_table.csv \
  --gtf           gencode.v49.annotation.gtf \
  --annotations   homo_sapiens.GRCh38.Regulatory_Build.regulatory_features.20190329.promoter.gff \
  --zoom_window   5000 \
  --out           umap_GH_J15.pdf \
  --raster
```

---

## Output layout

- One row of the grid per module.
- Number of columns per row is `1 + has(z_files) + has(lz_files)`
  (1–3 columns).
- Each LocusZoom panel is itself a stack of sub‑panels:
  1. `-log10(P)` scatter with the highlighted lead SNP,
  2. gene‑body track (genes stacked in lanes, focal eGene in red),
  3. zoom annotation panel (only if `--annotations` is supplied).

Sizing scales automatically: page width = `n_cols × 6"`, page height =
`n_modules × 4.5"` (minimum `5"`).

---

## Tips

- If the run is slow because of very large UMAPs, lower `--max_cells`
  (default `250 000`) or keep rasterisation on (default).
- The zero‑byte cache problem you can hit from a crashed GTF run is handled
  — passing an already‑built `*.coding_genes.tsv` to `--gtf` is auto‑detected.
- To change the order of annotation lanes, concatenate the files in the
  order you want in `--annotations`; the first‑seen feature types end up
  at the bottom of the zoom panel.
