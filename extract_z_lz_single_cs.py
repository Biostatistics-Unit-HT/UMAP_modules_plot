# Program to extract z-scores and p-values from an anndata

import scanpy as sc
import pandas as pd
import click
import cloup
import tiledb
import subprocess
import anndata as ad
import os
from dataclasses import dataclass
import numpy as np


@dataclass
class GenomicRegion:
    chrom_num: int
    cell:str
    gene:str
    start: int
    end: int
    credset: str

@click.command()
@cloup.option(
    "--dis_adata",
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True),
    help="Path of the anndata for the disease",
)
@cloup.option(
    "--qtl_cs_adata",
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True),
    help="Path of the anndata for the QTL",
)
@cloup.option("--qtl_module", type=str, help="Name of the module for the QTL")
@cloup.option(
    "--dis_cs", type=str, default=None, help="Name of the credible set of the disease"
)
@cloup.option(
    "--tiledb_path",
    type=click.Path(exists=True, file_okay=False, dir_okay=True, readable=True),
    help="TileDB where raw QTLs are stored",
)
@cloup.option(
    "--dis_gz",
    type=click.Path(exists=True, file_okay=True, dir_okay=False, readable=True),
    help="gz file of summary statistics",
)
@cloup.option(
    "--ld_file", type=str, default=None, help="Path plink genomics file to calculate LD"
)
@cloup.option("--out", type=str, help="Prefix for the output files")
def cli(
    dis_adata,
    dis_cs,
    dis_gz,
    qtl_cs_adata,
    qtl_module,
    tiledb_path,
    ld_file,
    out,
):
    """Extract all the files necessary to plot UMAPs from single credible set"""
    adata = parse_adata(adata=qtl_cs_adata, cs_name=qtl_module)
    coord = extract_region(adata)
    pd_cs_qtl = query_tiledb(
        chrom_num=coord.chrom_num, start=coord.start, end=coord.end, cell=coord.cell, gene=coord.gene, tiledb_path=tiledb_path
    )
    pd_cs_qtl["CS"] = coord.credset
    pd_cs_qtl.to_csv(f"{out}.csv", index=False, mode="a", header=True)
    # Make LD from Plink
    snps_plink = pd_cs_qtl[["SNPID"]]
    snps_plink = snps_plink.drop_duplicates()
    snps_plink["SNPID"] = snps_plink["SNPID"].str.replace("chr", "")
    calculate_ld(df_snps=snps_plink, out_file=out, plink_file=ld_file)

    if dis_adata:
        qtl_adata_lz = pd_cs_qtl.copy()
        qtl_adata_lz["z_qtl"] = qtl_adata_lz["BETA"] / qtl_adata_lz["SE"]
        adata_dis = parse_adata(adata=dis_adata, cs_name=dis_cs)
        disease_adata_zscore = extract_z_adata(adata=adata_dis, chrom=coord.chrom_num)
        disease_qtl_zscore_merge = disease_adata_zscore.merge(qtl_adata_lz, on="SNPID")
        # FIX: Drop NaNs across both columns simultaneously
        icd10_zscore_merge_notnull = disease_qtl_zscore_merge.dropna(
            subset=["z_disease", "z_qtl"]
        )
        icd10_zscore_merge_notnull.to_csv(
            f"{out}_zscores.csv", index=False, mode="a", header=True
        )


def parse_adata(adata: str, cs_name: str) -> GenomicRegion:
    """Subsetting the QTL anndata for the module only"""
    adata_cs = sc.read_h5ad(adata)
    if cs_name not in adata_cs.obs_names:
        raise ValueError(f"The module '{cs_name}' is not a column in adata.obs.")
    adata_subset = adata_cs[adata_cs.obs_names == cs_name, :].copy()
    return adata_subset


def extract_region(adata: ad.AnnData) -> GenomicRegion:
    """
    Extract region information from anndata
    """
    cs = adata.obs.index[0]
    chrom_str = adata.obs["chr"][0]
    chrom_num = int(chrom_str.split("chr")[1])
    start = adata.obs["start"].iloc[0]
    end = adata.obs["end"].iloc[0]
    cell, gene = adata.obs["phenotype_id"][0].split(":")
    region = GenomicRegion(credset=cs, chrom_num=chrom_num, start=start, end=end, cell=cell, gene=gene)
    return region


def query_tiledb(
    chrom_num: int, start: int, end: int, cell: str, gene: str, tiledb_path: str
) -> pd.DataFrame:
    """Query TileDB to get the stats for the region selected"""
    tdb = tiledb.open(tiledb_path)
    gene_extracted = tdb.query(attrs=["SNPID", "P", "BETA", "SE"]).df[
        chrom_num, cell, gene, start:end
    ]
    return gene_extracted


def calculate_ld(df_snps: pd.DataFrame, out_file: str, plink_file: str):
    """Extract teh LD information from genotypes using Plink"""
    tmp_file = f"{out_file}_tmp_snp_plink.txt"
    df_snps.to_csv(tmp_file, index=False, header=None)
    cmd = [
        "/lustre/scratch124/humgen/projects_v2/cardinal_analysis/analysis/core_dataset/Manuel_beta_test/plink2",
        "--pfile",
        f"{plink_file}",
        "--extract",
        tmp_file,
        "--export",
        "A",
        "include-alt",
        "--out",
        f"{out_file}",
    ]
    subprocess.run(cmd, check=True)
    # Clean up the unique temporary plink text file
    if os.path.exists(tmp_file):
        os.remove(tmp_file)


def extract_z_adata(adata: ad.AnnData, chrom: int) -> pd.DataFrame:
    start = adata.obs["start"].iloc[0]
    end = adata.obs["end"].iloc[0]
    mask_var = (
        (adata.var["chr"] == chrom)
        & (adata.var["pos"] > start)
        & (adata.var["pos"] < end)
    )
    disease_adata_obs_var = adata[:, mask_var]
    # FIX: Safe division to prevent 0-division errors
    beta_dis = disease_adata_obs_var.layers["beta"].toarray()[0]
    se_dis = disease_adata_obs_var.layers["se"].toarray()[0]
    with np.errstate(divide="ignore", invalid="ignore"):
        z_disease = np.true_divide(beta_dis, se_dis)
        z_disease[~np.isfinite(z_disease)] = (
            np.nan
        )  # Converts inf/zero-division results to NaN
    disease_adata_zscore = pd.DataFrame(
        {"SNPID": disease_adata_obs_var.var["snp"], "z_disease": z_disease}
    )
    return disease_adata_zscore


if __name__ == "__main__":
    cli()



# --- Environment & Paths ---
# (Keeping these paths the same as your script)
#EXTRACT_SCRIPT="/work/UMAP_modules_plot/extract_z_lz_single_cs.py"
#PLOT_SCRIPT="/work/UMAP_modules_plot/plot_umap_simplified_multimodules_one_side.R"
##TILEDB_DIR="/lustre/scratch124/humgen/projects_v2/cardinal_analysis/analysed_datasets/eqtls/F3_QTLs_updated_12_2025/TileDB_ingestion/TileDB_UKBB/celltype_2/TileDB/TileDB_tiledb_ukbb_celltype2_f3_16_12_25"
#TILEDB_DIR="/lustre/scratch124/humgen/projects_v2/cardinal_analysis/analysed_datasets/eqtls/F3_QTLs_saige/UKBB/celltype_3/TileDB/TileDB/TileDB_tiledb_ukbb_celltype3_f3_11_03_26"
#TILEDB_DIS="/lustre/scratch124/humgen/projects_v2/cardinal_analysis/analysed_datasets/eqtls/F3_QTLs_saige/UKBB/celltype_3/TileDB/TileDB/"
#PLINK_DIR="/lustre/scratch124/humgen/projects_v2/cardinal_analysis/analysis/core_dataset/genotypes/UKB/pgen/maf0.001/hard_calls_all"
#ADATA_FOLDER="/lustre/scratch124/humgen/projects_v2/cardinal_analysis/analysis/al37/finemapping/UKBB_cardinal_ct3_CMV_finemapping/flanders_output/results/annda