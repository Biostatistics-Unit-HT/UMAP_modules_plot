from dataclasses import dataclass
import tiledb
import anndata as ad
import scanpy as sc
import subprocess
import os
import pandas as pd 
import numpy as np 

@dataclass
class GenomicRegion:
    chrom_num: int
    cell:str
    gene:str
    start: int
    end: int
    credset: str

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
    chrom_str = adata.obs["chr"].iloc[0]
    chrom_num = int(chrom_str.split("chr")[1])
    start = adata.obs["start"].iloc[0]
    end = adata.obs["end"].iloc[0]
    cell, gene = adata.obs["phenotype_id"].iloc[0].split(":")
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
    adata.var["pos"] = adata.var["pos"].astype(int)
    chrom_str = f"chr{chrom}"
    mask_var = (
        (adata.var["chr"] == chrom_str)
        & (adata.var["pos"] > start)
        & (adata.var["pos"] < end)
    )
    disease_adata_obs_var = adata[:, mask_var]
    # FIX: Safe division to prevent 0-division errors
    beta_dis = disease_adata_obs_var.layers["beta"].toarray()[0]
    se_dis = disease_adata_obs_var.layers["se"].toarray()[0]
    with np.errstate(divide="ignore", invalid="ignore"):
        #z_disease = np.true_divide(beta_dis, se_dis)
        z_disease = beta_dis
        #z_disease[~np.isfinite(z_disease)] = (
        #    np.nan
        #)  # Converts inf/zero-division results to NaN
    disease_adata_zscore = pd.DataFrame(
        {"SNPID": disease_adata_obs_var.var["snp"], "z_disease": z_disease}
    )
    return disease_adata_zscore
