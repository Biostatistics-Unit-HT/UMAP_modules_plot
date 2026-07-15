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


def extract_beta_adata(adata: ad.AnnData, chrom: int) -> pd.DataFrame:
    """Extract the disease per-SNP BETA over the credible-set region from an anndata.

    The disease effect size used downstream is the raw BETA (not BETA/SE), so the
    coloc panel compares QTL beta vs disease beta.

    Example:
        # adata with layers["beta"] and var columns chr/pos/snp
        extract_beta_adata(adata, chrom=22)
        # ->            SNPID  beta_disease
        #   0  chr22:17604981:C:T          0.12
        #   1  chr22:17605100:G:A         -0.04
    """
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
    beta_dis = disease_adata_obs_var.layers["beta"].toarray()[0]
    disease_beta = pd.DataFrame(
        {"SNPID": disease_adata_obs_var.var["snp"], "beta_disease": beta_dis}
    )
    return disease_beta
