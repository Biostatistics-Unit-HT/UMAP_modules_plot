# Program to extract z-scores and p-values from an anndata

import scanpy as sc
import pandas as pd
import argparse
import tiledb
from scipy.stats import norm
import numpy as np
import subprocess
import os

parser = argparse.ArgumentParser(
    prog='extract_z-scores',
    description='This program extracts z-scores from anndata'
)

parser.add_argument("--dis_adata", default=False, help="Anndata containing the credible sets for the diseases")
parser.add_argument("--qtl_module_adata",  help="Anndata containing the modules for the QTLs")
parser.add_argument("--qtl_cs_adata",  help="Anndata containing the credible set for the QTLs")
parser.add_argument("--tiledb", help="TileDB where raw QTLs are stored")
parser.add_argument("--dis_cs", default=False, help="Name of the credible set of the disease")
parser.add_argument("--qtl_module", help="Name of the credible set (module) of the qtl")
parser.add_argument("--safeld", default=False, help="Path of SAFE-ld")
parser.add_argument("--out", help="Name of the output file")
args = parser.parse_args()

qtl_adata_module = sc.read_h5ad(args.qtl_module_adata)
qtl_adata_cs = sc.read_h5ad(args.qtl_cs_adata)

tdb = tiledb.open(args.tiledb)

# FIX: Simplified redundant indexing
qtl_adata_obs = qtl_adata_module[qtl_adata_module.obs_names == args.qtl_module, :].copy()
list_pheno = qtl_adata_obs.obs["phenotype_id"].iloc[0].split(", ")
list_cs = qtl_adata_obs.obs["list_of_cs"].iloc[0].split(", ")
chrom = qtl_adata_obs.obs['chr'].iloc[0]

if args.dis_adata and args.dis_cs:
    disease_adata = sc.read_h5ad(args.dis_adata)
    mask_obs = disease_adata.obs["cs_name"] == args.dis_cs
    disease_adata_obs = disease_adata[mask_obs, :].copy()
    start = disease_adata_obs.obs['start'].iloc[0]
    end = disease_adata_obs.obs['end'].iloc[0]
    mask_var = (disease_adata_obs.var["chr"] == chrom) & (disease_adata_obs.var["pos"] > start) & (disease_adata_obs.var["pos"] < end)
    disease_adata_obs_var = disease_adata_obs[:, mask_var]
    
    # FIX: Safe division to prevent 0-division errors
    beta_dis = disease_adata_obs_var.layers['beta'].toarray()[0]
    se_dis = disease_adata_obs_var.layers['se'].toarray()[0]
    with np.errstate(divide='ignore', invalid='ignore'):
        z_disease = np.true_divide(beta_dis, se_dis)
        z_disease[~np.isfinite(z_disease)] = np.nan # Converts inf/zero-division results to NaN

    disease_adata_zscore = pd.DataFrame({
        "snp": disease_adata_obs_var.var['snp'],
        "z_disease": z_disease
    })

snps_plink = pd.DataFrame()

for idx, (pheno, cs) in enumerate(zip(list_pheno, list_cs)):
    cell, gene = pheno.split(":")
    qtl_adata_cs_obs = qtl_adata_cs[cs, :].copy()
    start = qtl_adata_cs_obs.obs['start'].iloc[0]
    end = qtl_adata_cs_obs.obs['end'].iloc[0]
    mask_var = (qtl_adata_cs_obs.var["chr"] == chrom) & (qtl_adata_cs_obs.var["pos"] > start) & (qtl_adata_cs_obs.var["pos"] < end)
    qtl_adata_cs_obs_var = qtl_adata_cs_obs[:, mask_var]
    
    # FIX: Query TileDB only once to save time
    chrom_num = int(chrom.split("chr")[1])
    gene_extracted = tdb.query(attrs=["SNPID", "P"]).df[chrom_num, cell, gene, start:end]
    
    qtl_adata_lz = pd.DataFrame({
        "CS": cs, "CHR": chrom, "CELL": cell, "GENE": gene,  
        "POS": gene_extracted["POS"], "P": gene_extracted["P"]
    })
    
    # FIX: Conditionally write headers only on the very first loop iteration
    write_header = (idx == 0)
    
    qtl_adata_lz.to_csv(f"{args.out}.csv", index=False, mode='a', header=write_header)
    gene_extracted.to_csv(f"{args.out}_genelist.csv", index=False, mode='a', header=write_header)
    
    # FIX: Pass as a DataFrame (double brackets) to prevent column naming issues
    snps_plink = pd.concat([snps_plink, gene_extracted[["SNPID"]]])
    
    if args.dis_adata and args.dis_cs:
        mask_var_modules = (qtl_adata_cs_obs.var["chr"] == chrom) & (qtl_adata_cs_obs.var["pos"] > start) & (qtl_adata_cs_obs.var["pos"] < end)
        qtl_adata_modules_obs_var = qtl_adata_cs_obs[:, mask_var_modules]
        
        # FIX: Safe division for QTL as well
        beta_qtl = qtl_adata_modules_obs_var.layers['beta'].toarray()[0]
        se_qtl = qtl_adata_modules_obs_var.layers['se'].toarray()[0]
        with np.errstate(divide='ignore', invalid='ignore'):
            z_qtl = np.true_divide(beta_qtl, se_qtl)
            z_qtl[~np.isfinite(z_qtl)] = np.nan

        qtl_adata_zscore = pd.DataFrame({
            "cs_qtl": cs,  
            "snp": qtl_adata_modules_obs_var.var['snp'].values,
            "z_qtl": z_qtl
        })
        
        disease_qtl_zscore_merge = disease_adata_zscore.merge(qtl_adata_zscore, on="snp")
        
        # FIX: Drop NaNs across both columns simultaneously
        icd10_zscore_merge_notnull = disease_qtl_zscore_merge.dropna(subset=["z_disease", "z_qtl"])
        
        icd10_zscore_merge_notnull.to_csv(f"{args.out}_zscores.csv", index=False, mode="a", header=write_header)

snps_plink = snps_plink.drop_duplicates()
snps_plink["SNPID"] = snps_plink["SNPID"].str.replace("chr", "")

if args.safeld:
    # FIX: Make the temp file unique to this specific parallel job using args.out
    tmp_file = f"{args.out}_tmp_snp_plink.txt"
    snps_plink.to_csv(tmp_file, index=False, header=None)
    
    cmd = [
        'plink2', 
        '--pfile', f'{args.safeld}', 
        '--extract', tmp_file, 
        '--export', 'A', 'include-alt', 
        '--out', f'{args.out}'
    ]
    subprocess.run(cmd, check=True)
    
    # Clean up the unique temporary plink text file
    if os.path.exists(tmp_file):
        os.remove(tmp_file)
