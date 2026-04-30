#Program to extract z-score and p-vlaues from an anndata

import scanpy as sc
import pandas as pd
import argparse
import tiledb
from scipy.stats import norm
import numpy as np
import subprocess

parser = argparse.ArgumentParser(
                    prog='extract_z-scores',
                    description='This program extract z-scores from anndata'
)

parser.add_argument("--dis_adata", default=False, help="Anndata containing the credible sets for the diseases")
parser.add_argument("--qtl_module_adata",  help="Anndata containing the modules for the QTLs")
parser.add_argument("--qtl_cs_adata",  help="Anndata containing the credible set for the QTLs")
parser.add_argument("--tiledb", help = "TileDB were raw QTLs are stored")
parser.add_argument("--dis_cs", default=False, help="Name of the credible set of the disease")
parser.add_argument("--qtl_module", help="Name of the credible set (module) of the qtl")
parser.add_argument("--safeld", default = False, help="Path of SAFE-ld")
parser.add_argument("--out", help = "Name of the output file")
args = parser.parse_args()


#this program was tested on the following parameters
# /ssu/bsssu/anndata_finemapping_repository/GH:gwas:ICD10codes_anndata.h5ad
# /ssu/bsssu/anndata_finemapping_repository/GH_meta_cardinal_F3_qtl_cis_07_01_hypr_modules_modules_JF_PIP_isCS99.h5ad
# /ssu/bsssu/anndata_finemapping_repository/GH_meta_cardinal_F3_qtl_cis_07_01_anndata_PIP_isCS99.h5ad
# /project/cardinal/QTLs/TileDB_GH_f3_05_1_26/TileDB_tiledb_gh_meta_celltype2_f3_05_01_26
# chr18::GH:gwas:ICD10codes::J15::chr18:80004668:A:T::L1
# M_18448
# /project/cardinal/QTLs/TileDB_GH_f3_05_1_26/TileDB_tiledb_gh_meta_celltype2_f3_05_01_26
# M_3934
#/project/cardinal/safeld_storage/GH_5k/final_output

qtl_adata_module = sc.read_h5ad(args.qtl_module_adata)
qtl_adata_cs = sc.read_h5ad(args.qtl_cs_adata)

tdb = tiledb.open(args.tiledb)


selected_ids = qtl_adata_module.obs_names[qtl_adata_module.obs.index == args.qtl_module]
qtl_adata_obs = qtl_adata_module[qtl_adata_module.obs_names.isin(selected_ids), :].copy()
list_pheno = qtl_adata_obs.obs['phenotype_id'][0].split(", ")
list_cs = qtl_adata_obs.obs['list_of_cs'][0].split(", ")
chrom = qtl_adata_obs.obs['chr'].iloc[0]

if (args.dis_adata!=False & args.dis_cs!=False):
    disease_adata = sc.read_h5ad(args.dis_adata)
    mask_obs = disease_adata.obs["cs_name"]==args.dis_cs
    disease_adata_obs = disease_adata[mask_obs,:].copy()
    start = disease_adata_obs.obs['start'].iloc[0]
    end = disease_adata_obs.obs['end'].iloc[0]
    mask_var = (disease_adata_obs.var["chr"]==chrom) & (disease_adata_obs.var["pos"] > start) & (disease_adata_obs.var["pos"] < end)
    disease_adata_obs_var = disease_adata_obs[:, mask_var]
    disease_adata_zscore = pd.DataFrame({"snp":disease_adata_obs_var.var['snp'],"z_disease":disease_adata_obs_var.layers['beta'].toarray()[0]/disease_adata_obs_var.layers['se'].toarray()[0]})

snps_plink = pd.DataFrame()
for pheno,cs in zip(list_pheno,list_cs):
    cell,gene = pheno.split(":")
    qtl_adata_cs_obs = qtl_adata_cs[cs,:].copy()
    start = qtl_adata_cs_obs.obs['start'].iloc[0]
    end = qtl_adata_cs_obs.obs['end'].iloc[0]
    mask_var = (qtl_adata_cs_obs.var["chr"]==chrom) & (qtl_adata_cs_obs.var["pos"] > start) & (qtl_adata_cs_obs.var["pos"] < end)
    qtl_adata_cs_obs_var = qtl_adata_cs_obs[:, mask_var]
    #beta = qtl_adata_cs_obs_var.layers['beta'].toarray()[0]
    #se = qtl_adata_cs_obs_var.layers['se'].toarray()[0]
    #with np.errstate(divide='ignore', invalid='ignore'):
    #    z_scores = np.true_divide(beta, se)
    #    z_scores[~np.isfinite(z_scores)] = 0  # Replace NaNs or Infs with 0 or NaN as preferred
    #abs_z = np.abs(z_scores)
    #mlog10p = -(np.log10(2) + norm.logsf(abs_z) / np.log(10))
    gene_extracted_to_plot = tdb.query(attrs = ["SNPID","P"]).df[int(chrom.split("chr")[1]), cell, gene, start:end]
    gene_extracted = tdb.query(attrs = ["SNPID","P"]).df[int(chrom.split("chr")[1]), cell, gene, start:end]
    qtl_adata_lz = pd.DataFrame({"CS":cs ,"CHR": chrom, "CELL":cell, "GENE":gene,  "POS":gene_extracted["POS"],"P":gene_extracted["P"]})
    qtl_adata_lz.to_csv(f"{args.out}.csv", index = False, mode = 'a', header=None)
    gene_extracted_to_plot.to_csv(f"{args.out}_genelist.csv", index = False, mode = 'a', header=None)
    snps_plink = pd.concat([snps_plink, gene_extracted["SNPID"]])
    if (args.dis_adata!=False & args.dis_cs!=False):
        mask_var_modules = (qtl_adata_cs_obs.var["chr"]==chrom) & (qtl_adata_cs_obs.var["pos"] > start) & (qtl_adata_cs_obs.var["pos"] < end)
        qtl_adata_modules_obs_var = qtl_adata_cs_obs[:, mask_var_modules]
        qtl_adata_zscore = pd.DataFrame({"cs_qtl":qtl_adata_modules_obs_var.obs['cs_name'],"snp":qtl_adata_modules_obs_var.var['snp'],"z_qtl":qtl_adata_modules_obs_var.layers['beta'].toarray()[0]/qtl_adata_modules_obs_var.layers['se'].toarray()[0]})
        disease_qtl_zscore_merge = disease_adata_zscore.merge(qtl_adata_zscore, on = "snp")
        icd10_zscore_merge_notnull= disease_qtl_zscore_merge[~disease_qtl_zscore_merge["z_disease"].isna()]
        icd10_zscore_merge_notnull= disease_qtl_zscore_merge[~disease_qtl_zscore_merge["z_qtl"].isna()]
        icd10_zscore_merge_notnull.to_csv(f"{args.out}_zscores.csv", index = False, mode = "a")

snps_plink = snps_plink.drop_duplicates()
snps_plink = snps_plink.drop("cs_qtl")
snps_plink["SNPID"] = snps_plink["SNPID"].str.replace("chr","")


# 4. Run the command
# check=True ensures Python throws an error if the PLINK command fails
if args.safeld:
    # 2. Save to a text file (one SNP per line, no headers or row numbers)
    snps_plink.to_csv("tmp_snp_plink.txt", index=False, header=None)
    # 3. Build the PLINK2 command as a list
    cmd = [
    'plink2', 
    '--pfile', f'{args.safeld}', 
    '--extract', 'tmp_snp_plink.txt', 
    '--export', 'A', 'include-alt', # Each word is a separate string
    '--out', f'{args.out}'
    ]
    subprocess.run(cmd, check=True)

