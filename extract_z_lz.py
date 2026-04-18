#Program to extract z-score and p-vlaues from an anndata

import scanpy as sc
import pandas as pd
import argparse

parser = argparse.ArgumentParser(
                    prog='extract_z-scores',
                    description='This program extract z-scores from anndata'
)

parser.add_argument("--dis_adata", help="Anndata containing the credible sets for the diseases")
parser.add_argument("--qtl_adata", help="Anndata containing the credible set (or modules) for the QTLs")
parser.add_argument("--dis_cs", help="Name of the credible set of the disease")
parser.add_argument("--qtl_cs", help="Name of the credible set (module) of the qtl")
parser.add_argument("--out", help = "Name of the output file")
args = parser.parse_args()

#this program was tested on the following parameters
# /Users/bruno.ariano/work/HT/popspecific_analysis/GH_joint_ICD10_NK_j15_locus.h5ad
# /Users/bruno.ariano/work/HT/popspecific_analysis/GH_joint_Cardinal_NK_j15_locus.h5ad
#chr18::GH:gwas:ICD10codes::J15::chr18:80004668:A:T::L1
#M_18448


disease_adata = sc.read_h5ad(args.dis_adata)
qtl_adata = sc.read_h5ad(args.qtl_adata)

mask_obs = disease_adata.obs["cs_name"]==args.dis_cs

disease_adata_obs = disease_adata[mask_obs,:].copy()
start = disease_adata_obs.obs['start'].iloc[0]
end = disease_adata_obs.obs['end'].iloc[0]
chrom = args.dis_cs.split("::")[0]
mask_var = (disease_adata_obs.var["chr"]==chrom) & (disease_adata_obs.var["pos"] > start) & (disease_adata_obs.var["pos"] < end)
disease_adata_obs_var = disease_adata_obs[:, mask_var]


selected_ids = qtl_adata.obs_names[qtl_adata.obs.index == args.qtl_cs]
qtl_adata_obs = qtl_adata[qtl_adata.obs_names.isin(selected_ids), :].copy()
start = qtl_adata_obs.obs['start'].iloc[0]
end = qtl_adata_obs.obs['end'].iloc[0]
mask_var = (qtl_adata_obs.var["chr"]==chrom) & (qtl_adata_obs.var["pos"] > start) & (qtl_adata_obs.var["pos"] < end)
qtl_adata_obs_var = qtl_adata_obs[:, mask_var].copy()


disease_adata_zscore = pd.DataFrame({"snp":disease_adata_obs_var.var['snp'],"z_disease":disease_adata_obs_var.layers['beta'].toarray()[0]/disease_adata_obs_var.layers['se'].toarray()[0]})
qtl_adata_zscore = pd.DataFrame({"snp":qtl_adata_obs_var.var['snp'],"z_qtl":qtl_adata_obs_var.layers['beta'].toarray()[0]/qtl_adata_obs_var.layers['se'].toarray()[0]})
disease_qtl_zscore_merge = disease_adata_zscore.merge(qtl_adata_zscore, on = "snp")
icd10_j15_zscore_merge_notnull= disease_qtl_zscore_merge[~disease_qtl_zscore_merge["z_disease"].isna()]
icd10_j15_zscore_merge_notnull= disease_qtl_zscore_merge[~disease_qtl_zscore_merge["z_qtl"].isna()]

icd10_j15_zscore_merge_notnull.to_csv(args.out, index = False)