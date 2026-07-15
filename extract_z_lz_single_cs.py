# Program to extract effect sizes (BETA) and p-values from a single credible set.
import subprocess
import sys

import pandas as pd
import click
import cloup
from utils.functions import parse_adata, extract_region, query_tiledb, calculate_ld, extract_beta_adata


def _tabix_region(gz_path, chrom, start, end):
    """Yield raw lines from a bgzipped+tabix-indexed file for chrom:start-end.

    Tries pysam first; falls back to the `tabix` binary. The contig name in the
    index may or may not carry a 'chr' prefix, so both are attempted.
    """
    chrom = str(chrom).replace("chr", "")
    candidates = [chrom, f"chr{chrom}"]

    try:
        import pysam

        with pysam.TabixFile(gz_path) as tbx:
            available = set(tbx.contigs)
            for c in candidates:
                if c in available:
                    yield from tbx.fetch(c, int(start), int(end))
                    return
        raise ValueError(
            f"Neither '{chrom}' nor 'chr{chrom}' is a contig in {gz_path}. "
            f"Contigs look like: {sorted(available)[:5]}"
        )
    except ImportError:
        pass

    last_err = None
    for c in candidates:
        try:
            res = subprocess.run(
                ["tabix", gz_path, f"{c}:{int(start)}-{int(end)}"],
                capture_output=True, text=True, check=True,
            )
            if res.stdout.strip():
                yield from res.stdout.splitlines()
                return
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            last_err = e
    raise RuntimeError(
        f"tabix query failed for {gz_path} at {chrom}:{start}-{end} "
        f"(pysam unavailable, `tabix` binary error: {last_err})"
    )


def extract_beta_gz(gz_path, chrom, start, end):
    """Disease BETA (effect size) from a tabix-indexed sumstats file.

    Expected columns (tab-separated, header starting '#chrom'):
        chrom pos ref alt rsids nearest_genes pval mlogp beta sebeta af_alt ...

    Returns a DataFrame [SNPID, beta_disease] where SNPID is
    'chr<c>:<pos>:<ref>:<alt>' to match the TileDB QTL SNPIDs. Both allele
    orientations are emitted (the flipped one with a negated BETA) so the
    downstream merge matches regardless of which strand/orientation the
    sumstats used.

    Example:
        extract_beta_gz("gwas.tsv.gz", chrom=22, start=17600000, end=17610000)
        # ->            SNPID  beta_disease
        #   0  chr22:17604981:C:T          0.12
        #   1  chr22:17604981:T:C         -0.12
    """
    rows = []
    for line in _tabix_region(gz_path, chrom, start, end):
        if line.startswith("#"):
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 10:
            continue
        c, pos, ref, alt = f[0].replace("chr", ""), f[1], f[2], f[3]
        try:
            beta = float(f[8])
        except ValueError:
            continue
        if pd.isna(beta):
            continue
        
        rows.append((f"chr{c}:{pos}:{ref}:{alt}", beta))
        rows.append((f"chr{c}:{pos}:{alt}:{ref}", -beta))

    if not rows:
        raise ValueError(
            f"No usable rows in {gz_path} for region {chrom}:{start}-{end}"
        )

    df = pd.DataFrame(rows, columns=["SNPID", "beta_disease"])
    # A palindromic SNP can produce the same SNPID twice; keep the first.
    return df.drop_duplicates(subset=["SNPID"], keep="first")


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
    default=None,
    help="Bgzipped + tabix-indexed disease summary statistics (.gz with a .tbi "
         "alongside). Columns: chrom pos ref alt rsids nearest_genes pval mlogp "
         "beta sebeta af_alt ... Takes precedence over --dis_adata when both are "
         "given. The disease effect size used is the raw BETA.",
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
        chrom_num=coord.chrom_num, start=coord.start, end=coord.end,
        cell=coord.cell, gene=coord.gene, tiledb_path=tiledb_path
    )
    pd_cs_qtl["CS"] = coord.credset
    # mode="w": re-running a job must overwrite, not append a second header.
    pd_cs_qtl.to_csv(f"{out}.csv", index=False, mode="w", header=True)

    # Make LD from Plink
    snps_plink = pd_cs_qtl[["SNPID"]].drop_duplicates()
    snps_plink["SNPID"] = snps_plink["SNPID"].str.replace("chr", "")
    calculate_ld(df_snps=snps_plink, out_file=out, plink_file=ld_file)

    if not (dis_gz or dis_adata):
        print("No --dis_gz or --dis_adata given; skipping beta extraction.")
        return

    qtl_adata_lz = pd_cs_qtl.copy()
    qtl_adata_lz["beta_qtl"] = qtl_adata_lz["BETA"]

    if dis_gz:
        print(f"Disease beta from tabix sumstats: {dis_gz}")
        disease_beta = extract_beta_gz(
            dis_gz, coord.chrom_num, coord.start, coord.end
        )
    else:
        print(f"Disease beta from anndata: {dis_adata} (cs={dis_cs})")
        adata_dis = parse_adata(adata=dis_adata, cs_name=dis_cs)
        disease_beta = extract_beta_adata(adata=adata_dis, chrom=coord.chrom_num)

    merged = disease_beta.merge(qtl_adata_lz, on="SNPID")
    print(
        f"Region {coord.chrom_num}:{coord.start}-{coord.end} | "
        f"QTL SNPs={len(qtl_adata_lz)} disease SNPs={len(disease_beta)} "
        f"merged={len(merged)}"
    )
    if merged.empty:
        print(f"  sample QTL SNPID:     {qtl_adata_lz['SNPID'].iloc[0]}")
        print(f"  sample disease SNPID: {disease_beta['SNPID'].iloc[0]}")
        sys.exit("ERROR: no SNPs shared between QTL and disease (ID format mismatch?)")

    keep = merged[
        merged["beta_disease"].notna() & (merged["beta_disease"] != 0)
        & merged["beta_qtl"].notna() & (merged["beta_qtl"] != 0)
    ].copy()
    print(f"  {len(keep)} SNPs with non-zero beta in both.")
    keep.to_csv(f"{out}_betas.csv", index=False, mode="w", header=True)


if __name__ == "__main__":
    cli()
