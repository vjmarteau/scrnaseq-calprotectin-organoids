#!/usr/bin/env python3

"""
Usage:
  Filter_mito_outliers.py --adata=<adata> [options]

Mandatory arguments:
  --adata=<adata>       adata

Optional arguments:
  --resDir=<resDir>     Output directory [default: ./]
"""

from docopt import docopt
import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd

args = docopt(__doc__)
adata = args["--adata"]
resDir = args["--resDir"]

adata = sc.read_h5ad(adata)

# Filter outliers
def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * M.mad()) | (
        np.median(M) + nmads * M.mad() < M
    )
    return outlier

sample_d = dict()
for s in adata.obs['sample'].values.unique():
    adata_f = adata[adata.obs["sample"] == s, :].copy()
    adata_f.obs["value"] = 0

    # mitochondrial genes
    adata_f.var["mt"] = adata_f.var_names.str.startswith("MT-")
    # ribosomal genes
    adata_f.var["ribo"] = adata_f.var_names.str.startswith(("RPS", "RPL"))
    # hemoglobin genes.
    adata_f.var["hb"] = adata_f.var_names.str.contains(("^HB[^(P)]"))

    sc.pp.calculate_qc_metrics(
        adata_f, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
    )

    adata_f.obs["outlier"] = (
        is_outlier(adata_f, "log1p_total_counts", 5)
        | is_outlier(adata_f, "log1p_n_genes_by_counts", 5)
        | is_outlier(adata_f, "pct_counts_in_top_20_genes", 5)
    )
    adata_f.obs.outlier.value_counts()

    adata_f.obs["mt_outlier"] = is_outlier(adata_f, "pct_counts_mt", 3) | (
        adata_f.obs["pct_counts_mt"] > 12
    )
    adata_f.obs.mt_outlier.value_counts()

    adata_f = adata_f[(~adata_f.obs.outlier) & (~adata_f.obs.mt_outlier)].copy()
    sample_d[s] = adata_f

# After filtering individual samples concat to make one adata object
adatas_new_l = []
for s in adata.obs['sample'].values.unique():
    adatas_new_l.append(sample_d[s])
adata_filtered = ad.concat(adatas_new_l)

adata = adata[adata_filtered.obs.index, adata_filtered.var.index]

# Save adata
adata.write(f"{resDir}/adata_filtered.h5ad", compression="gzip")
