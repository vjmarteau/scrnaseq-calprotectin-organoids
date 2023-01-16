#!/usr/bin/env python3

"""
Usage:
  Load_raw.py --samplesheet=<samplesheet> [options]

Mandatory arguments:
  --samplesheet=<samplesheet>       Samplesheet from nfcore/scrnaseq pipeline

Optional arguments:
  --resDir=<resDir>   Output directory [default: ./]
"""

from docopt import docopt
import anndata as ad
import scanpy as sc
import numpy as np
import pandas as pd
import scipy.io

args = docopt(__doc__)
samplesheet = args["--samplesheet"]
resDir = args["--resDir"]

# Get metadata from samplesheet
meta = pd.read_csv(samplesheet)
meta.drop(axis="columns", labels=["fastq_1", "fastq_2"], inplace=True)
#meta = pd.read_csv(meta, index_col=0)

# Reorder by Sample ID, drop double columns, and update index
sample_idx = []
for s in meta["sample"].values:
    sample_idx.append(int(s[2:]))

meta["sample_idx"] = sample_idx
meta = meta.sort_values("sample_idx")
meta.index = meta.sample_idx
meta = meta.drop(columns=["sample_idx"])

meta["sample_idx"] = meta["sample"]
meta = meta.set_index("sample_idx")
meta.index.name = None

# Drop patient P163
meta = meta[meta["patient"] != "P163"]

# load cellranger .h5 feature matrices
adatas_d = dict()
cnvan_key_l = []
for ind, sample in zip(meta.index, meta.to_dict(orient="records")):
    # nf finds the path because we give it ch_input_files as input - 
    # and the dot below completes to /data/projects/2022/Adolph-scRNA-organoids/01_nfcore_scrnaseq
    path_h5ads = f"./01_nfcore_scrnaseq/cellranger/sample-{ind}/outs/raw_feature_bc_matrix.h5"
    tmp_adata = sc.read_10x_h5(path_h5ads)

    # save gene conversion key and switch index to ensembl ids before making unique
    cnvan_key_l.append(tmp_adata.var.copy())
    tmp_adata.var = tmp_adata.var.drop(columns=["feature_types", "genome"])
    tmp_adata.var_names_make_unique()
    assert tmp_adata.obs_names.is_unique
    tmp_adata.obs = tmp_adata.obs.assign(**sample)
    adatas_d[sample["sample"]] = tmp_adata  # assign sample_id to barcodes

for k in cnvan_key_l[1:]:
    assert np.all(k == cnvan_key_l[0])

cnvan_key = cnvan_key_l[-1]
raw = ad.concat(adatas_d, index_unique="_")

# Use conversion key to re-assign symbols to ensembl ids
raw.var.loc[cnvan_key.index, "gene_ids"] = cnvan_key.gene_ids

raw.X = scipy.sparse.csr_matrix(raw.X)

raw.obs["group"] = pd.Categorical(
    raw.obs["group"], categories=["Ctrl", "A8", "A9", "A8A9"]
)
raw.obs["sample"] = pd.Categorical(
    raw.obs["sample"], categories=raw.obs["sample"].unique()
)
raw.obs["patient"] = pd.Categorical(raw.obs["patient"])
raw.obs["batch"] = pd.Categorical(raw.obs["batch"])

# Save adata
raw.write(f"{resDir}/raw.h5ad", compression="gzip")
