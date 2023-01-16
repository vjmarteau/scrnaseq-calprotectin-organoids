#!/usr/bin/env python3

"""
Usage:
  Convert_extract_adata.py --adata=<adata> [options]

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
import scipy.io
import os
import gzip

args = docopt(__doc__)
adata = args["--adata"]
resDir = args["--resDir"]

adata = sc.read_h5ad(adata)

sc.tl.leiden(adata, key_added="leiden", resolution=0.5)

# convert to CSC if possible. See https://github.com/MarioniLab/scran/issues/70
def convert_to_CSC(data_mat):
    if scipy.sparse.issparse(data_mat):
        if data_mat.nnz > 2**31 - 1:
            data_mat = data_mat.tocoo()
        else:
            data_mat = data_mat.tocsc()
    
    return(data_mat)
  
scipy.io.mmwrite(
    os.path.join(resDir, "counts_matrix.mtx"),
    convert_to_CSC(
        scipy.sparse.csr_matrix(adata.layers["counts"].T),
    ),
)

scipy.io.mmwrite(
    os.path.join(resDir, "denoised_matrix.mtx"),
    convert_to_CSC(
        scipy.sparse.csr_matrix(adata.layers["denoised"].T),
    ),
)

features = pd.DataFrame(adata.var.gene_ids).reset_index()
features = features[["gene_ids", "index"]].rename(
    columns={"gene_ids": "0", "index": "1"}
)
features.to_csv(
    os.path.join(resDir, "features.tsv"),
    sep="\t",
    index=False,
    header=False,
    compression="gzip",
)

pd.DataFrame(adata.obs.index).to_csv(
    os.path.join(resDir, "barcodes.tsv"),
    sep="\t",
    index=False,
    header=False,
    compression="gzip",
)

adata.obs.to_csv(os.path.join(resDir, "metadata.tsv"), sep="\t", index=True)
adata.var.to_csv(os.path.join(resDir, "metadata_var.tsv"), sep="\t", index=True)

pd.DataFrame(adata.obsm["X_umap"]).to_csv(
    os.path.join(resDir, "umap.tsv"), sep="\t", index=False, compression="gzip"
)