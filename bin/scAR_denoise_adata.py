#!/usr/bin/env python3

"""
Usage:
  scAR_denoise_adata.py --raw_adata=<raw_adata> --filtered_adata=<filtered_adata> --cpus=${task.cpus} [options]

Mandatory arguments:
  --raw_adata=<raw_adata>             raw_adata
  --filtered_adata=<filtered_adata>   filtered_adata
  --cpus=${task.cpus}                 Number of cpus
  --output_file=<output_file>
"""

from docopt import docopt
import scvi
import anndata
import scanpy as sc

from threadpoolctl import threadpool_limits
import multiprocessing

def set_all_seeds(seed=0):
    import os
    import random
    import numpy as np
    import torch
    scvi.settings.seed = seed
    os.environ["PYTHONHASHSEED"] = str(seed)  # Python general
    os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"
    np.random.seed(seed)  # Numpy random
    random.seed(seed)  # Python random
    torch.manual_seed(seed)
    torch.use_deterministic_algorithms(True)
    if torch.cuda.is_available():
        torch.cuda.manual_seed(seed)
        torch.cuda.manual_seed_all(seed)  # For multiGPU

args = docopt(__doc__)
raw_adata = args["--raw_adata"]
filtered_adata = args["--filtered_adata"]
cpus = int(args["--cpus"])
output_file = args["--output_file"]

threadpool_limits(cpus)
set_all_seeds()

raw_adata = sc.read_h5ad(raw_adata)
filtered_adata = sc.read_h5ad(filtered_adata)


# Ambient RNA removal (scAR)

scvi.external.SCAR.setup_anndata(filtered_adata)
scvi.external.SCAR.get_ambient_profile(
adata=filtered_adata, raw_adata=raw_adata, prob=0.995
)

model = scvi.external.SCAR(filtered_adata)
model.train(early_stopping=True, use_gpu=True)

filtered_adata.obsm["X_scAR"] = model.get_latent_representation()
filtered_adata.layers["denoised"] = model.get_denoised_counts()

# Save raw counts to layer counts
filtered_adata.layers["counts"] = filtered_adata.X.copy()

filtered_adata.write_h5ad(output_file)