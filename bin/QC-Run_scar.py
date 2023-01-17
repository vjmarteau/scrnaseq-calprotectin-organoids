#!/usr/bin/env python3

"""
Usage:
  Run_scar.py --raw_adata=<raw_adata> --adata=<adata> --cell_cycle_genes=<cell_cycle_genes> [options]

Mandatory arguments:
  --raw_adata=<raw_adata>   raw_adata
  --adata=<adata>           adata
  --cell_cycle_genes=<cell_cycle_genes>

Optional arguments:
  --resDir=<resDir>     Output directory [default: ./]
"""

from docopt import docopt
import scvi
import anndata as ad
import scanpy as sc
import pandas as pd

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
adata = args["--adata"]
cell_cycle_genes = args["--cell_cycle_genes"]
resDir = args["--resDir"]

threadpool_limits(16)
set_all_seeds()

raw_adata = sc.read_h5ad(raw_adata)
adata = sc.read_h5ad(adata)


# Ambient RNA removal (scAR)
sample_d = dict()
for s in adata.obs['sample'].unique():

    raw_adata_s = raw_adata[raw_adata.obs["sample"] == s, :].copy()
    raw_adata_s.obs["value"] = 0

    adata_s = adata[adata.obs["sample"] == s, :].copy()
    adata_s.obs["value"] = 0

    scvi.external.SCAR.setup_anndata(adata_s)
    scvi.external.SCAR.get_ambient_profile(adata=adata_s, raw_adata=raw_adata_s, prob=0.995)

    model = scvi.external.SCAR(adata_s)
    model.train(early_stopping=True, use_gpu=True)

    adata_s.obsm["X_scAR"] = model.get_latent_representation()
    adata_s.layers['denoised'] = model.get_denoised_counts()

    sample_d[s] = adata_s

# Integrate back
adata = ad.concat(sample_d, label='sample' , merge="unique")

#adata.write(f"{resDir}/denoised_adata.h5ad", compression="gzip")

# Save raw counts to layer counts and use denoised counts for subsequent analysis
adata.layers['counts'] = adata.X.copy()
adata.X = adata.layers['denoised']

# very basic cell/gene filtering on denoised counts
sc.pp.filter_cells(adata, min_counts=800)
sc.pp.filter_cells(adata, min_genes=200)
sc.pp.filter_cells(adata, max_counts=100000)

#sc.pp.filter_genes(adata, min_cells=10)
#sc.pp.filter_genes(adata, min_counts=10)

# annotate the group of mitochondrial genes as 'mito'
adata.var['mito'] = adata.var_names.str.startswith('MT-')
# ribosomal genes as ribo
adata.var['ribo'] = adata.var_names.str.startswith(('RPS', 'RPL'))

sc.pp.calculate_qc_metrics(adata, qc_vars=['mito', 'ribo'], percent_top=None, log1p=False, inplace=True)

adata = adata[adata.obs["pct_counts_mito"] < 25].copy()


# Annotate cell cycle phase
cell_cycle_genes = [
    x.strip() for x in open(cell_cycle_genes)
]

# Split into 2 lists
s_genes = cell_cycle_genes[:43]
g2m_genes = cell_cycle_genes[43:]
cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]

adata_cycle = adata.copy()
sc.pp.normalize_total(adata_cycle)
sc.pp.log1p(adata_cycle)
sc.pp.scale(adata_cycle)

sc.tl.score_genes_cell_cycle(adata_cycle, s_genes=s_genes, g2m_genes=g2m_genes)

adata.obs['S_score'] = adata_cycle.obs['S_score']
adata.obs['G2M_score'] = adata_cycle.obs['G2M_score']
adata.obs['phase'] = adata_cycle.obs['phase']

adata.write(f"{resDir}/denoised_adata.h5ad", compression="gzip")