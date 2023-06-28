#!/usr/bin/env python3

"""
Usage:
  scAR_load_adata.py --samplesheet=<samplesheet> --cpus=<cpus> [options]

Mandatory arguments:
  --samplesheet=<samplesheet>   Samplesheet from nfcore/scrnaseq pipeline
  --cpus=<cpus>                 Number of cpus
"""
from docopt import docopt
import anndata
import pandas as pd
import scanpy as sc
import scipy.sparse
from tqdm.contrib.concurrent import process_map
import itertools

args = docopt(__doc__)
samplesheet = args["--samplesheet"]
cpus = int(args["--cpus"])

# Get metadata from samplesheet
meta = pd.read_csv(samplesheet)
meta.rename(columns={"sample": "sample_id", "patient": "patient_id"}, inplace=True)
meta.drop_duplicates(subset=['sample_id'], inplace=True)

# Sort by id
meta['tmp_int'] = meta['sample_id'].str.extract('(\d+)', expand=False).astype(int)
meta = meta.sort_values(['batch', 'tmp_int'])
meta = meta.drop('tmp_int', axis=1)
meta = meta.reset_index(drop=True)

# Save metadata as csv
meta.to_csv(f"metadata.csv", index=False)


def load_counts(sample_meta, mat_type="filtered"):
    file_path = f"./30_nfcore_scrnaseq/cellranger/{sample_meta['sample_id']}/outs/{mat_type}_feature_bc_matrix.h5"
    adata = sc.read_10x_h5(file_path)
    adata.obs = adata.obs.assign(**sample_meta)
    adata.var_names_make_unique()
    adata.X = scipy.sparse.csr_matrix(adata.X)
    adata.obs_names = (
        adata.obs["sample_id"] + "_" + adata.obs_names.str.replace("-1", "")
    )
    assert adata.obs_names.is_unique
    return adata

for mat_type in ["filtered", "raw"]:
    adatas = process_map(
        load_counts,
        [r for _, r in meta.iterrows()],
        itertools.repeat(mat_type),
        max_workers=cpus,
    )
    if mat_type == "filtered":
        adata = anndata.concat(adatas, join="outer")
        pd.DataFrame(adata.var_names).to_csv("var_names.csv", index=False, header=False)
        
        sc.pp.filter_genes(adata, min_cells=10)
        sc.pp.filter_genes(adata, min_counts=10)
        sc.pp.filter_cells(adata, min_counts=500)
        sc.pp.filter_cells(adata, min_genes=200)

        adatas = [
            adata[adata.obs["sample_id"] == sample, :].copy()
            for sample in adata.obs["sample_id"].unique()
        ]

    for adata in adatas:
        sample_id = adata.obs["sample_id"].unique().tolist()[0]
        adata.write_h5ad(f"{sample_id}_{mat_type}_adata.h5ad")
