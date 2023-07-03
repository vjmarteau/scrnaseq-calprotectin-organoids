#!/usr/bin/env python3

import argparse
import sys

import scanpy as sc
import scvi


def parse_args(args=None):
    Description = "Run Ambient RNA removal (scAR) algorithm per sample"
    Epilog = "Example usage: python scAR_denoise_adatas.py <raw_adata> <filtered_adata> <cpus>  <output_file>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("--raw_adata")
    parser.add_argument("--filtered_adata")
    parser.add_argument("--cpus", type=int)
    parser.add_argument("--output_file")
    return parser.parse_args(args)


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


def Run_scAR(raw_adata, filtered_adata, cpus, output_file, set_all_seeds=set_all_seeds):
    import multiprocessing
    from threadpoolctl import threadpool_limits

    threadpool_limits(cpus)
    set_all_seeds()

    # Ambient RNA removal (scAR)
    scvi.external.SCAR.setup_anndata(filtered_adata)
    scvi.external.SCAR.get_ambient_profile(
        adata=filtered_adata, raw_adata=raw_adata, prob=0.995
    )
    model = scvi.external.SCAR(filtered_adata)
    model.train(early_stopping=True)

    filtered_adata.obsm["X_scAR"] = model.get_latent_representation()
    filtered_adata.layers["denoised"] = model.get_denoised_counts()
    filtered_adata.layers["counts"] = filtered_adata.X.copy()

    filtered_adata.write_h5ad(output_file, compression="lzf")


def main(args=None):
    args = parse_args(args)

    Run_scAR(
        raw_adata=sc.read_h5ad(args.raw_adata),
        filtered_adata=sc.read_h5ad(args.filtered_adata),
        cpus=args.cpus,
        output_file=args.output_file,
    )


if __name__ == "__main__":
    sys.exit(main())