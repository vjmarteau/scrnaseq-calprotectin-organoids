#!/usr/bin/env python3

import argparse
import pathlib
import sys

import pandas as pd
import scanpy as sc


def parse_args(args=None):
    Description = "Load 10X raw and filtered sample counts and perform basic filtering. Output can be passed to scAR denoising algorithm. Assumes unique sample ids"
    Epilog = "Example usage: python check_samplesheet.py <cellranger_path> <samplesheet> <cpus>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("--cellranger_path", type=pathlib.Path)
    parser.add_argument("--samplesheet")
    parser.add_argument("--cpus", type=int)
    return parser.parse_args(args)


def tidy_samplesheet(samplesheet):
    samplesheet = pd.read_csv(samplesheet)
    samplesheet.rename(
        columns={"sample": "sample_id", "patient": "patient_id"}, inplace=True
    )
    samplesheet.drop_duplicates(subset=["sample_id"], inplace=True)

    # Sort by batch and sample_id
    samplesheet["tmp_int"] = (
        samplesheet["sample_id"].str.extract("(\d+)", expand=False).astype(int)
    )
    samplesheet = samplesheet.sort_values(["batch", "tmp_int"])
    samplesheet = samplesheet.drop("tmp_int", axis=1)
    samplesheet = samplesheet.reset_index(drop=True)
    samplesheet.to_csv(f"samplesheet.csv", index=False)
    return samplesheet


def load_counts(sample_meta, cellranger_path, mat_type="filtered"):
    import scipy.sparse

    file_path = f"{cellranger_path}/{sample_meta['sample_id']}/outs/{mat_type}_feature_bc_matrix.h5"
    adata = sc.read_10x_h5(file_path)
    adata.obs = adata.obs.assign(**sample_meta)
    adata.var_names_make_unique()
    adata.X = scipy.sparse.csr_matrix(adata.X)
    adata.obs_names = (
        adata.obs["sample_id"] + "_" + adata.obs_names.str.replace("-1", "")
    )
    assert adata.obs_names.is_unique
    return adata


def save_adatas_by_sample_and_mat_type(samplesheet, cellranger_path, cpus, load_counts=load_counts):
    import anndata
    import itertools
    from tqdm.contrib.concurrent import process_map


    for mat_type in ["filtered", "raw"]:
        adatas = process_map(
            load_counts,
            [r for _, r in samplesheet.iterrows()],
            itertools.repeat(cellranger_path),
            itertools.repeat(mat_type),
            max_workers=cpus,
        )
        if mat_type == "filtered":
            adata = anndata.concat(adatas, join="outer")
            pd.DataFrame(adata.var_names).to_csv(
                "var_names.csv", index=False, header=False
            )

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
            adata.write_h5ad(f"{sample_id}_{mat_type}_adata.h5ad", compression="lzf")


def main(args=None):
    args = parse_args(args)
    samplesheet = tidy_samplesheet(args.samplesheet)
    
    save_adatas_by_sample_and_mat_type(
        samplesheet=samplesheet,
        cellranger_path=args.cellranger_path,
        cpus=args.cpus,
    )


if __name__ == "__main__":
    sys.exit(main())