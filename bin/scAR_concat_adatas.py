#!/usr/bin/env python3

import argparse
import sys

import pandas as pd


def parse_args(args=None):
    Description = "Concatenates individually denoised scAR output samples by study"
    Epilog = "Example usage: python scAR_concat_adatas.py <adatas_path> <samplesheet>  <var_names>  <gtf_file> <hgnc_file> <cpus>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("--adatas_path")
    parser.add_argument("--samplesheet")
    parser.add_argument("--var_names")
    parser.add_argument("--gtf_file")
    parser.add_argument("--hgnc_file")
    parser.add_argument("--cpus", type=int)
    return parser.parse_args(args)


def get_ordered_h5ad_file_paths(samplesheet):
    from pathlib import Path

    # Get file paths and reorder by study and sample id
    files = [str(x) for x in Path(".").glob("*.h5ad")]
    names = [x.stem.replace("_denoised_adata", "") for x in Path(".").glob("*.h5ad")]

    files = pd.merge(
        samplesheet,
        pd.DataFrame({"files": files, "sample_id": names}),
        how="left",
        on=["sample_id"],
        validate="m:1",
    )["files"].tolist()

    # Remove float("nan") that are generated if input is only a subset of samples in samplesheet, sc.read_h5ad() will raise TypeError if given nan!
    files = [x for x in files if not isinstance(x, float)]
    return files


def get_gtf_annotation(gtf_file, hgnc_file, var_names):
    var_names = pd.read_csv(var_names, names=["var_names"], index_col=0)

    gtf = pd.read_csv(gtf_file)
    gtf["ensembl"] = gtf["Geneid"].str.split(".").str[0]
    # Append gene info from gtf file and keep original adata var_names as index (should be equivalent to "GeneSymbol" column in gtf file + suffix if not unique)
    gtf["var_names"] = var_names.index
    # Remove sex chromosome ids if ensembl id is duplicated
    gtf = gtf[~(gtf["ensembl"].duplicated() & gtf["Geneid"].str.contains("_PAR_Y"))]

    hgnc = pd.read_csv(hgnc_file, delimiter="\t")
    hgnc = hgnc[
        [
            "ensembl_gene_id",
            "entrez_id",
            "hgnc_id",
            "uniprot_ids",
            "alias_symbol",
            "symbol",
        ]
    ]
    hgnc = hgnc.dropna(subset=["ensembl_gene_id"])
    # Drop single duplicated ensembl_gene_id, keep symbol from gtf
    hgnc = hgnc[hgnc["symbol"] != "LINC00856"]
    hgnc = hgnc[
        ["ensembl_gene_id", "entrez_id", "hgnc_id", "uniprot_ids", "alias_symbol"]
    ]

    gtf = pd.merge(
        gtf,
        hgnc,
        how="left",
        left_on="ensembl",
        right_on="ensembl_gene_id",
        validate="m:1",
    )

    return gtf


def concat_denoised_adatas(files, gtf, cpus):
    import anndata
    import scanpy as sc
    from tqdm.contrib.concurrent import process_map

    # Concat adatas from specified directory
    adatas = process_map(sc.read_h5ad, files, max_workers=cpus)
    adata = anndata.concat(adatas, join="outer", index_unique="_")
    assert adata.obs_names.is_unique

    adata.var = pd.merge(
        pd.DataFrame({"var_names": adata.var_names}),
        gtf,
        how="left",
        on=["var_names"],
        validate="m:1",
    ).set_index("var_names")
    assert adata.var_names.is_unique
    adata.var_names.name = None

    for var in ["n_counts", "n_genes", "_scvi_batch", "_scvi_labels"]:
        del adata.obs[var]

    adata.write_h5ad("denoised_adata.h5ad", compression="lzf")


def main(args=None):
    args = parse_args(args)
    files = get_ordered_h5ad_file_paths(samplesheet=pd.read_csv(args.samplesheet))
    gtf = get_gtf_annotation(
        gtf_file=args.gtf_file, hgnc_file=args.hgnc_file, var_names=args.var_names
    )
    concat_denoised_adatas(files, gtf, cpus=args.cpus)


if __name__ == "__main__":
    sys.exit(main())