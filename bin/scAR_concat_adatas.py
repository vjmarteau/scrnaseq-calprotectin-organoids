#!/usr/bin/env python3

"""
Usage:
  scAR_concat_adatas.py --adatas_path=<adatas_path> --metadata=<metadata> --var_names=<var_names> --gtf_file=<gtf_file> --cpus=<cpus> [options]

Mandatory arguments:
  --adatas_path=<adatas_path>   Path to adatas that should be concatenated
  --metadata=<metadata>
  --var_names=<var_names>
  --gtf_file=<gtf_file>         Gtf file that was used to annotate genes in nf-core/scrnaseq pipeline
  --cpus=<cpus>                 Number of cpus
"""
from docopt import docopt
import anndata
import scanpy as sc
import pandas as pd
from pathlib import Path
from tqdm.contrib.concurrent import process_map

args = docopt(__doc__)
metadata = args["--metadata"]
var_names = args["--var_names"]
gtf_file = args["--gtf_file"]
cpus = int(args["--cpus"])

# Get metadata
meta = pd.read_csv(metadata)

meta["prefix"] = meta["sample_id"].str.split(r'(\d+)', expand=True)[0]
meta["no"] = meta["sample_id"].str.split(r'(\d+)', expand=True)[1].astype(float)

meta = meta.sort_values(by = ["batch", "prefix", "no"])
meta = meta.drop(columns=["prefix", "no"]).reset_index(drop=True)

# Get file paths and reorder by sample id
files = [str(x) for x in Path(".").glob("*.h5ad")]
names = [x.stem.replace("_denoised_adata", "") for x in Path(".").glob("*.h5ad")]

sample_idx = pd.DataFrame({"files": files, "sample_id": names})

files = pd.merge(
    meta,
    sample_idx,
    how="left",
    on=["sample_id"],
    validate="m:1",
)["files"].tolist()

# Get original var_names with all genes
var_names = pd.read_csv(var_names, delimiter="\t", names=['var_names'], index_col=0)

# Get gene annotations
gtf = pd.read_csv(gtf_file, delimiter="\t", skipinitialspace = True)
# Remove Ensembl version number
gtf.insert(0, "Ensembl", gtf["Geneid"].str.replace(r"\.[^.]+$", ""))
# Append back sex chromosome info to new column "Ensembl"
condi = gtf["Geneid"].str.contains("_PAR_Y")
gtf.loc[condi, "Ensembl"] = gtf.loc[condi, "Ensembl"] + "_PAR_Y"

# Concat adatas from specified directory
adatas = process_map(sc.read_h5ad, files, max_workers=cpus)
adata = anndata.concat(adatas, join="outer")
assert adata.obs_names.is_unique

# Append gene info from gtf file and keep original adata var_names as index (should be equivalent to "GeneSymbol" column in gtf file + suffix if not unique)
#gtf.index = adata.var_names
#adata.var = gtf
#assert adata.var_names.is_unique

gtf["var_names"] = var_names.index

adata.var = pd.merge(
    pd.DataFrame({"var_names": adata.var_names}),
    gtf,
    how="left",
    on=["var_names"],
    validate="m:1",
).set_index("var_names")

assert adata.var_names.is_unique
adata.var_names.name = None

adata.write_h5ad("adata.h5ad")