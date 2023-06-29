#!/usr/bin/env python3

"""
Usage:
  scAR_concat_adatas.py --adatas_path=<adatas_path> --metadata=<metadata> --var_names=<var_names> --gtf_file=<gtf_file> --hgnc_file=<hgnc_file> --cpus=<cpus> [options]

Mandatory arguments:
  --adatas_path=<adatas_path>   Path to adatas that should be concatenated
  --metadata=<metadata>
  --var_names=<var_names>
  --gtf_file=<gtf_file>         Gtf file that was used to annotate genes in nf-core/scrnaseq pipeline
  -hgnc_file=<hgnc_file>
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
hgnc_file = args["--hgnc_file"]
cpus = int(args["--cpus"])


# Get metadata
meta = pd.read_csv(metadata)

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
var_names = pd.read_csv(var_names, names=['var_names'], index_col=0)

# Get gene annotations
gtf = pd.read_csv(gtf_file)
gtf["ensembl"] = gtf["Geneid"].str.split(".").str[0]
gtf["var_names"] = var_names.index
# Remove sex chromosome ids if ensembl id is duplicated
gtf = gtf[~(gtf["ensembl"].duplicated() & gtf["Geneid"].str.contains("_PAR_Y"))]

# Get hgnc annotations
hgnc = pd.read_csv(hgnc_file, delimiter="\t")
hgnc = hgnc[["ensembl_gene_id", "entrez_id", "hgnc_id", "uniprot_ids", "alias_symbol", "symbol"]]
hgnc = hgnc.dropna(subset=['ensembl_gene_id'])
# Drop single duplicated ensembl_gene_id, keep symbol from gtf
hgnc = hgnc[hgnc['symbol'] != 'LINC00856']
hgnc = hgnc[["ensembl_gene_id", "entrez_id", "hgnc_id", "uniprot_ids", "alias_symbol"]]

gtf = pd.merge(
    gtf,
    hgnc,
    how="left",
    left_on="ensembl",
    right_on="ensembl_gene_id",
    validate="m:1",
)

# Concat adatas from specified directory
adatas = process_map(sc.read_h5ad, files, max_workers=cpus)
adata = anndata.concat(adatas, join="outer")
assert adata.obs_names.is_unique

# Append gene info from gtf file and keep original adata var_names as index (should be equivalent to "GeneSymbol" column in gtf file + suffix if not unique)
adata.var = pd.merge(
    pd.DataFrame({"var_names": adata.var_names}),
    gtf,
    how="left",
    on=["var_names"],
    validate="m:1",
).set_index("var_names")

assert adata.var_names.is_unique
adata.var_names.name = None

# Filter genes on denoised counts
adata.X = adata.layers["denoised"]
sc.pp.filter_genes(adata, min_cells=10)
sc.pp.filter_genes(adata, min_counts=10)
adata.X = adata.layers["counts"]

adata.write_h5ad("adata.h5ad")