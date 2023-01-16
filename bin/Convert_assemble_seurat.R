#!/usr/bin/env Rscript
'
Usage:
  Convert_assemble_seurat.R --counts_mtx=<counts_mtx> --denoised_mtx=<denoised_mtx> --features=<features> --barcodes=<barcodes> --metadata=<metadata> --metadata_var=<metadata_var> --umap_tsv=<umap_tsv>[options]

Mandatory arguments:
  --counts_mtx=<counts_mtx>   Count matrix
  --denoised_mtx=<denoised_mtx>     Experiment metadata
  --features=<features>               Genes of interest in txt format
  --barcodes=<barcodes>
  --metadata=<metadata>
  --metadata_var=<metadata_var>
  --umap_tsv=<umap_tsv>

  --output_file=<output_file>   Path to output file
' -> doc

library(conflicted)
library(docopt)
arguments <- docopt(doc, version = "0.1")
print(arguments)

library(readr)
library(dplyr)
library(forcats)
library(tibble)
library(Seurat)

# Load parameters
output_file <- arguments$output_file

# Get the expression matrix
exprs <- ReadMtx(
  mtx = arguments$counts_mtx,
  features = arguments$features,
  cells = arguments$barcodes
)

# Create the Seurat object
seurat <- CreateSeuratObject(exprs)

# Set the expression assay
seurat <- SetAssayData(seurat, "counts", exprs)

# Add observation metadata
metadata <- read_tsv(arguments$metadata) |> as_data_frame() |> column_to_rownames(var = "...1")
seurat <- AddMetaData(seurat, metadata)

# Add fetaure metadata
metadata_var <- read_tsv(arguments$metadata_var) |> as_data_frame() |> column_to_rownames(var = "...1")
seurat[["RNA"]][["gene_ids"]] <- metadata_var$gene_ids
seurat[["RNA"]][["n_cells"]] <- metadata_var$n_cells

# Get denoised expression matrix
#exprs_d <- ReadMtx(
#  mtx = arguments$denoised_mtx,
#  features = arguments$features,
#  cells = arguments$barcodes
#)

# Add denoised expression assay
#seurat[["denoised"]] <- SetAssayData(seurat, "data", exprs_d)

# Add embedding
embedding <- read_tsv(arguments$umap_tsv) |> as.matrix()
rownames(embedding) <- colnames(x = seurat)
colnames(embedding) <- c("umap_1", "umap_2")
seurat[["umap"]] <- CreateDimReducObject(embedding, key = "umap_")

seurat <- SetIdent(seurat, value = "cell_type")

saveRDS(seurat, file = file.path(paste0(output_file, "seurat_annotated.rds")))
