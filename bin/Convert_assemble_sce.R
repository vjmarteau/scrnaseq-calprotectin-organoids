#!/usr/bin/env Rscript
'
Usage:
  Convert_assemble_sce.R --counts_mtx=<counts_mtx> --denoised_mtx=<denoised_mtx> --features=<features> --barcodes=<barcodes> --metadata=<metadata> --metadata_var=<metadata_var> --umap_tsv=<umap_tsv>[options]

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
library(tibble)
library(Matrix)
library(SingleCellExperiment)

# Load parameters
output_file <- arguments$output_file

load_mat <- function(mat_path, features_path, barcodes_path) {

  mat <- readMM(mat_path)
  genes <- read_tsv(features_path, col_names = FALSE)
  gene_ids <- genes$X2
  cell_ids <- read_tsv(barcodes_path, col_names = FALSE)$X1
  rownames(mat) <- gene_ids
  colnames(mat) <- cell_ids

  return(mat)
}

counts <- load_mat(
  mat_path = arguments$counts_mtx,
  features_path = arguments$features,
  barcodes_path = arguments$barcodes
)

denoised <- load_mat(
  mat_path = arguments$denoised_mtx,
  features_path = arguments$features,
  barcodes_path = arguments$barcodes
)

metadata <- read_tsv(arguments$metadata) |>
  as_data_frame() |>
  column_to_rownames(var = "...1")

metadata_var <- read_tsv(arguments$metadata_var) |>
  as_data_frame() |>
  column_to_rownames(var = "...1")

# Add embedding
embedding <- read_tsv(arguments$umap_tsv) |> as.matrix()
rownames(embedding) <- colnames(counts)
colnames(embedding) <- c("umap_1", "umap_2")

#sce <- SingleCellExperiment(
#  assays      = list(counts = counts,
#                     denoised = denoised),
#  colData     = metadata,
#  rowData     = metadata_var,
#  reducedDims = list(UMAP = embedding)
#)

counts_sce <- SingleCellExperiment(
  assays      = list(counts = counts),
  colData     = metadata,
  rowData     = metadata_var,
  reducedDims = list(UMAP = embedding)
)

saveRDS(counts_sce, file = file.path(paste0(output_file, "counts_sce.rds")))

denoised_sce <- SingleCellExperiment(
  assays      = list(counts = denoised),
  colData     = metadata,
  rowData     = metadata_var,
  reducedDims = list(UMAP = embedding)
)

saveRDS(denoised_sce, file = file.path(paste0(output_file, "denoised_sce.rds")))