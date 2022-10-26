#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { Load_adata } from "./modules/Load_adata"
include { H5AD_TO_SEURAT } from "./modules/H5AD_TO_SEURAT"

workflow {
    // Retrieve and validate parameters
    assert params.samplesheet != null : "Please specify the `samplesheet` parameter"
    samplesheet = file(params.samplesheet, checkIfExists: true)
    ch_input_files = Channel.fromPath(params.input_path)

    // start workflow
    Load_adata(samplesheet, ch_input_files)
    H5AD_TO_SEURAT(Load_adata.out.adata)
}