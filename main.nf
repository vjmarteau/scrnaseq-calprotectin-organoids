#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { LOAD_ADATA } from "./modules/LOAD_ADATA"
include { LOAD_RAW } from "./modules/LOAD_RAW"
include { FILTER_MITO_OUTLIERS } from "./modules/FILTER_MITO_OUTLIERS"
include { H5AD_TO_SEURAT } from "./modules/H5AD_TO_SEURAT"

workflow {
    // Retrieve and validate parameters
    assert params.samplesheet != null : "Please specify the `samplesheet` parameter"
    samplesheet = file(params.samplesheet, checkIfExists: true)
    ch_input_files = Channel.fromPath(params.input_path)

    // start workflow
    LOAD_ADATA(samplesheet, ch_input_files)
    LOAD_RAW(LOAD_ADATA.out.meta, ch_input_files)
    FILTER_MITO_OUTLIERS(LOAD_ADATA.out.adata)
    H5AD_TO_SEURAT(LOAD_ADATA.out.adata)
    }