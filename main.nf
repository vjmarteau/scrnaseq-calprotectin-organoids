#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { LOAD_ADATA } from "./modules/LOAD_ADATA"

workflow {
    // Retrieve and validate parameters
    assert params.samplesheet != null : "Please specify the `samplesheet` parameter"
    samplesheet = file(params.samplesheet, checkIfExists: true)
    ch_input_files = Channel.fromPath(params.input_path)

    // start workflow
    LOAD_ADATA(samplesheet, ch_input_files)
}