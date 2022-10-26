#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { Load_adata } from "./modules/Load_adata"

workflow {
    // Retrieve and validate parameters
    assert params.samplesheet != null : "Please specify the `samplesheet` parameter"
    samplesheet = file(params.samplesheet, checkIfExists: true)
    ch_input_files = Channel.fromPath(params.input_path)

    // start workflow
    Load_adata(samplesheet, ch_input_files)
}