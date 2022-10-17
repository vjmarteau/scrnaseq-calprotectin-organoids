#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { Load_adata } from "./modules/Load_adata"

workflow {
    // Retrieve and validate parameters
    assert params.samplesheet != null : "Please specify the `samplesheet` parameter"
    samplesheet = file(params.samplesheet, checkIfExists: true)

    // start workflow
    Load_adata(samplesheet)
}
