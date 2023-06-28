#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { PRE_PROCESS } from "./workflows/Pre-process"

workflow {
    // Retrieve and validate parameters
    assert params.samplesheet != null : "Please specify the `samplesheet` parameter"
    samplesheet = file(params.samplesheet, checkIfExists: true)
    ch_input_files = Channel.fromPath(params.input_path)
    cell_cycle_genes = file(params.cell_cycle_genes, checkIfExists: true)

    // start workflow
    PRE_PROCESS(samplesheet, ch_input_files, cell_cycle_genes)
    }