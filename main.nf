#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { Remove_ambient_RNA_scAR } from "${baseDir}/workflows/Remove_ambient_RNA_scAR"
include { PRE_PROCESS } from "./workflows/Pre-process"

workflow {
    // Retrieve and validate parameters
    assert params.samplesheet != null : "Please specify the `samplesheet` parameter"
    samplesheet = file(params.samplesheet, checkIfExists: true)
    ch_input_files = Channel.fromPath(params.cellranger)
    gtf_file = file(params.gencode_gtf, checkIfExists: true)
    hgnc_file = file(params.hgnc, checkIfExists: true)
    //cell_cycle_genes = file(params.cell_cycle_genes, checkIfExists: true)
    
    // start workflow
    Remove_ambient_RNA_scAR(samplesheet, ch_input_files, gtf_file, hgnc_file)
    // PRE_PROCESS(samplesheet, ch_input_files, cell_cycle_genes)
    }