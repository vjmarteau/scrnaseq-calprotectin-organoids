nextflow.enable.dsl=2

out_dir = file(params.outdir)
mode = params.publish_dir_mode

process LOAD_ADATA {
    publishDir "${out_dir}", mode: "$mode"

    input:
        path(samplesheet)
        path(ch_input_files)

    output:
        path("metadata.csv"), emit: meta
        path("adata.h5ad"), emit: adata

	script:
	"""
    QC-Load_adata.py \\
    --samplesheet=${samplesheet}
	"""
}

process LOAD_RAW {
    publishDir "${out_dir}", mode: "$mode"

    input:
        path(samplesheet)
        path(ch_input_files)

    output:
        path("raw.h5ad"), emit: raw_adata

	script:
	"""
    QC-Load_raw.py \\
    --samplesheet=${samplesheet}
	"""
}

process RUN_SCAR {
    publishDir "${out_dir}", mode: "$mode"
    label "gpu"

    input:
        path(raw_adata)
        path(adata)
        path(cell_cycle_genes)

    output:
        path("denoised_adata.h5ad"), emit: denoised_adata

	script:
	"""
    QC-Run_scar.py \\
    --raw_adata=${raw_adata} \\
    --adata=${adata} \\
    --cell_cycle_genes=${cell_cycle_genes}
	"""
}

process RUN_SCVI_AND_SOLO {
    publishDir "${out_dir}", mode: "$mode"
    label "gpu"

    input:
        path(adata)

    output:
        path("scVI_model"), emit: scVI_model
        path("adata_nodoublet.h5ad"), emit: adata_nodoublet
        path("is_doublet.png")

	script:
	"""
    QC-Run_scVI_and_SOLO.py \\
    --adata=${adata}
	"""
}

workflow PRE_PROCESS {
    take:
        samplesheet
        ch_input_files
        cell_cycle_genes

    main:
        LOAD_ADATA(samplesheet, ch_input_files)
        LOAD_RAW(samplesheet, ch_input_files)
        RUN_SCAR(LOAD_RAW.out.raw_adata, LOAD_ADATA.out.adata, cell_cycle_genes)
        RUN_SCVI_AND_SOLO(RUN_SCAR.out.denoised_adata)

    emit:
        raw_adata = LOAD_RAW.out.raw_adata
        adata = LOAD_ADATA.out.adata
        denoised_adata = RUN_SCAR.out.denoised_adata
        adata_nodoublet = RUN_SCVI_AND_SOLO.out.adata_nodoublet
}