nextflow.enable.dsl=2

out_dir = file(params.resDir)
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

    output:
        path("denoised_adata.h5ad"), emit: denoised_adata

	script:
	"""
    QC-Run_scar.py \\
    --raw_adata=${raw_adata} \\
    --adata=${adata}    
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

    main:
        LOAD_ADATA(samplesheet, ch_input_files)
        LOAD_RAW(samplesheet, ch_input_files)
        RUN_SCAR(LOAD_RAW.out.raw_adata, LOAD_ADATA.out.adata)
        RUN_SCVI_AND_SOLO(RUN_SCAR.out.denoised_adata)

    emit:
        raw_adata = LOAD_RAW.out.raw_adata
        adata = LOAD_ADATA.out.adata
        denoised_adata = RUN_SCAR.out.denoised_adata
        adata_nodoublet = RUN_SCVI_AND_SOLO.out.adata_nodoublet
}