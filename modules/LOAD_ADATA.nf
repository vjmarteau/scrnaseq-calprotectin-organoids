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
    Load_adata.py \\
    --samplesheet=${samplesheet}
	"""
}
