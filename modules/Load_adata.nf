nextflow.enable.dsl=2

out_dir = file(params.resDir)
mode = params.publish_dir_mode

process Load_adata {
    publishDir "${out_dir}", mode: "$mode"

    input:
        path(samplesheet)
        path(ch_input_files)

    output:
        path("metadata.csv"), emit: metadata
        path("adata.h5ad"), emit: adata

	script:
	"""
    01-Load_adata.py \\
    --samplesheet=${samplesheet}
	"""
}
