nextflow.enable.dsl=2

out_dir = file(params.resDir)
mode = params.publish_dir_mode

process Load_adata {
    publishDir "${out_dir}", mode: "$mode"

    input:
        path(samplesheet)

    output:
        path("metadata.csv"), emit: metadata

	script:
	"""
    Load_adata.py \\
    --samplesheet=${samplesheet}
	"""
}
