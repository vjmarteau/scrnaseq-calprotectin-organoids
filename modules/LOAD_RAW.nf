nextflow.enable.dsl=2

out_dir = file(params.resDir)
mode = params.publish_dir_mode

process LOAD_RAW {
    publishDir "${out_dir}", mode: "$mode"

    input:
        path(meta)
        path(ch_input_files)

    output:
        path("raw.h5ad"), emit: raw

	script:
	"""
    Load_raw.py \\
    --meta=${meta}
	"""
}
