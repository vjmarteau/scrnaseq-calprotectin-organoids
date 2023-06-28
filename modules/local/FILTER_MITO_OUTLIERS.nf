nextflow.enable.dsl=2

out_dir = file(params.resDir)
mode = params.publish_dir_mode

process FILTER_MITO_OUTLIERS {
    publishDir "${out_dir}", mode: "$mode"

    input:
        path(adata)

    output:
        path("adata_filtered.h5ad"), emit: adata_filtered

	script:
	"""
    Filter_mito_outliers.py \\
    --adata=${adata}
	"""
}
