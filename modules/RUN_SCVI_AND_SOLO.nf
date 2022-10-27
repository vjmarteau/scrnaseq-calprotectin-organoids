nextflow.enable.dsl=2

out_dir = file(params.resDir)
mode = params.publish_dir_mode

process RUN_SCVI_AND_SOLO {
    publishDir "${out_dir}", mode: "$mode"

    input:
        path(adata)

    output:
        path("scVI_model"), emit: scVI_model
        path("adata_nodoublet.h5ad"), emit: adata_nodoublet

	script:
	"""
    Run_scVI_and_SOLO.py \\
    --adata=${adata}
	"""
}
