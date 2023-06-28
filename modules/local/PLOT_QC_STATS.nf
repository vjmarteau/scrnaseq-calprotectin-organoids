nextflow.enable.dsl=2

out_dir = file(params.resDir)
mode = params.publish_dir_mode

process PLOT_QC_STATS {
    publishDir "${out_dir}", mode: "$mode"

    input:
        path(adata)
        path(adata_filtered)
        path(adata_nodoublet)

    output:
        path("stats.csv"), emit: stats
        path("*_highest_expr_genes.png"), emit: highest_expr_genes, optional: true
        path("*_Violin_total_counts.png"), emit: Violin_total_counts, optional: true
        path("*_Violin_pct_counts_mito.png"), emit: Violin_pct_counts_mito, optional: true

	script:
	"""
    Plot_QC_stats.py \\
    --adata=${adata} \\
    --adata_filtered=${adata_filtered} \\
    --adata_nodoublet=${adata_nodoublet}
	"""
}
