nextflow.enable.dsl=2

out_dir = file(params.resDir)
mode = params.publish_dir_mode

process EXTRACT_ADATA {
    publishDir "${out_dir}", mode: "$mode"

    input:
        path(adata)

    output:
        tuple path("counts_matrix.mtx"),
        path("denoised_matrix.mtx"),
        path("features.tsv"),
        path("barcodes.tsv"),
        path("metadata.tsv"),
        path("metadata_var.tsv"),
        path("umap.tsv"), emit: convert
    

	script:
	"""
    Convert_extract_adata.py \\
    --adata=${adata}
    
	"""
}

process ASSEMBLE_SEURAT {
    publishDir "${out_dir}", mode: "$mode"
    label "Rscript"

    input:
         tuple path(counts_matrix),
         path(denoised_matrix),
         path(features),
         path(barcodes),
         path(metadata),
         path(metadata_var),
         path(umap)

    output:
        path("*.rds"), emit: seurat

	script:
	"""
    Convert_assemble_seurat.R \\
    --counts_mtx=${counts_matrix} \\
    --denoised_mtx=${denoised_matrix} \\
    --features=${features} \\
    --barcodes=${barcodes} \\
    --metadata=${metadata} \\
    --metadata_var=${metadata_var} \\
    --umap_tsv=${umap}
    
	"""
}

process ASSEMBLE_SCE {
    publishDir "${out_dir}", mode: "$mode"
    label "Rscript"

    input:
         tuple path(counts_matrix),
         path(denoised_matrix),
         path(features),
         path(barcodes),
         path(metadata),
         path(metadata_var),
         path(umap)

    output:
        path("*.rds"), emit: sce

	script:
	"""
    Convert_assemble_sce.R \\
    --counts_mtx=${counts_matrix} \\
    --denoised_mtx=${denoised_matrix} \\
    --features=${features} \\
    --barcodes=${barcodes} \\
    --metadata=${metadata} \\
    --metadata_var=${metadata_var} \\
    --umap_tsv=${umap}
    
	"""
}

workflow CONVERT_SC_FORMAT {
    take:
        adata

    main:
        EXTRACT_ADATA(adata)
        ASSEMBLE_SEURAT(EXTRACT_ADATA.out.convert)
        ASSEMBLE_SCE(EXTRACT_ADATA.out.convert)

    emit:
        seurat = ASSEMBLE_SEURAT.out.seurat
        sce = ASSEMBLE_SCE.out.sce
}
