process H5AD_TO_SEURAT {
    container "${baseDir}/envs/seuratdisk.sif"
    stageInMode 'link'

    input:
        tuple path(input_adata)

    output:
        tuple path("*.h5seurat"), emit: h5seurat

    script:
    """
    #!/usr/bin/env Rscript
    SeuratDisk::Convert("${input_adata}", dest="h5seurat")
    """
}