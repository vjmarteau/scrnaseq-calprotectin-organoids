process SEURAT_TO_H5AD {
    container "${baseDir}/envs/seuratdisk.sif"
    stageInMode 'link'

    input:
        tuple path(input_seurat)

    output:
        tuple path("*.h5ad"), emit: h5ad

    script:
    """
    #!/usr/bin/env Rscript
    SeuratDisk::Convert("${input_seurat}", dest="h5ad")
    """
}