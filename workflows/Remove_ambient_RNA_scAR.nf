nextflow.enable.dsl=2

out_dir = file(params.outdir)
mode = params.publish_dir_mode

process LOAD_ADATA {
    publishDir "${out_dir}", mode: "$mode"
    label "extra_cpus"

    input:
        path(samplesheet)
        path(cellranger_path)

    output:
        path("samplesheet.csv"), emit: samplesheet
        path("*_filtered_adata.h5ad"), emit: filtered_adata
        path("*_raw_adata.h5ad"), emit: raw_adata
        path("var_names.csv"), emit: var_names

	script:
	"""
    scAR_load_adata.py \\
    --cellranger_path=${cellranger_path} \\
    --samplesheet=${samplesheet} \\
    --cpus=${task.cpus}
	"""
}

process RUN_SCAR {
    publishDir "${out_dir}", mode: "$mode"
    label "gpu"

    input:
        tuple val(id), path(filtered_adata), path(raw_adata)

    output:
        tuple val(id), path("${id}_denoised_adata.h5ad"), emit: denoised_adata

	script:
	"""
    scAR_denoise_adata.py \\
    --raw_adata=${raw_adata} \\
    --filtered_adata=${filtered_adata} \\
    --cpus=${task.cpus} \\
    --output_file="${id}_denoised_adata.h5ad"
	"""
}

process CONCAT_ADATAS {
    publishDir "${out_dir}", mode: "$mode"
    label "extra_cpus"

    input:
        path(adatas_path)
        path(samplesheet)
        path(var_names)
        path(gtf_file)
        path(hgnc_file)

    output:
        path("denoised_adata.h5ad"), emit: denoised_adata

	script:
	"""
    scAR_concat_adatas.py \\
    --adatas_path=. \\
    --samplesheet=${samplesheet} \\
    --var_names=${var_names} \\
    --gtf_file=${gtf_file} \\
    --hgnc_file=${hgnc_file} \\
    --cpus=${task.cpus}
	"""
}

workflow Remove_ambient_RNA_scAR {
    take:
        samplesheet
        ch_input_files
        gtf_file
        hgnc_file

    main:
        LOAD_ADATA(samplesheet, ch_input_files)

        ch_filtered_adata = LOAD_ADATA.out.filtered_adata.flatten().map { 
        it -> [it.baseName.replace("_filtered_adata", ""), it]
        }

        ch_raw_adata = LOAD_ADATA.out.raw_adata.flatten().map { 
            it -> [it.baseName.replace("_raw_adata", ""), it]
        }
        
        ch_scAR_input = ch_filtered_adata.join(ch_raw_adata)

        RUN_SCAR(ch_scAR_input)

        ch_concat_adatas = RUN_SCAR.out.denoised_adata.map{ id, path -> path }.collect()

        CONCAT_ADATAS(ch_concat_adatas, LOAD_ADATA.out.samplesheet, LOAD_ADATA.out.var_names, gtf_file, hgnc_file)
        
    emit:
        denoised_adata = CONCAT_ADATAS.out.denoised_adata

}