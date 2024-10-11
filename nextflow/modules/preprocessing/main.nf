process PREPROCESSING {
    label 'process_medium'
    container "library://mamie_wang/nf-scrnaseq/postprocessing.sif:latest"
    publishDir "${params.outdir}/preprocessing/", mode: 'copy'

    input:
    tuple val(name), path(xenium_folder), val(sc_path)

    output:
    path "${name}.h5ad", emit: preprocessing_h5ad
    val "${sc_path}", emit: sc_h5ad
    val "${name}", emit: name

    script:
    """
    export NUMBA_CACHE_DIR=${workDir}
    export MPLCONFIGDIR=${workDir}
    python ${baseDir}/bin/preprocessing.py \
        ${xenium_folder} \
        ${name}.h5ad \
        --sample_name ${name}
    """
}
