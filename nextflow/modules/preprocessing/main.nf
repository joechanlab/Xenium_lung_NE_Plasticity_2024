process PREPROCESSING {
    label 'process_medium'
    container "library://mamie_wang/nf-scrnaseq/postprocessing.sif:latest"
    publishDir "${params.outdir}/preprocessing/", mode: 'copy'

    input:
    tuple val(name), path(xenium_folder)

    output:
    path "${name}.h5ad", emit: preprocessing_h5ad
    val "${name}", emit: name

    script:
    """
    export NUMBA_CACHE_DIR=${workDir}
    export MPLCONFIGDIR=${workDir}
    python ${baseDir}/bin/preprocessing.py \
        ${xenium_folder} \
        --output ${name}.h5ad
    """
}
