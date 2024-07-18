process SVG {
    label 'process_medium'
    container "library://mamie_wang/nf-scrnaseq/squidpy.sif:latest"
    publishDir "${params.outdir}/SVG/", mode: 'copy'

    input:
    path preprocessing_h5ad
    val name

    output:
    path "${name}_SVG.h5ad", emit: SVG_h5ad
    val "${name}", emit: name

    script:
    """
    export NUMBA_CACHE_DIR=${workDir}
    export MPLCONFIGDIR=${workDir}
    python ${baseDir}/bin/SVG.py \
        ${preprocessing_h5ad} \
        --output ${name}_SVG.h5ad
    """
}
