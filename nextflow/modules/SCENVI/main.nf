process SCENVI {
    label 'process_medium'
    container "library://mamie_wang/nf-scrnaseq/ENVI.sif:latest"
    publishDir "${params.outdir}/preprocessing/", mode: 'copy'

    input:
    path xenium_h5ad
    path sc_h5ad
    val name

    output:
    path "${name}_ENVI_ST.h5ad", emit: envi_st_h5ad
    path "${name}_ENVI_SC.h5ad", emit: envi_sc_h5ad
    val "${name}", emit: name

    script:
    """
    export NUMBA_CACHE_DIR=${workDir}
    export MPLCONFIGDIR=${workDir}
    python ${baseDir}/bin/ENVI.py \
        ${xenium_h5ad} \
        ${sc_h5ad} \
        ${name}_ENVI_ST.h5ad \
        ${name}_ENVI_SC.h5ad \
        --downsample 10000 \
        --patient ${name}
    """
}
