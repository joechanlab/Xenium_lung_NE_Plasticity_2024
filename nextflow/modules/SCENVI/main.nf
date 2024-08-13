process SCENVI {
    label 'gpus'
    container "library://mamie_wang/xenium/scenvi.sif:latest"
    containerOptions "--bind ${params.mount}"
    publishDir "${params.outdir}/SCENVI/", mode: 'copy'

    input:
    path st_h5ad
    path sc_h5ad
    val name

    output:
    path "${name}_ENVI_ST.h5ad", emit: envi_st_h5ad
    path "${name}_ENVI_SC.h5ad", emit: envi_sc_h5ad
    path "${name}_ENVI_model.pkl", emit: envi_model_pkl
    path "${name}_ENVI_ST_imputation.pkl", emit: envi_st_imputation_pkl
    val "${name}", emit: name

    script:
    """
    export NUMBA_CACHE_DIR=${workDir}
    export MPLCONFIGDIR=${workDir}
    python ${baseDir}/bin/ENVI.py \
        ${st_h5ad} \
        ${sc_h5ad} \
        ${name}_ENVI_ST.h5ad \
        ${name}_ENVI_SC.h5ad \
        ${name}_ENVI_model.pkl \
        ${name}_ENVI_ST_imputation.h5ad \
        --downsample ${params.SCENVI.downsample} \
        --HVG ${params.SCENVI.HVG}
    """
}
