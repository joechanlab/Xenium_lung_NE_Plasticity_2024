process LABELTRANSFER {
    label 'process_medium'
    container "library://mamie_wang/xenium/scenvi.sif:latest"
    containerOptions "--bind ${params.mount}"
    publishDir "${params.outdir}/LABELTRANSFER/", mode: 'copy'

    input:
    path st_h5ad
    path sc_h5ad
    val name

    output:
    path "${name}_labeled_ST.h5ad", emit: labeled_st_h5ad
    val "${name}", emit: name

    script:
    """
    export NUMBA_CACHE_DIR=${workDir}
    export MPLCONFIGDIR=${workDir}
    python ${baseDir}/bin/label_transfer.py \
        ${st_h5ad} \
        ${sc_h5ad} \
        ${params.labeltransfer.labelfile} \
        ${name}_labeled_ST.h5ad \
        --celltype ${params.labeltransfer.celltype}
    """
}
