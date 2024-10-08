process QCREPORT {
    label 'process_medium'
    container "library://mamie_wang/nf-scrnaseq/postprocessing.sif:latest"
    containerOptions "--bind ${params.mount}"
    publishDir "${params.outdir}/report/", mode: 'copy'

    input:
    path st_h5ad
    val name

    output:
    path "${name}_report.ipynb", emit: report_ipynb
    path "${name}_report.html", emit: report_html

    script:
    """
    export NUMBA_CACHE_DIR=${workDir}
    export MPLCONFIGDIR=${workDir}
    papermill ${baseDir}/bin/QC.ipynb ${name}_report.ipynb \
        -p h5ad ${st_h5ad} \
        -p gene_annotation_file ${params.QC.gene_information_file} \
        -p extra_markers ${params.QC.extra_markers}
    jupyter nbconvert --to html ${name}_report.ipynb
    """
}
