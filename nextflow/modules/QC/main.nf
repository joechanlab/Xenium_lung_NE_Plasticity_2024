process QCREPORT {
    label 'process_low'
    container "library://mamie_wang/nf-scrnaseq/postprocessing.sif:latest"
    publishDir "${params.outdir}/report/", mode: 'copy'

    input:
    path preprocessing_h5ad
    val name

    output:
    path "${name}_report.html", emit: report_html

    script:
    """
    export NUMBA_CACHE_DIR=/tmp/numba_cache
    papermill ${baseDir}/bin/QC.ipynb ${name}_report.ipynb \
        -p h5ad ${preprocessing_h5ad} \
        -p gene_annotation_file ${params.gene_information_file}
    jupyter nbconvert --to html ${name}_report.ipynb
    """
}
