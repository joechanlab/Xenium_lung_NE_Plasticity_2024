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
    if [ ! -d "/tmp/ipykernel" ]; then
        mkdir -p "/tmp/ipykernel"
    fi
    export HOME=/tmp/ipykernel
    python -m ipykernel install --user --name postprocessing
    papermill ${baseDir}/bin/QC.ipynb ${name}_report.ipynb
    jupyter nbconvert --to html ${name}_report.ipynb
    """
}
