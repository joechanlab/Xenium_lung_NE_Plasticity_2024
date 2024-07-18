#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {PREPROCESSING} from './modules/preprocessing'
include {SVG} from './modules/SVG'
include {QCREPORT} from './modules/QC'

workflow {
    // access the samplesheet
    sample_sheet = file(params.samplesheet)

    // read the sample sheet as a CSV
    sample_sheet_data = sample_sheet.text.readLines().drop(1).collect { it.split(',') }

    // create a channel from the paths
    ch_input = Channel.from(sample_sheet_data).map { row ->
        def name = row[0]
        def xenium_path = file(row[1])
        return tuple(name, xenium_path)
    }

    // run preprocessing
    PREPROCESSING(ch_input)

    // run SVG
    SVG(PREPROCESSING.out.preprocessing_h5ad, PREPROCESSING.out.name)

    // generate QC report
    QCREPORT(SVG.out.SVG_h5ad, SVG.out.name)
}
