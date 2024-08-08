#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {PREPROCESSING} from './modules/preprocessing'
include {SVG} from './modules/SVG'
include {SCENVI} from './modules/SCENVI'
include {QCREPORT} from './modules/QC_report'

workflow {
    // access the samplesheet
    sample_sheet = file(params.samplesheet)

    // read the sample sheet as a CSV
    sample_sheet_data = sample_sheet.text.readLines().drop(1).collect { it.split(',') }

    // create a channel from the paths
    ch_input = Channel.from(sample_sheet_data).map { row ->
        def name = row[0]
        def st_path = file(row[1])
        def sc_path = file(row[2])
        return tuple(name, st_path, sc_path)
    }

    // run preprocessing
    PREPROCESSING(ch_input)

    // run SVG
    SVG(PREPROCESSING.out.preprocessing_h5ad, 
        PREPROCESSING.out.sc_h5ad,
        PREPROCESSING.out.name
        )

    // run ENVI
    SCENVI(SVG.out.SVG_h5ad, 
        SVG.out.sc_h5ad,
        SVG.out.name
        )

    // generate QC report
    QCREPORT(SCENVI.out.envi_st_h5ad, SCENVI.out.name)
}
