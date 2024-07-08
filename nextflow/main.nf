#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include {PREPROCESSING} from './modules/preprocessing'

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

    // run Cellbender
    PREPROCESSING(ch_input)
}
