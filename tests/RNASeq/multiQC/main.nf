#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.accessionNumber = "testmultiQC"
params.outdir ="data/${params.accessionNumber}"


include { multiQC } from '../../../modules/RNASeq/multiQC'

workflow test_multiQC {
   
    multiQC(file('tests/testData/multiQC'),
    file('tests/testData/multiQC'),
    file('tests/testData/multiQC'))
}
