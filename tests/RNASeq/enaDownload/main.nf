#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


params.accessionNumber = "testDownload"
params.ENA = "PRJNA607948"
params.sampleTable = "$baseDir/../../testData/testDownload_sampleTable.txt"
params.outdir = "data/${params.accessionNumber}"
params.scriptDir = "../../../scripts"

include { enaDownload } from '../../../modules/RNASeq/enaDownload'

workflow test_enaDownload {
    enaDownload()
}
