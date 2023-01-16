#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


params.accessionNumber = "testDownload"
params.sampleTable = "$baseDir/../../testData/testDownload_sampleTable.txt"
params.outdir = "data/${params.accessionNumber}"
params.scriptDir = "../../../scripts"
params.cpuCores = 1

include { sraDownload } from '../../../modules/RNASeq/sraDownload'

workflow test_sraDownload {
    sraDownload(file(params.sampleTable))
}
