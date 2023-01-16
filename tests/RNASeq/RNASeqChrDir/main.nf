#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


params.accessionNumber = "testRNASeq"
params.comparisonsTable = "$baseDir/../../testData/testRNASeq_comparisons.txt"
params.sampleTable = "$baseDir/../../testData/testRNASeq_sampleTable.txt"
params.outdir = "data/${params.accessionNumber}"
params.scriptDir = "../../../scripts"
params.species = "Human"
params.foldChangeOnly = "FALSE"
params.technicalReplicates = "FALSE"
params.batchCorrect = "FALSE"
params.supervised="FALSE"

include { RNASeqChrDir } from '../../../modules/RNASeq/RNASeqChrDir'
inputData = file("tests/testData/testRNASeq_txi.RDS")


workflow test_RNASeqChrDir {
    RNASeqChrDir(inputData)
}
