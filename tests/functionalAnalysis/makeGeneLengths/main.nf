#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.species = "Human"
species="human"
params.scriptDir = "../../../scripts"
params.referenceDir = "$baseDir/../../../references"
geneLengthsFile = "${species}GeneLengths.txt"
params.geneLengths= "test_$geneLengthsFile"
params.GTF="$baseDir/../../testData/test_human.gtf"

include { makeGeneLengths } from '../../../modules/functionalAnalysis/makeGeneLengths'

workflow test_makeGeneLengths {

    makeGeneLengths(params.geneLengths)
}
