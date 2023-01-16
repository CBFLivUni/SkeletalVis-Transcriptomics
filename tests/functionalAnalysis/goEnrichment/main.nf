#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


params.accessionNumber = "E-GEOD-68760"
params.outdir ="data/${params.accessionNumber}"
params.scriptDir = "../../../scripts"
params.species = "Human"
species="human"
params.foldChangeOnly = "FALSE"
params.referenceDir = "$baseDir/../../../references"
params.homology= "$params.referenceDir/homology/Homology.${species}.txt"
geneLengthsFile = "${species}GeneLengths.txt"
params.geneLengths= "$params.referenceDir/$geneLengthsFile"


include { goEnrichment } from '../../../modules/functionalAnalysis/goEnrichment'

workflow test_goEnrichment{
 def input = []
    input = [1,file("tests/testData/E-GEOD-68760_cutTable_1.txt") ] 
    

    goEnrichment(input,file(params.geneLengths),0.05,[1.5,1])
}
