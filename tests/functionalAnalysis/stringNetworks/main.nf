#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


params.accessionNumber = "E-GEOD-68760"
params.outdir ="data/${params.accessionNumber}"
params.scriptDir = "../../../scripts"
params.species = "Human"
species="human"
params.foldChangeOnly = "FALSE"
params.referenceDir = "$baseDir/../../../references"
params.string = "$params.referenceDir/networks/${species}StringNetwork.txt"

include { stringNetworks } from '../../../modules/functionalAnalysis/stringNetworks'

workflow test_stringNetworks {
 def input = []
    input = [1,file("tests/testData/E-GEOD-68760_cutTable_1.txt") ] 
    

    stringNetworks(input)
}
