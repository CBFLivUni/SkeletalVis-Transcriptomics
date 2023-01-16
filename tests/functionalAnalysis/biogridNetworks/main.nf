#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


params.accessionNumber = "E-GEOD-68760"
params.outdir ="data/${params.accessionNumber}"
params.scriptDir = "../../../scripts"
params.species = "Human"
species="human"
params.foldChangeOnly = "FALSE"
params.referenceDir = "$baseDir/../../../references"
params.biogrid = "$params.referenceDir/networks/${species}BioGridNetwork.txt"

include { biogridNetworks } from '../../../modules/functionalAnalysis/biogridNetworks'

workflow test_biogridNetworks {
 def input = []
    input = [1,file("tests/testData/E-GEOD-68760_cutTable_1.txt") ] 
    

    biogridNetworks(input)
}
