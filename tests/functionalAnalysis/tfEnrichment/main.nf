#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


params.accessionNumber = "E-GEOD-68760"
params.outdir ="data/${params.accessionNumber}"
params.scriptDir = "../../../scripts"
params.species = "Human"
species="human"
params.foldChangeOnly = "FALSE"
params.referenceDir = "$baseDir/../../../references"
params.TF_CHEA3_LIB = "$params.referenceDir/enrichment/TFs.RDS"
params.homology= "$params.referenceDir/homology/Homology.${species}.txt"

include { tfEnrichment } from '../../../modules/functionalAnalysis/tfEnrichment'

workflow test_tfEnrichment {
 def input = []
    input = [1,file("tests/testData/E-GEOD-68760_cutTable_1.txt") ] 

    tfEnrichment(input)
}
