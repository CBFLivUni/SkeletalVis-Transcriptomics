#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


params.accessionNumber = "E-GEOD-68760"
params.outdir ="data/${params.accessionNumber}"
params.scriptDir = "../../../scripts"

include { cutFoldChangeTables } from '../../../modules/functionalAnalysis/cutFoldChangeTables'
inputData = file("tests/testData/E-GEOD-68760_foldChangeTable.txt")
numComps = 1

workflow test_cutFoldChangeTables {

    cutFoldChangeTables(inputData,numComps)
}
