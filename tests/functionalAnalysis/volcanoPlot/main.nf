#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


params.accessionNumber = "E-GEOD-68760"
params.outdir ="data/${params.accessionNumber}"
params.scriptDir = "../../../scripts"
params.foldChangeOnly = "FALSE"

include { volcanoPlot } from '../../../modules/functionalAnalysis/volcanoPlot'

workflow test_volcanoPlot {
 def input = []
    input = [1,file("tests/testData/E-GEOD-68760_cutTable_1.txt") ] 

    volcanoPlot(input)
}
