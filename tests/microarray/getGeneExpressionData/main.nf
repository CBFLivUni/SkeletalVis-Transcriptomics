#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


params.accessionNumber = "E-GEOD-68760"
params.outdir ="data/${params.accessionNumber}"
params.scriptDir = "../../../scripts"
params.platform = "Affy"
params.site = "ArrayExpress"
params.grepExpression = "FALSE"
params.grepString = ""
params.numberLinesSkip = 0
params.numberFileRemove = 0
params.split = "FALSE" 
params.splitField = ""
params.splitSep = ""
params.splitPos = ""
params.remove = "FALSE"
params.removeSample = ""
params.offset = 0

include { getGeneExpressionData } from '../../../modules/microarray/getGeneExpressionData'

workflow test_getGeneExpressionData {

    getGeneExpressionData ()
}
