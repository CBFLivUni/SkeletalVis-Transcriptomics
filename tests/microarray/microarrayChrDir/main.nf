#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


params.accessionNumber = "E-GEOD-68760"
params.comparisonsTable = "$baseDir/../../testData/E-GEOD-68760_comparisons.txt"
params.outdir ="data/${params.accessionNumber}"
params.scriptDir = "../../../scripts"
params.platform = "Affy"
params.site = "ArrayExpress"
params.annotationFile = "primeview.db"
params.foldChangeOnly = "FALSE"
params.offset = 0


include { microarrayChrDir } from '../../../modules/microarray/microarrayChrDir'
inputData = file("tests/testData/E-GEOD-68760.RData")


workflow test_microarrayChrDir {

    microarrayChrDir(inputData)
}
