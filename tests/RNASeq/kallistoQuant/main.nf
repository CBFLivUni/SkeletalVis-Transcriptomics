#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.accessionNumber = "test"
params.outdir ="data/${params.accessionNumber}"
params.referenceDir = "$baseDir/../../../references"
params.stranded = false
params.single = true


include { kallistoQuant } from '../../../modules/RNASeq/kallistoQuant'

workflow test_kallistoQuant {

  def input = []
    input = ["testSample",file("https://github.com/pachterlab/kallisto/raw/master/test/reads_1.fastq.gz") ] 

    kallistoQuant(input,file("tests/testData/testIndex.idx"))
}
