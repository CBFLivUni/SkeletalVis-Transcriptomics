#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.accessionNumber = "test"
params.outdir ="data/${params.accessionNumber}"



include { fastQC } from '../../../modules/RNASeq/fastQC'

workflow test_fastQC {

  def input = []
    input = ["testSample",file("https://github.com/pachterlab/kallisto/raw/master/test/reads_1.fastq.gz") ] 

    fastQC(input)
}
