#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.accessionNumber = "test"
params.outdir ="data/${params.accessionNumber}"
params.referenceDir = "$baseDir/../../../references"
params.fastq_screen_conf = "$params.referenceDir/fastqscreen/FastQ_Screen_Genomes/fastq_screen.conf"


include { fastScreen } from '../../../modules/RNASeq/fastScreen'

workflow test_fastScreen {

  def input = []
    input = ["testSample",file("https://github.com/pachterlab/kallisto/raw/master/test/reads_1.fastq.gz") ] 

    fastScreen(input)
}
