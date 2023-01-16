#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


params.referenceDir = "$baseDir/../../../references"

include { kallistoIndex } from '../../../modules/RNASeq/kallistoIndex'

workflow test_kallistoIndex {

    kallistoIndex("testIndex.idx",file("https://github.com/pachterlab/kallisto/raw/master/test/transcripts.fasta.gz"))
}
