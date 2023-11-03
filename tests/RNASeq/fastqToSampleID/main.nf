#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


params.accessionNumber = "testfastqToSampleID"
params.sampleTable = "$baseDir/../../testData/testDownload_sampleTable.txt"
params.scriptDir = "../../../scripts"

include { fastqToSampleID } from '../../../modules/RNASeq/fastqToSampleID'



workflow test_fastqToSampleID {
    def input = []
    input = ["SRR11140744.fastq.gz",file("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR111/044/SRR11140744/SRR11140744.fastq.gz")]
    out = fastqToSampleID(input,params.sampleTable).view()
}
