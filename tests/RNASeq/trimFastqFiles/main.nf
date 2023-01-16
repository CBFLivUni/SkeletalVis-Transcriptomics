#!/usr/bin/env nextflow

nextflow.enable.dsl = 2


params.accessionNumber = "testDownload"
params.outdir = "data/${params.accessionNumber}"
params.scriptDir = "../../../scripts"
params.single = "True"

include { trimFastqFiles } from '../../../modules/RNASeq/trimFastqFiles'

workflow test_trimFastqFiles {
def input = []
    input = ["SRR11140744",file("ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR111/044/SRR11140744/SRR11140744.fastq.gz") ] 
    trimFastqFiles(input)
}
