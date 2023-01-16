process sraDownload {

    label 'multi_cpu'

    publishDir path: "${params.outdir}/fastqFiles", mode: 'copy'

    input:
    path sampleTable

    output:
    path ('*.fastq.gz'), emit: fastqFiles

    script:
    
    """
    python3 ${params.scriptDir}/RNASeq/getRNASeqExpressionDataFromSRA.py $sampleTable "${params.cpuCores}"
    """

}

