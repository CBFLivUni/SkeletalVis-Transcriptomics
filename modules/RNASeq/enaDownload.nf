process enaDownload {
  label 'download'
  publishDir path: "${params.outdir}/fastqFiles", mode: 'copy'

  output:
  path ('*.fastq.gz'), emit: fastqFiles


  script:
  
  """
  Rscript ${params.scriptDir}/RNASeq/getRNASeqExpressionData.R --accessionNumber "${params.ENA}" --sampleTable "${params.sampleTable}" --key ${params.asperaKey} --aspera "TRUE"
  """

}

