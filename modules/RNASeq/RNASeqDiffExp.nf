process RNASeqDiffExp {

   label 'multi_big_mem'
   publishDir path: "${params.outdir}", saveAs: { filename -> "${params.accessionNumber}_$filename" } , mode: 'copy'

   input:
   path txiData

   output:
   path "foldChangeTable.txt"


   script:
   """
   Rscript ${params.scriptDir}/RNASeq/differentialExpressionRNASeq.R --txiData=$txiData --sampleTable="${params.sampleTable}" --comparisonsTable="${params.comparisonsTable}" --species="${params.species}" --technicalReplicates="${params.technicalReplicates}" --batchCorrect="${params.batchCorrect}" --supervised="${params.supervised}" --foldChangeOnly="${params.foldChangeOnly}" --foldChangeTable="foldChangeTable.txt"

   """

}

