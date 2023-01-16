process RNASeqChrDir {

   label 'big_mem'
   publishDir path: "${params.outdir}", saveAs: { filename -> "${params.accessionNumber}_$filename" } , mode: 'copy'

   input:
   val chrDirData

   output:
   path "chrDirTable.txt" optional true

   script:
   """
   Rscript ${params.scriptDir}/RNASeq/characteristicDirectionRNASeq.R --txiData=$chrDirData --sampleTable="${params.sampleTable}" --comparisonsTable="${params.comparisonsTable}" --species="${params.species}" --technicalReplicates="${params.technicalReplicates}" --batchCorrect="${params.batchCorrect}" --supervised="${params.supervised}" --foldChangeOnly="${params.foldChangeOnly}" --chrDirTable="chrDirTable.txt"
   """

}

