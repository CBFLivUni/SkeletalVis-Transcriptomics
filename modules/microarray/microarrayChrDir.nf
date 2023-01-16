process microarrayChrDir {

   label 'big_mem'
   publishDir path: "${params.outdir}", saveAs: { filename -> "${params.accessionNumber}_$filename" } , mode: 'copy'

   input:
   val chrDirData

   output:
   path "chrDirTable.txt" optional true

   script:
   """
   Rscript ${params.scriptDir}/microarray/characteristicDirectionMicroarray.R --inputfile=$chrDirData --comparisonsTable="${params.comparisonsTable}" --platform="${params.platform}" --annotationFile="${params.annotationFile}" --foldChangeOnly="${params.foldChangeOnly}" --offset="${params.offset}" --chrDirTable="chrDirTable.txt"
   """

}

