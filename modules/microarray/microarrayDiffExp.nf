process microarrayDiffExp {

      publishDir path: "${params.outdir}", saveAs: { filename -> "${params.accessionNumber}_$filename" } , mode: 'copy'

      input:
      val diffData

      output:
      path "foldChangeTable.txt"


      script:
      """
      Rscript ${params.scriptDir}/microarray/differentialExpressionMicroarray.R --inputfile=$diffData --comparisonsTable="${params.comparisonsTable}" --platform="${params.platform}" --annotationFile="${params.annotationFile}" --foldChangeOnly="${params.foldChangeOnly}" --offset="${params.offset}" --foldChangeTable="foldChangeTable.txt"

      """
}

