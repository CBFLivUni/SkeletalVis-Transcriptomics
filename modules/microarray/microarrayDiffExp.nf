process microarrayDiffExp {

      publishDir path: "${params.outdir}", saveAs: { filename -> "${params.accessionNumber}_$filename" } , mode: 'copy'

      input:
      val diffData

      output:
      path "sampleSizes.txt"
      path "foldChangeTable.txt", emit: diffTables

      script:
      """
      Rscript ${params.scriptDir}/microarray/differentialExpressionMicroarray.R --inputfile=$diffData --comparisonsTable="${params.comparisonsTable}" --platform="${params.platform}" --annotationFile="${params.annotationFile}" --foldChangeOnly="${params.foldChangeOnly}" --offset="${params.offset}" --covariates="${params.covariates}" --covariate_types="${params.covariate_types}" --foldChangeTable="foldChangeTable.txt"

      """
}

