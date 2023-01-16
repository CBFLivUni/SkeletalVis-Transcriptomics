process tfEnrichment {

  label 'short_big_mem'
  tag "chea3 $num"
  publishDir path: "${params.outdir}", saveAs: { filename -> "${params.accessionNumber}_$filename" } , mode: 'copy'

  input:
  tuple val(num), path (comparison)

  output:
  path ("*.txt")

  script:
  """
  Rscript ${params.scriptDir}/functionalAnalysis/chea3.R --differentialExpression=$comparison --homology="${params.homology}" --foldChangeOnly="${params.foldChangeOnly}" --species="${params.species}" --TFs="${params.TF_CHEA3_LIB}" --TFUp="CHEA3_Up_${num}.txt" --TFDown="CHEA3_Down_${num}.txt"

  """

}

