process stringNetworks {

  label 'big_mem'
  tag "subnet $num"
  publishDir path: "${params.outdir}", saveAs: { filename -> "${params.accessionNumber}_$filename" } , mode: 'copy'

  input:
  tuple val(num), path (comparison)

  
  output:
  path ("*.{txt,RDS}") optional true

  script:
  """
  Rscript ${params.scriptDir}/functionalAnalysis/GigaSubnetworks.R --differentialExpression=$comparison --interactome="${params.string}" --foldChangeOnly="${params.foldChangeOnly}" --species="${params.species}" --visNetworks="networks_string_${num}.RDS" --summaryTable="summaryTable_string_${num}.txt"

  """

}

