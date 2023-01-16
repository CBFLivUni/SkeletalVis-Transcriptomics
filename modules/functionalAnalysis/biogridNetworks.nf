process biogridNetworks {

  label 'big_mem'
  tag "subnet $num"
  publishDir path: "${params.outdir}", saveAs: { filename -> "${params.accessionNumber}_$filename" } , mode: 'copy'

  input:
  tuple val(num), path (comparison)

  
  output:
  path ("*.{txt,RDS,pdf}") optional true

  script:
  """
  Rscript ${params.scriptDir}/functionalAnalysis/GigaSubnetworks.R --differentialExpression=$comparison --interactome="${params.biogrid}" --foldChangeOnly="${params.foldChangeOnly}" --species="${params.species}" --visNetworks="networks_biogrid_${num}.RDS" --summaryTable="summaryTable_biogrid_${num}.txt"

  """

}

