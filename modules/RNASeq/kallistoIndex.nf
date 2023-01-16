process kallistoIndex {
  publishDir path: "${params.referenceDir}", mode: 'copy'
  tag "kallisto index"
  label 'big_mem'

  input:

  val(indexName)
  path(fasta)

  output:

  file indexName

  script:

  """
  kallisto index -i $indexName $fasta

  """
}
