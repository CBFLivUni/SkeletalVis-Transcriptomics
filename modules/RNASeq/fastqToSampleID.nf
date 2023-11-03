process fastqToSampleID {

  label 'very_short_job'

  input:
  tuple val(name), path(sample)
  path sampleTable

  output:
  tuple stdout , path(sample)

  script:
  """
  python ${params.scriptDir}/RNASeq/fastqToSampleID.py $sampleTable $name > out.txt

  if [ -s "out.txt" ]; then
    cat out.txt
  else
    echo "$sample not matched"
    exit 1
  fi	
  """
}

