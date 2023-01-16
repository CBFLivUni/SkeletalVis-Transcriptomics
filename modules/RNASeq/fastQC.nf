process fastQC {
  
  label 'standard'
  tag "$sampleId"
  publishDir "${params.outdir}/fastqc", mode: 'copy',
  saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

  input:
  tuple val(sampleId), path(reads)

  output:
  path "*_fastqc.{zip,html}", emit: stats

  shell:
  '''
  fastqc -t 1 -q $(echo !{reads}|tr " " "\n"|sort|tr "\n" " ")
  '''

}
