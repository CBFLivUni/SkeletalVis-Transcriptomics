process fastqScreen {

  label 'standard'
  tag "$sampleId"
  publishDir "${params.outdir}/fastqscreen", mode: 'copy'

  input:
  tuple val(sampleId), path(reads)

  output:
  path "*_screen.{txt,html}", emit: stats

  shell:
  '''
  fastq_screen --aligner bowtie2 --threads 1 --conf !{params.fastq_screen_conf} $(echo !{reads}|tr " " "\n"|sort|tr "\n" " ") --force
  '''

}
