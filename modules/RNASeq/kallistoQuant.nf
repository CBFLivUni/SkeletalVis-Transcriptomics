process kallistoQuant {

  tag "$sampleId"
  publishDir path: "${params.outdir}/kallisto", mode: 'copy'

  label 'multi_cpu'

  input:
  tuple val(sampleId), path(reads)
  path index

  output:

  path ("${sampleId}/*abundances"), emit: abundances
  path ("${sampleId}/*.log"), emit: stats

  shell:
  stranded = ""
  if (params.stranded) {
     '''
     stranded = "--fr-stranded"
     '''
  }

  if (params.single) {
     '''
     mkdir !{sampleId}
     kallisto quant --plaintext -i !{index} --single -l 200 -s 20 -t 12 -o !{sampleId} $(echo !{reads}|tr " " "\n"|sort|tr "\n" " ") &>  !{sampleId}/!{sampleId}_kallisto.log
     mv !{sampleId}/abundance.tsv !{sampleId}/!{sampleId}.abundances

     '''
     } else {
       '''
       mkdir !{sampleId}
       kallisto quant --plaintext -i !{index} !{stranded} -t 12 -o !{sampleId} $(echo !{reads}|tr " " "\n"|sort|tr "\n" " ") &> !{sampleId}/!{sampleId}_kallisto.log
       mv !{sampleId}/abundance.tsv !{sampleId}/!{sampleId}.abundances
       '''
    }

 }
