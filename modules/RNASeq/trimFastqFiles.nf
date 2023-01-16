process trimFastqFiles {
  
    label 'multi_cpu'
    tag "$sampleId"

    input:
    tuple val(sampleId), path(reads)

    output:
    tuple val(sampleId), path("trimmed_*")

    shell:
    '''
    python2 !{params.scriptDir}/RNASeq/Trimmomatic.py !{params.single} $(echo !{reads}|tr " " "\n"|sort|tr "\n" " ")
    '''
}
