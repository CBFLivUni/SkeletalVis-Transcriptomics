process SETUP {

    label 'bigmem'

    publishDir path: "${params.referenceDir}", mode: 'move'
    
    output:
    path "homology/*"
    path "networks/*"
    path "transcriptomes/*"
    path "enrichment/*"
    path "fastqscreen/*"

    script:
    """
    #make the subdirectories
    mkdir transcriptomes
    mkdir homology
    mkdir networks
    mkdir enrichment

    #Rscript ${params.scriptDir}/install/makeHomologyTables.R

    #download the fasta and gtf files needed
    #bash ${params.scriptDir}/install/downloadReferences.sh

    #make all the supporting files needed for the pipeline
    
    Rscript ${params.scriptDir}/install/makeBiogridNetwork.R
    Rscript ${params.scriptDir}/install/makeStringNetwork.R
    Rscript ${params.scriptDir}/install/makeReactomePathways.R

    mkdir -p references/fastqscreen
    bash ${params.scriptDir}/install/getFastqScreenGenomes.sh	
    """

}

