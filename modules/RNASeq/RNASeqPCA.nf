process RNASeqPCA {

    label 'bigmem'

    publishDir path: "${params.outdir}", saveAs: { filename -> "${params.accessionNumber}_$filename" } , mode: 'copy'
    
    input:
    val pcaData

    output:
    path "pcaPlot.png"
    path "geneInfluence.txt"
    path "normalisedExp.txt"

    script:
    """
    Rscript ${params.scriptDir}/RNASeq/PCARNASeq.R --txiData=$pcaData --sampleTable="${params.sampleTable}" --comparisonsTable="${params.comparisonsTable}" --technicalReplicates="${params.technicalReplicates}" --species="${params.species}" --batchCorrect="${params.batchCorrect}" --supervised="${params.supervised}" --pcaPlot="pcaPlot.png" --geneInfluence="geneInfluence.txt" --normalisedExp="normalisedExp.txt"

    """

}

