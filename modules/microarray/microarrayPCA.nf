process microarrayPCA {

    label 'short_job'

    publishDir path: "${params.outdir}", saveAs: { filename -> "${params.accessionNumber}_$filename" } , mode: 'copy'
    
    input:
    val pcaData

    output:
    path "pcaPlot.png"
    path "geneInfluence.txt"
    path "normalisedExp.txt"

    script:
    """
    Rscript ${params.scriptDir}/microarray/PCAMicroarray.R --inputfile=$pcaData --comparisonsTable="${params.comparisonsTable}" --platform="${params.platform}" --annotationFile="${params.annotationFile}" --foldChangeOnly="${params.foldChangeOnly}" --offset="${params.offset}" --pcaPlot="pcaPlot.png" --geneInfluence="geneInfluence.txt" --normalisedExp="normalisedExp.txt"
    """

}

