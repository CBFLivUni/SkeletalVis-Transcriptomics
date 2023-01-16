process cutFoldChangeTables {

    label 'short_job'

    tag "cutTable"
    publishDir path: "${params.outdir}", saveAs: { filename -> "${params.accessionNumber}_$filename" } , mode: 'copy'

    input:
    path foldChange
    each numberComps
    
    
    output:
    tuple val (numberComps), path ("*.txt")

    script:
    """
    Rscript ${params.scriptDir}/functionalAnalysis/cutFoldChangeTable.R --foldChangeTable=$foldChange --comparisonNumber $numberComps --cutTable=cutTable_${numberComps}.txt

    """
}

