process getGeneExpressionData {

    label 'download'

    publishDir path: "${params.outdir}/rawMicroarrayData", mode: 'copy'

    output:
    path ('*.RDS')

    script:
    
    """
    Rscript ${params.scriptDir}/microarray/getMicroarrayExpressionData.R --accessionNumber "${params.accessionNumber}" --platform "${params.platform}" --numberFileRemove "${params.numberFileRemove}" --grepExpression "${params.grepExpression}" --grepString "${params.grepString}" --numberLinesSkip "${params.numberLinesSkip}" --split "${params.split}" --splitField "${params.splitField}"  --splitSep "${params.splitSep}"  --splitPos "${params.splitPos}" --remove "${params.remove}" --removeSample "${params.removeSample}" --site "${params.site}" --expressionData "${params.accessionNumber}.RDS"
    """

}

