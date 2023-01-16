process volcanoPlot {

    label 'short_big_mem'
    tag "volcano $num"
    publishDir path: "${params.outdir}", saveAs: { filename -> "${params.accessionNumber}_$filename" } , mode: 'copy'

    input:
    tuple val(num), path (comparison)
    
    output:
    path ("*.png") optional true

    script:

    if (params.foldChangeOnly=="FALSE") {
        """
        Rscript ${params.scriptDir}/functionalAnalysis/volcanoPlot.R --differentialExpression=$comparison --foldChangeOnly="${params.foldChangeOnly}" --volcanoOutput="volcanoPlot_${num}.png"

        """
        }  else {

           """
           echo "Not creating a volcano plot as no significance values"
           """
       }

   }

