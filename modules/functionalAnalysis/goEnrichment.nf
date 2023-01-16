process goEnrichment {

  label 'standard'


  tag "GO $num $padj $foldchange"
  publishDir path: "${params.outdir}", saveAs: { filename -> "${params.accessionNumber}_$filename" } , mode: 'copy'
  errorStrategy { task.exitStatus == 3 ? 'ignore' : 'terminate' }
  
  input:
  tuple val(num), path (comparison)
  path geneLengths
  each padj
  each foldchange

  output:
  file ("*.{txt,html}") optional true

  script:

  if (!(padj==1 & foldchange==1)){
    if (params.foldChangeOnly=="TRUE") {
       """
       Rscript ${params.scriptDir}/functionalAnalysis/GOEnrichment.R --differentialExpression=$comparison --geneLengths=$geneLengths --foldChangeOnly="${params.foldChangeOnly}" --species="${params.species}" --foldchange=$foldchange --padj=$padj --enrichedTerms="goterms_${num}_${foldchange}.txt" --enrichedTermsReduced="goterms_reduced_${num}_${foldchange}.txt" --mdsPlot="GO.MDS_${num}_${foldchange}.html"

       """
       } else{
        """
        Rscript ${params.scriptDir}/functionalAnalysis/GOEnrichment.R --differentialExpression=$comparison --geneLengths=$geneLengths --foldChangeOnly="${params.foldChangeOnly}" --species="${params.species}" --foldchange=$foldchange --padj=$padj --enrichedTerms="goterms_${num}_${foldchange}_${padj}.txt" --enrichedTermsReduced="goterms_reduced_${num}_${foldchange}_${padj}.txt" --mdsPlot="GO.MDS_${num}_${foldchange}_${padj}.html"
        """
     }
     }  else {

        """
        echo "excluding parameter combination"
        """
     }

  }

