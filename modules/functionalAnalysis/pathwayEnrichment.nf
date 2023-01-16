process pathwayEnrichment {
	label 'short_big_mem'

	tag "pathways $num $padj $foldchange"
	publishDir path: "${params.outdir}", saveAs: { filename -> "${params.accessionNumber}_$filename" } , mode: 'copy'
	errorStrategy { task.exitStatus == 3 ? 'ignore' : 'terminate' }

	input:
	tuple val(num), path (comparison)
	path geneLengths
	each padj
	each foldchange

	output:
	file ("*.txt") optional true

	script:

	if (!(padj==1 & foldchange==1)){
		if (params.foldChangeOnly=="TRUE") {

			"""
			Rscript ${params.scriptDir}/functionalAnalysis/pathwayEnrichment.R --differentialExpression=$comparison --reactomePathways="${params.reactomePathways}" --geneExonLengths=$geneLengths --foldChangeOnly="${params.foldChangeOnly}" --foldchange=$foldchange --padj=$padj --enrichedPathways="pathways_${num}_${foldchange}.txt"

			"""

			} else{
				"""
				Rscript ${params.scriptDir}/functionalAnalysis/pathwayEnrichment.R --differentialExpression=$comparison --reactomePathways="${params.reactomePathways}" --geneExonLengths=$geneLengths --foldChangeOnly="${params.foldChangeOnly}" --foldchange=$foldchange --padj=$padj --enrichedPathways="pathways_${num}_${foldchange}_${padj}.txt"

				"""
			}
			}  else {

				"""
				echo "excluding parameter combination"
				"""
			}

		}

