process makeGeneLengths {

	publishDir path: "${params.referenceDir}", mode: 'copy'
	tag "genelengths"
	label 'big_mem'

	input:
	val (geneLengthsFile)

	output: 
	path ("*.txt")
	script:

	"""
	Rscript ${params.scriptDir}/functionalAnalysis/GTFtoGeneLength.R --gtffile="${params.GTF}" --species="${params.species}" --geneExonLengths=$geneLengthsFile
	"""
}
