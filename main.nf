//!/usr/bin/env nextflow
nextflow.enable.dsl=2


def helpMessage() {
    log.info"""
    =========================================
     Skeletalvis-Transcriptomics
    =========================================
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow run main.nf -profile slurm -params-file params/<accessionNumber>.yaml -with-singularity library://jsoul/default/skeletalvis-pipeline -entry <setup,skeletalvis>

    Required arguments:
         -profile                      Configuration profile to use. <local, slurm>
         -params-file                  Yaml file containing parameters for the analysis
         --with-singularity            Recommend to run the pipeline within the provided singularity container
         -entry			       Workflow to run for inital setup or analysis

""".stripIndent()
}

params.help = false
if (params.help){
    helpMessage()
    exit 0
}



//define the parameters

//outdirectory for the results
params.outdir = "$baseDir/data/${params.accessionNumber}"

//site to download the sra/fastq files if needed
params.downloadSite = "ENA"

//accession number for ENA
params.ENA = ""

//directory with the fastqFiles
//ignoring R0 fastq files(usually discarded reads in trimming) by default
params.fastqFileDir = "$params.outdir/fastqFiles/*[!_R0].fastq.gz"
downloadFiles = file(params.fastqFileDir).isEmpty()


//make sure the species is lower case to match the downloaded reference data
params.species = "Human"
species = params.species.toLowerCase()

//define the ensembl release used for the fasta and the gtf file
params.ensemblRelease = "103"


params.referenceDir = "$baseDir/references"

//GTF file to generate the gene lengths
params.GTF = "$params.referenceDir/transcriptomes/${species}_release${params.ensemblRelease}.gtf"

//fasta file for generating the kallisto index
params.fasta = "$params.referenceDir/transcriptomes/${species}_release${params.ensemblRelease}.fa.gz"

//directory with the scripts
params.scriptDir = "$baseDir/scripts"

//directory with the fastqscreen genomes
params.fastq_screen_conf = "$params.referenceDir/fastqscreen/FastQ_Screen_Genomes/fastq_screen.conf"

//index for kallisto
params.indexName = "${params.species}_release${params.ensemblRelease}IndexFile"
params.index = "$params.referenceDir/$params.indexName"

//stranded library option for kallisto?
params.stranded = false

//should trimming be skipped?
params.skipTrimming = false

//sample table for the experiment
params.sampleTable = "$baseDir/params/sampleTables/${params.accessionNumber}_sampleTable.txt"

//comparison table for the experiment
params.comparisonsTable = "$baseDir/params/comparisons/${params.accessionNumber}_comparisons.txt"

//number of cpu cores to use for sra to fastq parsing. Check nextflow.config for local cpu executor usage.
params.cpuCores = 1

//does the differential expression analysis need supervised batch effect correction (uses svaseq in supervised mode)?
params.supervised = "FALSE"

//protein-protein interaction network files
params.biogrid = "$params.referenceDir/networks/${species}BioGridNetwork.txt"
params.string = "$params.referenceDir/networks/${species}StringNetwork.txt"

//enrichment libraries
params.TF_CHEA3_LIB = "$params.referenceDir/enrichment/TFs.RDS"
params.reactomePathways = "$params.referenceDir/enrichment/reactomePathways${params.species}.RDS"

//gene length files - will create from a provided GTF file if not provided
geneLengthsFile = "${species}GeneLengths.txt"
params.geneLengths= "$params.referenceDir/$geneLengthsFile"

//other species to human homology files for mapping genes to pathway databases
params.homology= "$params.referenceDir/homology/Homology.${species}.txt"


//covariates to adjust for in the linear modelling
params.covariates=""
params.covariate_types=""


//thresholds to interate through to define differentially expressed genes for the downstream analysis
params.foldchange = [1,1.5,2]

params.foldChangeOnly = "TRUE"

if (params.foldChangeOnly=="FALSE"){
params.padj = [1,0.05,0.01]
} else {
params.padj = [1]
}

//microarray workflow specific parameters

params.localData = false
params.grepExpression = "FALSE"
params.grepString = ""
params.numberLinesSkip = 0
params.numberFileRemove = 0
params.split = "FALSE" 
params.splitField = ""
params.splitSep = ""
params.splitPos = ""
params.remove = "FALSE"
params.removeSample = ""
params.offset = 80
params.newColumns = ""

params.setup = false


params.asperaKey = "$params.referenceDir/asperaweb_id_dsa.openssh"

params.txiDataPath = "$params.outdir/${params.accessionNumber}_txi.RDS"

//setup module
include { SETUP } from './modules/install/setup'


//Microarray modules
include { getGeneExpressionData } from './modules/microarray/getGeneExpressionData'
include { microarrayQC } from './modules/microarray/microarrayQC'
include { microarrayChrDir } from './modules/microarray/microarrayChrDir'
include { microarrayPCA } from './modules/microarray/microarrayPCA'
include { microarrayDiffExp } from './modules/microarray/microarrayDiffExp'

//RNA-seq modules
include { sraDownload } from './modules/RNASeq/sraDownload'
include { enaDownload } from './modules/RNASeq/enaDownload'
include { fastqToSampleID } from './modules/RNASeq/fastqToSampleID'
include { trimFastqFiles } from './modules/RNASeq/trimFastqFiles'
include { fastQC } from './modules/RNASeq/fastQC'
include { fastqScreen } from './modules/RNASeq/fastqScreen'
include { kallistoQuant } from './modules/RNASeq/kallistoQuant'
include { RNASeqPCA } from './modules/RNASeq/RNASeqPCA'
include { RNASeqChrDir } from './modules/RNASeq/RNASeqChrDir'
include { RNASeqDiffExp } from './modules/RNASeq/RNASeqDiffExp'
include { multiQC } from './modules/RNASeq/multiQC'
include { abundancesToTxi } from './modules/RNASeq/abundancesToTxi'
include { kallistoIndex } from './modules/RNASeq/kallistoIndex'


//Functional analysis modules
include { goEnrichment } from './modules/functionalAnalysis/goEnrichment'
include { tfEnrichment } from './modules/functionalAnalysis/tfEnrichment'
include { cutFoldChangeTables } from './modules/functionalAnalysis/cutFoldChangeTables'
include { pathwayEnrichment } from './modules/functionalAnalysis/pathwayEnrichment'
include { biogridNetworks } from './modules/functionalAnalysis/biogridNetworks'
include { stringNetworks } from './modules/functionalAnalysis/stringNetworks'
include { volcanoPlot } from './modules/functionalAnalysis/volcanoPlot'
include { makeGeneLengths } from './modules/functionalAnalysis/makeGeneLengths'

if(!params.setup){
  //number of comparisons in the differential expression analysis
  params.numComps=file(params.comparisonsTable).countLines()-1
  numberComps= (1..params.numComps)
}


//set up the directory with the required data
workflow setup {

   SETUP()    
}


//main entry point for the transcriptomics analysis
workflow skeletalvis {

	if(params.Type=="Microarray") {

	MicroarrayData()

	} else {

	RNASeqData()
	}

}




//subworkflow to perform pathway, GO and transcription factor analysis
workflow functionalAnalysis {


	take: diffExpTable
	      geneLengthFile
	      padj
	      foldchange
	      
	main:
		goEnrichment(diffExpTable,geneLengthFile,padj,foldchange)
		pathwayEnrichment(diffExpTable,geneLengthFile,padj,foldchange)
		tfEnrichment(diffExpTable)
		biogridNetworks(diffExpTable)
		stringNetworks(diffExpTable)
		volcanoPlot(diffExpTable)

}

//subworkflow to download the fastq files from ENA/SRA if needed
workflow downloadRNASeqData {

 	main:
  	if(params.downloadSite=="ENA") {
  	fastqFiles = enaDownload()
  	} else {
  	fastqFiles = sraDownload()
  	}


  emit:
    fastqFiles


}


//workflow to perform RNA-seq analysis including differential expression and functional analysis
workflow RNASeqData {

  if(downloadFiles){

    rawData = downloadRNASeqData().fastqFiles.flatten()

  } else {

    rawData = channel.fromPath( params.fastqFileDir).view()

  }

  //group the data if paired reads
  groupedReads = rawData.map { file ->
          def key = file.name.toString().split("_R\\d*")[0]
          return tuple(key, file)
       }
      .groupTuple(sort:true).view()
          

  if(!file(params.geneLengths).exists()){

  	geneLengths = makeGeneLengths(geneLengthsFile)

  	} else {

    geneLengths = file(params.geneLengths)
  }


  if(!file(params.index).exists()){

  	index = kallistoIndex(params.indexName,params.fasta)

  } else {
    index = file(params.index)
  }

 
  groupedSamples = fastqToSampleID(groupedReads,params.sampleTable)
    .groupTuple()
    .map{id,fastq -> tuple(id,fastq.flatten().sort())}.view()

	if(params.skipTrimming){
		trimmedSamples = groupedSamples
	} else {
		trimmedSamples = trimFastqFiles(groupedSamples)
	}

  fastQC(trimmedSamples)
  fastqScreen(trimmedSamples)
  trimmedSamples.view()
  kallistoQuant(trimmedSamples,index)
  multiQC(fastQC.out.stats.collect(),fastqScreen.out.stats.collect(),kallistoQuant.out.stats.collect())

  abundanceMatrix = abundancesToTxi(kallistoQuant.out.abundances.collect())

  RNASeqPCA(abundancesToTxi.out.txiData)

  diffTables = RNASeqDiffExp(abundancesToTxi.out.txiData)
  cutTables = cutFoldChangeTables(diffTables,numberComps)

}


//workflow to perform microarray analysis including differential expression and functional analysis
workflow MicroarrayData {

  if(!file(params.geneLengths).exists()){

 		geneLengths = makeGeneLengths(geneLengthsFile)

 	} else {

    geneLengths = file(params.geneLengths)
 }

	if( params.localData ){
	  rawData = Channel
		.fromPath("$baseDir/params/localData/${params.accessionNumber}.RDS")
		.ifEmpty { exit 1, "local data not found: $baseDir/params/localData/${params.accessionNumber}.RData" }
	} else {
		rawData = getGeneExpressionData()
	}
 
   microarrayQC(rawData)

  microarrayPCA(rawData)
  diffTables = microarrayDiffExp(rawData)
  cutTables = cutFoldChangeTables(microarrayDiffExp.out.diffTables,numberComps)

}
