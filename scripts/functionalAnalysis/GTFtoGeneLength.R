#!/usr/bin/env Rscript

 
  suppressPackageStartupMessages(library("optparse"))
  
  option_list <- list()
  
  option_list$gtffile <- make_option('--gtffile', type='character')
  option_list$species <- make_option('--species', type='character')
  option_list$geneExonLengths <- make_option('--geneExonLengths', type='character')
  
  opt <- parse_args(OptionParser(option_list=option_list))

  
GTFtoGeneLength <- function (gtffile = GalaxyInputFile(required = T, formatFilter = "gft"), species = GalaxySelectParam(c("Human","Mouse","Rat","Pig","Cow","Horse","Zebrafish")),
                             geneExonLengths = GalaxyOutput("geneExonLengths", "tabular")) 
{
  
  suppressPackageStartupMessages(library("EnsDb.Hsapiens.v103"))
  suppressPackageStartupMessages(library("EnsDb.Mmusculus.v103"))
  suppressPackageStartupMessages(library("EnsDb.Rnorvegicus.v103"))
  suppressPackageStartupMessages(library("EnsDb.Ecaballus.v103"))
  suppressPackageStartupMessages(library("EnsDb.Btaurus.v103"))
  suppressPackageStartupMessages(library("EnsDb.Sscrofa.v103"))
  suppressPackageStartupMessages(library("EnsDb.Drerio.v103"))
  suppressPackageStartupMessages(library("GenomicFeatures"))

  
  txdb <- makeTxDbFromGFF(gtffile, format = "gtf", 
                          circ_seqs = character() )
  ebg <- exonsBy(txdb, by = "gene")
  exonic.gene.sizes <- lapply(ebg, function(x){sum(width(reduce(x)))})
  exonic.gene.sizes <-stack(exonic.gene.sizes)
  colnames(exonic.gene.sizes)<-c("length","gene_id")
  
  #convert IDs to gene symbol
  if (species=="Human"){
    endf <- GenomicFeatures::genes(EnsDb.Hsapiens.v103, return.type="DataFrame")
  } else if (species=="Rat"){
    endf <- GenomicFeatures::genes(EnsDb.Rnorvegicus.v103, return.type="DataFrame")
  } else if (species=="Cow"){
    endf <- GenomicFeatures::genes(EnsDb.Btaurus.v103, return.type="DataFrame")
  } else if (species=="Horse"){
    endf <- GenomicFeatures::genes(EnsDb.Ecaballus.v103, return.type="DataFrame")
  } else if (species=="Zebrafish"){
    endf <- GenomicFeatures::genes(EnsDb.Drerio.v103, return.type="DataFrame")
  } else if (species=="Pig"){
    endf <- GenomicFeatures::genes(EnsDb.Sscrofa.v103, return.type="DataFrame")
  } else {
    endf <- GenomicFeatures::genes(EnsDb.Mmusculus.v103, return.type="DataFrame")
  }
  en2gene <- na.omit(as.data.frame(endf[, c("gene_id", "gene_name")]))
  exonic.gene.sizes <- merge(exonic.gene.sizes, en2gene,by = "gene_id")
  exonic.gene.sizes <- exonic.gene.sizes[ ,c(3,2)]
  
  write.table(exonic.gene.sizes,file=geneExonLengths,col.names=T,row.names=F,sep="\t",quote=F)
  
}


opt <- opt[names(opt) != "help"]
do.call(GTFtoGeneLength, opt)
