suppressPackageStartupMessages(library("EnsDb.Hsapiens.v103"))
suppressPackageStartupMessages(library("EnsDb.Mmusculus.v103"))
suppressPackageStartupMessages(library("EnsDb.Rnorvegicus.v103"))
suppressPackageStartupMessages(library("EnsDb.Ecaballus.v103"))
suppressPackageStartupMessages(library("EnsDb.Btaurus.v103"))
suppressPackageStartupMessages(library("EnsDb.Sscrofa.v103"))
suppressPackageStartupMessages(library("EnsDb.Drerio.v103"))

suppressPackageStartupMessages(library("tximport"))


KallistoAbundancestoGeneCountMatrix<-function(filepaths,filenames,species=GalaxySelectParam(c("Human","Mouse","Rat","Pig","Cow","Horse","Zebrafish")),output,outputTPM){
  
  #extract the ensembl transcripts
  if (species=="Human"){
    txdf <- transcripts(EnsDb.Hsapiens.v103, return.type="DataFrame")
  } else if (species =="Cow") {
    txdf <- transcripts(EnsDb.Btaurus.v103, return.type="DataFrame")
  } else if (species =="Zebrafish") {
    txdf <- transcripts(EnsDb.Drerio.v103, return.type="DataFrame")
  } else if (species =="Pig") {
    txdf <- transcripts(EnsDb.Sscrofa.v103, return.type="DataFrame")
  }  else if (species =="Horse") {
    txdf <- transcripts(EnsDb.Ecaballus.v103, return.type="DataFrame")
  }  else if (species =="Rat") {
    txdf <- transcripts(EnsDb.Rnorvegicus.v103, return.type="DataFrame")
  }  else {
	  txdf <- transcripts(EnsDb.Mmusculus.v103, return.type="DataFrame")
  }

  #make a dataframe as required for tximport
  tx2gene <- as.data.frame(txdf[,c("tx_id","gene_id")])
  
  #summarise the transcript level counts into gene level counts
  filenames<-gsub(pattern = ".abundances",x = filenames,replacement = "")
  names(filepaths)<-filenames
  txi <- tximport(filepaths, type = "kallisto", tx2gene = tx2gene,ignoreTxVersion=T)
  txiScaledTPM <- tximport(filepaths, type = "kallisto", tx2gene = tx2gene,ignoreTxVersion=T,countsFromAbundance = "lengthScaledTPM")
  saveRDS(txi,file=output)
  saveRDS(txiScaledTPM,file=outputTPM)

}


#need to parse the data collection list of file paths which is of unknown length
args <- commandArgs(trailingOnly = TRUE)

hh <- paste(unlist(args),collapse=' ')
listoptions <- unlist(strsplit(hh,'--'))[-1]

options.args <- sapply(listoptions,function(x){
  unlist(strsplit(x, ' '))[-1]
})
options.names <- sapply(listoptions,function(x){
  option <-  unlist(strsplit(x, ' '))[1]
})
names(options.args) <- unlist(options.names)

KallistoAbundancestoGeneCountMatrix(options.args$filepaths,options.args$filenames,options.args$species,options.args$output,options.args$outputTPM)

