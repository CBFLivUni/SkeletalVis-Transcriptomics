suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("curl"))
suppressPackageStartupMessages(library("tidyverse"))

option_list <- list()

option_list$accessionNumber <- make_option('--accessionNumber', type='character', help="Project study accession number e.g PRJEB12345")
option_list$sampleTable <- make_option('--sampleTable', type='character', help="path to sample table with a File column containing SRR accessions")
option_list$aspera <- make_option('--aspera', type='logical',default=FALSE, help="use aspera to download?")
option_list$key <- make_option('--key', type='character', help="absolute path to aspera key")

opt <- parse_args(OptionParser(option_list=option_list,description = "download fastq files from ENA given project and sample accession numbers"))


downloadFastqFile <- function(fastq,md5){
  
  fileName <- basename(fastq)
  print(sprintf("downloading %s",fileName))
  curl_fetch_disk(fastq,fileName)
  #check the md5sum is correct - if not try downloading again
  if(tools::md5sum(fileName)!=md5) {
    print(sprintf("retrying %s",fileName))
    curl_fetch_disk(fastq,fileName)
    #if still wrong inform the user
    if(tools::md5sum(fileName)!=md5) {
      stop(sprintf("Error: MD5 sum of %s incorrect",fileName))
    }
  }
}


downloadAspera <- function(fastq,md5,key){
  
  fileName <- basename(fastq)
  print(sprintf("downloading %s",fileName))
  
  cmd <- sprintf("ascp -QT -l 300m -P 33001 -i %s era-fasp@%s %s",key,fastq,fileName)
  print(cmd)
  system(cmd)
  
  #check the md5sum is correct - if not try downloading again
  if(tools::md5sum(fileName)!=md5) {
    print(sprintf("retrying %s",fileName))
    system(cmd)
    #if still wrong inform the user
    if(tools::md5sum(fileName)!=md5) {
      stop(sprintf("Error: MD5 sum of %s incorrect",fileName))
    }
  }
}

getRNASeqExpressionData <- function (accessionNumber, sampleTable,key,aspera) {
  
  #read in the SRR ids
  sampleTable <- read.delim(sampleTable)
  
  #split the accession number in case the fastq files have been split into sub-series
  accessionNumber <- strsplit(accessionNumber,split = "|",fixed = TRUE)[[1]]
  
  for (i in accessionNumber) {
    
    #download the table of download urls using the accession number
    url <- paste0("http://www.ebi.ac.uk/ena/portal/api/filereport?accession=",i,
                  "&result=read_run&fields=study_accession,run_accession,fastq_ftp,fastq_aspera,fastq_md5&format=tsv")
    
    results <- read.delim(url,stringsAsFactors = FALSE)
    
    
    if(aspera){
      
      results <- results[results$fastq_aspera!= "", ] %>% separate_rows(fastq_aspera,fastq_md5,sep=";") %>% as.data.frame()
      results <- results[grep(paste(sampleTable$File, collapse = "|"),results$run_accession),]
      
      mapply(downloadAspera,results$fastq_aspera,results$fastq_md5,MoreArgs = list(key=key))
      
      
    } else{
      
      results <- results[results$fastq_ftp!= "", ] %>% separate_rows(fastq_ftp,fastq_md5,sep=";") %>% as.data.frame()
      results <- results[grep(paste(sampleTable$File, collapse = "|"),results$run_accession),]
      results$fastq_ftp <- paste0("ftp://", results$fastq_ftp)
      
      mapply(downloadFastqFile,results$fastq_ftp,results$fastq_md5)
    }
  }
  
}


opt <- opt[names(opt) != "help"]
do.call(getRNASeqExpressionData, opt)


