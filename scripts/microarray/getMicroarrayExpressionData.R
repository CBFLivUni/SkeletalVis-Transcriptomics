#!/usr/bin/env Rscript

options("timeout" = 1200)

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("GEOquery"))
suppressPackageStartupMessages(library("limma"))
suppressPackageStartupMessages(library("oligo"))
suppressPackageStartupMessages(library("affy"))
suppressPackageStartupMessages(library("lumi"))
suppressPackageStartupMessages(library("curl"))
suppressPackageStartupMessages(library("stringr"))

.libPaths(".")

option_list <- list()

option_list$accessionNumber <- make_option('--accessionNumber', type='character',help="accession number dataset to download")
option_list$platform <- make_option('--platform', type='character',help="The type of microarray platform")
option_list$numberFileRemove <- make_option('--numberFileRemove', type='integer',default=0,help="The number of files to remove")
option_list$grepExpression <- make_option('--grepExpression', type='logical',help="Should illumina data processed by searching for a string")
option_list$grepString <- make_option('--grepString', type='character',help="string to search for with illumina data")
option_list$numberLinesSkip <- make_option('--numberLinesSkip', type='integer',default=0,help="number of lines to skip for illumina data")
option_list$split <- make_option('--split', type='logical', default=FALSE,help="If a metadata field by split")
option_list$splitField <- make_option('--splitField', type='character',help="Metadata field to be split")
option_list$splitSep <- make_option('--splitSep', type='character',help="seperator to split metdata field on")
option_list$splitPos <- make_option('--splitPos', type='integer',help="position of string to keep")
option_list$remove <- make_option('--remove', type='logical',help="If a sample should be excluded")
option_list$removeSample <- make_option('--removeSample', type='character',help="Identifier of sample to be removed")
option_list$site <- make_option('--site', type='character',help="Repository to download data from")
option_list$expressionData <- make_option('--expressionData', type='character',help="name of download data")


opt <- parse_args(OptionParser(option_list=option_list,description="Download raw microarray data from GEO or ArrayExpress"))


# getArrayExpressSuppFiles <- function(accessionNumber) {
#   baseurl <- "https://www.ebi.ac.uk/arrayexpress/files/"
#   endurl <- ".raw.1.zip"
#   downloadUrl <- paste(baseurl, accessionNumber, "/", accessionNumber, 
#                        endurl, sep = "")
#   destFile <- paste(accessionNumber, ".zip", sep = "")
#   attempt <- 1
#   while ((is.na(file.size(destFile)) || file.size(destFile) < 
#           1000) && attempt <= 5) {
#     if (attempt > 1) {
#       Sys.sleep(60)
#     }
#     attempt <- attempt + 1
#     try(download.file(downloadUrl, destFile))
#   }
#   if (is.na(file.size(destFile)) || file.size(destFile) < 
#       1000) {
#     stop(sprintf("Failed to access raw data file for accession number %s", 
#                  accessionNumber))
#   }
#   return(destFile)
# }

# getArrayExpressSuppFiles <- function(accessionNumber) {
#   baseurl <- "https://www.ebi.ac.uk/arrayexpress/files/"
#   endurl <- ".sdrf.txt"
#   downloadUrl <- paste(baseurl, accessionNumber, "/", accessionNumber, 
#                        endurl, sep = "")
#   destFile <- paste(accessionNumber, ".txt", sep = "")
#   attempt <- 1
#   while ((is.na(file.size(destFile)) || file.size(destFile) < 
#           1000) && attempt <= 5) {
#     if (attempt > 1) {
#       Sys.sleep(60)
#     }
#     attempt <- attempt + 1
#     try(download.file(downloadUrl, destFile))
#   }
#   if (is.na(file.size(destFile)) || file.size(destFile) < 
#       1000) {
#     stop(sprintf("Failed to access phenotype data for accession number %s", 
#                  accessionNumber))
#   }
#   pdata <- read.delim(destFile, stringsAsFactors = F)
#   return(pdata)
# }
getGEOPhenotypeData <- function(accessionNumber,platform) {
  eset <- try(getGEO(accessionNumber)[[1]])
  attempt <- 2
  while (inherits(eset, "try-error") && attempt <= 5) {
    Sys.sleep(60)
    attempt <- attempt + 1
    eset <- try(getGEO(accessionNumber)[[1]])
  }
  if (inherits(eset, "try-error")) {
    stop(sprintf("Failed to access GEO phenotype data for accession number %s", 
                 accessionNumber))
  }
  if (platform != "2C-Agilent") {
    columnNames <- sapply(strsplit(x = as.character(colnames(pData(eset))), 
                                   split = ":"), "[[", 1)
    colnames(pData(eset)) <- make.names(columnNames, 
                                        unique = T)
    pData(eset) <- as.data.frame(apply(pData(eset), 2, 
                                       trimws))
  }
  return(eset)
}

getArrayExpressSuppFiles <- function(accessionNumber,pData=FALSE) {
  
  if(grepl("E-GEOD",accessionNumber)){
    type <- "E-GEOD-"
  } else if (grepl("E-MTAB",accessionNumber)){
    type <- "E-MTAB-"
  } else{
    type <- "E-MEXP-"
  }
  
  digits <- str_sub(accessionNumber,-3,-1)
  
  url <- sprintf("ftp://ftp.ebi.ac.uk/biostudies/nfs/%s/%s/%s/Files/",type,digits,accessionNumber)
  h <- new_handle(dirlistonly=TRUE)
  con <- try(curl(url, "r", h))
  print(url)
  if(inherits(con,"try-error")){
    url <- sprintf("ftp://ftp.ebi.ac.uk/biostudies/fire/%s/%s/%s/Files/",type,digits,accessionNumber)
    h <- new_handle(dirlistonly=TRUE)
    con <- curl(url, "r", h)
  }
  tbl <- read.table(con, stringsAsFactors=TRUE, fill=TRUE)
  
  
  if(pData) {
  tbl <- tbl[ grepl("sdrf.txt",tbl[,1]),]
  } else{
  tbl <- tbl[ !grepl("idf.txt|sdrf.txt",tbl[,1]),]
  }
  
  dir.create(accessionNumber)
  urls <- paste0(url, tbl)
  fls <- paste0(accessionNumber,"/",basename(urls))
  mapply(curl_download,urls, fls)
  
  if(pData){
    pdata <- read.delim(fls, stringsAsFactors = F)
    return(pdata)
  }
}

getMicroarrayExpressionData <- function (accessionNumber,
                                         platform = c("Affy","Affy-ST", "1C-Agilent", "2C-Agilent", "Illumina"),
                                          numberFileRemove=0, grepExpression =FALSE, grepString,
    numberLinesSkip = 0, split = FALSE, splitField , splitSep = " ", 
    splitPos = 1, remove = FALSE, 
    removeSample, site = c("GEO","ArrayExpress")
    ,expressionData ="expressionData.RDS") 
{
  
  #remove the added marker on the accession if needed
    accessionNumber <- gsub("_B$","",accessionNumber)
    
    if (site == "GEO") {
        filePaths <- try(getGEOSuppFiles(accessionNumber))
        attempt <- 2
        while (inherits(filePaths, "try-error") && attempt <= 
            5) {
            Sys.sleep(60)
            attempt <- attempt + 1
            filePaths <- try(getGEOSuppFiles(accessionNumber))
        }
        if (inherits(filePaths, "try-error")) {
            stop(sprintf("Failed to access GEO raw data for accession number %s", 
                accessionNumber))
        }
        untar(tarfile = grep(pattern = ".tar", x = rownames(filePaths), 
            value = T), exdir = accessionNumber)
        file.remove(grep(pattern = ".tar", x = rownames(filePaths), 
            value = T))
    }   else {
        getArrayExpressSuppFiles(accessionNumber)

    }
    if (platform == "1C-Agilent") {
        if (site == "GEO") {
            files <- list.files(accessionNumber, pattern = ".txt")
            if (numberFileRemove > 0) {
                files <- files[-1:-numberFileRemove]
            }
            setwd(accessionNumber)
            sapply(files, gunzip)
            files <- list.files(pattern = ".txt")
            if (numberFileRemove > 0) {
                files <- files[-1:-numberFileRemove]
            }
            expData <- read.maimages(path = ".", files = files, 
                source = "agilent", green.only = T)
            eset <- getGEOPhenotypeData(accessionNumber, platform)
            targets <- pData(eset)
            expData$targets <- targets
            setwd("..")
        } else {
            setwd(accessionNumber)
            files <- list.files(pattern = ".txt")
            if (numberFileRemove > 0) {
                files <- files[-1:-numberFileRemove]
            }
            print(files)
            expData <- read.maimages(path = ".", files = files, 
                source = "agilent", green.only = T)
            pdata <- getArrayExpressSuppFiles(accessionNumber,pData=TRUE)
            expData <- expData[, sort(colnames(expData$E))]
            pdata <- pdata[order(as.character(pdata$Array.Data.File)), 
                ]
            expData$targets <- pdata
            setwd("..")
        }
    } else if (platform == "2C-Agilent") {
        setwd(accessionNumber)
        files <- list.files(pattern = ".txt")
        if (numberFileRemove > 0) {
            files <- files[-1:-numberFileRemove]
        }
        expData <- read.maimages(path = ".", files = files, source = "agilent", 
            green.only = F)
        if (site == "GEO") {
            eset <- getGEOPhenotypeData(accessionNumber, platform)
            targets <- pData(eset)
            expData$targets <- targets
        }  else {
            pdata <- getArrayExpressSuppFiles(accessionNumber,pData=TRUE)
            expData$targets <- pdata
            n <- nrow(expData$targets)
            expData$targets <- merge(expData$targets[1:(n/2), 
                ], expData$targets[(n/2) + 1:n, ], by = "Hybridization.Name")
        }
        setwd("..")
    }  else if (platform == "Affy-ST") {
        files <- list.files(accessionNumber, pattern = ".CEL", 
            ignore.case = T)
        setwd(accessionNumber)
        expData <- oligo::read.celfiles(files)
        if (site == "GEO") {
            eset <- getGEOPhenotypeData(accessionNumber, platform)
            rownames(pData(eset)) <- rownames(pData(expData))
            pData(expData) <- pData(eset)
        }  else {
            pdata <- getArrayExpressSuppFiles(accessionNumber,pData=TRUE)
            expData <- expData[, sort(sampleNames(expData))]
            pdata <- pdata[order(as.character(pdata$Array.Data.File)), ]
            rownames(pdata) <- rownames(pData(expData))
            pData(expData) <- pdata
        }
        setwd("..")
    } else if (platform == "Illumina") {
        if (site == "GEO") {
            files <- list.files(accessionNumber, pattern = accessionNumber)
            setwd(accessionNumber)
            sapply(files, gunzip)
            file <- list.files(pattern = accessionNumber)
            if (length(file) > 1) {
            		if(any(grepl("non-normalized",file))){
            		file <- file[grep("non-normalized",file)]
            		}
              file <- file[length(file)]
            }
            data <- read.delim(file, skip = numberLinesSkip)
            if (grepExpression == T) {
                exprs <- data[, grep(grepString, colnames(data))]
            }  else {
                exprs <- data[, seq(from = 2, to = ncol(data), by = 2)]
                exprs <- exprs[, colSums(is.na(exprs)) < nrow(exprs)]
            }
            rownames(exprs) <- data[, 1]
            eset <- getGEOPhenotypeData(accessionNumber, platform)
            rownames(pData(eset)) <- colnames(exprs)
            expData <- ExpressionSet(assayData = as.matrix(exprs), 
                phenoData = phenoData(eset))
            setwd("..")
        } else {
            file <- list.files(accessionNumber, pattern = ".txt", 
                full.names = T)
            expData <- lumiR(file)
            pdata <- getArrayExpressSuppFiles(accessionNumber,pData=TRUE)
            expData <- expData[, sort(sampleNames(expData))]
            pdata <- pdata[order(as.character(pdata$Array.Data.File)), ]
            pData(expData) <- pdata
            rownames(pData(expData)) <- colnames(exprs(expData))
            expData <- ExpressionSet(assayData = as.matrix(exprs(expData)), 
                phenoData = phenoData(expData))
        }
    } else if (platform == "Affy") {
        files <- list.files(accessionNumber, pattern = ".CEL", 
            ignore.case = T)
        setwd(accessionNumber)
        expData <- ReadAffy(filenames = files)
        if (site == "GEO") {
            eset <- getGEOPhenotypeData(accessionNumber, platform)
            rownames(pData(eset)) <- rownames(pData(expData))
            pData(expData) <- pData(eset)
        }   else {
            pdata <- getArrayExpressSuppFiles(accessionNumber,pData=TRUE)
            expData <- expData[, sort(sampleNames(expData))]
            pdata <- pdata[order(as.character(pdata$Array.Data.File)), ]
            rownames(pdata) <- rownames(pData(expData))
            pData(expData) <- pdata
        }
        setwd("..")
    }
    if (split == TRUE) {
        if (platform == "1C-Agilent") {
            expData$targets[, splitField] <- trimws(sapply(strsplit(as.character(expData$targets[, 
                splitField]), split = splitSep), "[[", splitPos))
        }   else {
            pData(expData)[splitField] <- trimws(sapply(strsplit(as.character(pData(expData)[, 
                splitField]), split = splitSep), "[[", splitPos))
        }
    }
    if (remove == TRUE) {
        removeSample <- strsplit(removeSample, "_AND_", fixed = F)[[1]]
        expData <- expData[, !sampleNames(expData) %in% removeSample]
    }
    saveRDS(expData, file = expressionData)
}

opt <- opt[names(opt)!="help"]
do.call(getMicroarrayExpressionData, opt)