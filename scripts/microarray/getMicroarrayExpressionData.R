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
option_list$splitPos <- make_option('--splitPos', type='character',help="position(s) of string to keep (single integer or comma-separated)")
option_list$newColumns <- make_option('--newColumns', type='character', default=NULL, help="Comma-separated names for new columns when splitting into multiple positions")
option_list$remove <- make_option('--remove', type='logical',help="If a sample should be excluded")
option_list$removeSample <- make_option('--removeSample', type='character',help="Identifier of sample to be removed")
option_list$site <- make_option('--site', type='character',help="Repository to download data from")
option_list$expressionData <- make_option('--expressionData', type='character',help="name of download data")


opt <- parse_args(OptionParser(option_list=option_list,description="Download raw microarray data from GEO or ArrayExpress"))

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

getArrayExpressSuppFiles <- function(accessionNumber, pData = FALSE) {
  
  # Determine experiment type
  type <- if(grepl("E-GEOD", accessionNumber)) {
    "E-GEOD-"
  } else if (grepl("E-MTAB", accessionNumber)) {
    "E-MTAB-"
  } else {
    "E-MEXP-"
  }
  
  digits <- str_sub(accessionNumber, -3, -1)
  
  # Build URL for fire directory
  base_url <- sprintf("https://ftp.ebi.ac.uk/biostudies/fire/%s/%s/%s/Files/",
                      type, digits, accessionNumber)
  print(paste("Accessing:", base_url))
  
  # Get directory listing via HTTPS
  response <- curl_fetch_memory(base_url)
  content <- rawToChar(response$content)
  
  # Extract filenames from HTML directory listing
  links <- regmatches(content, gregexpr('href="([^"]+)"', content))[[1]]
  file_list <- gsub('href="|"', '', links)
  # Filter out parent directory and query strings
  file_list <- file_list[!grepl("^\\?|^\\.\\./|^/", file_list)]
  
  if(length(file_list) == 0) {
    stop("No files found in directory")
  }
  
  # Filter files based on pData parameter
  if(pData) {
    file_list <- file_list[grepl("sdrf.txt", file_list)]
  } else {
    file_list <- file_list[!grepl("idf.txt|sdrf.txt", file_list)]
  }
  
  if(length(file_list) == 0) {
    stop("No matching files found")
  }
  
  # Create directory and download files
  dir.create(accessionNumber, showWarnings = FALSE)
  
  urls <- paste0(base_url, file_list)
  local_files <- paste0(accessionNumber, "/", basename(file_list))
  
  # Download files
  for(i in seq_along(urls)) {
    cat(sprintf("Downloading: %s\n", basename(urls[i])))
    curl_download(urls[i], local_files[i])
  }
  
  # Return pData if requested
  if(pData) {
    pdata <- read.delim(local_files[1], stringsAsFactors = FALSE)
    return(pdata)
  }
  
  invisible(local_files)
}


getMicroarrayExpressionData <- function (accessionNumber,
                                         platform = c("Affy","Affy-ST", "1C-Agilent", "2C-Agilent", "Illumina"),
                                         numberFileRemove=0, grepExpression =FALSE, grepString,
                                         numberLinesSkip = 0, split = FALSE, splitField , splitSep = " ", 
                                         splitPos = "1", newColumns = NULL, remove = FALSE, 
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
    tar_files <- grep(pattern = "\\.tar$", x = rownames(filePaths), value = TRUE)
    
    if (length(tar_files) > 0) {
      untar(tarfile = tar_files, exdir = accessionNumber)
      file.remove(tar_files)
    }
    
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
      
      files <- list.files( accessionNumber, pattern = paste0(accessionNumber, ".*\\.gz$"))
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
  
  # Handle splitting metadata field
  if (split == TRUE) {
    # Parse split positions - can be single or multiple positions
    splitPos_vec <- as.integer(trimws(strsplit(splitPos, ",")[[1]]))
    print(paste("Split positions:", paste(splitPos_vec, collapse = ", ")))
    
    # Get the target data frame (depends on platform)
    if (platform == "1C-Agilent") {
      target_df <- expData$targets
    } else {
      target_df <- pData(expData)
    }
    
    # Check if literal or not
    use_fixed <- !grepl("\\|", splitSep) 
    
    split_values <- strsplit(
      as.character(target_df[, splitField]),
      split = splitSep,
      fixed = use_fixed
    )
    
    
    # Check if we're creating multiple columns
    if (length(splitPos_vec) > 1) {
      # Multiple positions: create new columns
      if (is.null(newColumns) || newColumns == "") {
        stop("When splitting into multiple positions, newColumns must be provided")
      }
      
      # Parse new column names
      newColumn_names <- strsplit(newColumns, ",")[[1]]
      newColumn_names <- trimws(newColumn_names)
      
      # Validate that number of positions matches number of new column names
      if (length(splitPos_vec) != length(newColumn_names)) {
        stop("Number of split positions must match number of new column names")
      }
      
      print(paste("Creating", length(newColumn_names), "new columns:", 
                  paste(newColumn_names, collapse = ", ")))
      
      # Extract each position and create new columns
      for (i in seq_along(splitPos_vec)) {
        pos <- splitPos_vec[i]
        col_name <- newColumn_names[i]
        print(paste("Extracting position", pos, "into column", col_name))
        
        target_df[[col_name]] <- trimws(sapply(split_values, function(x) {
          if (length(x) >= pos) {
            return(x[pos])
          } else {
            return(NA_character_)
          }
        }))
        
        # Debug: show first few values
        print(paste("First values in", col_name, ":", 
                    paste(head(target_df[[col_name]], 3), collapse = ", ")))
      }
      
    } else {
      # Single position: replace the original field
      pos <- splitPos_vec[1]
      target_df[, splitField] <- trimws(sapply(split_values, function(x) {
        if (length(x) >= pos) x[pos] else NA
      }))
    }
    
    # Update the data back
    if (platform == "1C-Agilent") {
      expData$targets <- target_df
    } else {
      pData(expData) <- target_df
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
