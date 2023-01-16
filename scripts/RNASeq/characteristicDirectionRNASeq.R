#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("genefilter"))
suppressPackageStartupMessages(library("EnsDb.Hsapiens.v103"))
suppressPackageStartupMessages(library("EnsDb.Mmusculus.v103"))
suppressPackageStartupMessages(library("EnsDb.Rnorvegicus.v103"))
suppressPackageStartupMessages(library("EnsDb.Ecaballus.v103"))
suppressPackageStartupMessages(library("EnsDb.Btaurus.v103"))
suppressPackageStartupMessages(library("EnsDb.Sscrofa.v103"))
suppressPackageStartupMessages(library("EnsDb.Drerio.v103"))
suppressPackageStartupMessages(library("sva"))
suppressPackageStartupMessages(library("GeoDE"))

option_list <- list()

option_list$txiData <- make_option('--txiData', type='character')
option_list$sampleTable <- make_option('--sampleTable', type='character')
option_list$comparisonsTable <- make_option('--comparisonsTable', type='character')
option_list$species <- make_option('--species', type='character')
option_list$foldChangeOnly <- make_option('--foldChangeOnly', type='logical')
option_list$technicalReplicates <- make_option('--technicalReplicates', type='logical')
option_list$batchCorrect <- make_option('--batchCorrect', type='logical')
option_list$supervised <- make_option('--supervised', type='logical')
option_list$chrDirTable <- make_option('--chrDirTable', type='character')


opt <- parse_args(OptionParser(option_list=option_list))

getResultsDataFrame <- function(dds, comparisonNum, comparisonsTable, correctedMatrix) {
  
  contrast <- comparisonsTable[comparisonNum, 1]
  numerator <- comparisonsTable[comparisonNum, 2]
  denominator <- comparisonsTable[comparisonNum, 3]
  keep <- which(colData(dds)[[contrast]] %in% c(numerator,denominator))
  
  if (is.null(correctedMatrix)) {
    exprs <- as.matrix(assay(varianceStabilizingTransformation(dds)))
  }  else {
    exprs <- correctedMatrix
  }
  
  exprs <- exprs[, keep]
  dds <- dds[, keep]
  classes <- as.factor(ifelse(colData(dds)[[contrast]] %in% numerator, 1, 2))
  exprs <- as.data.frame(exprs)
  exprs <- na.omit(cbind(rownames(exprs), exprs))
  chdir <- chdirAnalysis(exprs, classes, 1, CalculateSig = F)
  results <- stack(chdir$chdirprops$chdir[[1]][,1])[, 2:1]
  results <- results[order(abs(results[,2]), decreasing = T), ]
  #results <- results[1:500, ]
  results$comparison <- paste(numerator, denominator, sep = "vs")
  results$comparisonNumber <- comparisonNum
  colnames(results)[1:2] <- c("ID", "chrDir")
  return(results)
}

regressSVs <- function(y, mod, svaobj, P = ncol(mod)) {
  
  X <- cbind(mod, svaobj$sv)
  Hat <- solve(t(X) %*% X) %*% t(X)
  beta <- (Hat %*% t(y))
  cleany <- y - t(as.matrix(X[, -c(1:P)]) %*% beta[-c(1:P),])
  return(cleany)
}


characteristicDirectionRNASeq <- function (txiData, 
    sampleTable, 
    comparisonsTable, 
    species = c("Human", "Mouse", "Rat", "Pig", 
        "Cow", "Horse", "Zebrafish"), foldChangeOnly =FALSE, 
    technicalReplicates = FALSE, batchCorrect = FALSE, 
    supervised = FALSE, chrDirTable = "chrDirTable.txt") {

 
    if (foldChangeOnly == TRUE) {
        return("need replicates to perfrom characteristic direction analysis")
    }
    
    txi <- readRDS(txiData)
    sampleTable <- read.delim(sampleTable)
    comparisonsTable <- read.delim(comparisonsTable, stringsAsFactors = F)
    sampleTable <- sampleTable[, !grepl("File",colnames(sampleTable))]
    txi$counts <- txi$counts[, match(as.character(sampleTable[, 1]), colnames(txi$counts))]
    txi$abundance <- txi$abundance[, match(as.character(sampleTable[, 1]), colnames(txi$abundance))]
    txi$length <- txi$length[, match(as.character(sampleTable[, 1]), colnames(txi$length))]
    factors <- colnames(sampleTable)[-1]
    designFormula <- as.formula(paste("~", paste(factors, collapse = " + ")))
    
    if (technicalReplicates) {
        dds <- DESeqDataSetFromTximport(txi = txi, colData = sampleTable, 
            design = designFormula, ignoreRank = T)
        dds <- collapseReplicates(dds, paste0(sampleTable[, factors[1]], 
            sampleTable[, factors[2]]))
        design(dds) <- as.formula(paste("~", paste(factors[-1], 
            collapse = " + ")))
    }   else {
        dds <- DESeqDataSetFromTximport(txi = txi, colData = sampleTable, 
            design = designFormula)
    }
    
    if (batchCorrect) {
        idx <- rowMeans(counts(dds)) > 1
        dat <- counts(dds)[idx, ]
        mod <- model.matrix(designFormula, colData(dds))
        mod0 <- model.matrix(~1, colData(dds))
        if (supervised) {
            dds <- DESeq(dds[idx, ])
            resultsTable <- as.data.frame(results(dds, contrast = c(comparisonsTable[1, 1], comparisonsTable[1, 2],
                                                                    comparisonsTable[1,3])))
            resultsTable <- resultsTable[order(abs(resultsTable$pvalue)),]
            controls <- !rownames(dds) %in% rownames(resultsTable)[1:5000]
            svseq <- svaseq(dat, mod, mod0, controls = controls)
        }   else {
            svseq <- svaseq(dat, mod, mod0)
        }
        
        correctedMatrix <- regressSVs(dat, mod, svseq)
        correctedMatrix <- correctedMatrix + abs(min(correctedMatrix))
        dds <- DESeqDataSetFromMatrix(round(correctedMatrix), 
            colData = sampleTable, design = designFormula)
        correctedMatrix <- as.matrix(assay(varianceStabilizingTransformation(dds)))
        resultsTable <- lapply(1:nrow(comparisonsTable), function(x) getResultsDataFrame(dds, 
            comparisonNum = x, comparisonsTable = comparisonsTable, 
            correctedMatrix = correctedMatrix))
    }  else {
        resultsTable <- lapply(1:nrow(comparisonsTable), function(x) getResultsDataFrame(dds, 
            comparisonNum = x, comparisonsTable = comparisonsTable, 
            correctedMatrix = NULL))
    }
    
    resultsTable <- do.call(rbind, resultsTable)
    
    if (species == "Human") {
        endf <- GenomicFeatures::genes(EnsDb.Hsapiens.v103, return.type = "DataFrame")
    }    else if (species == "Rat") {
        endf <- GenomicFeatures::genes(EnsDb.Rnorvegicus.v103, 
            return.type = "DataFrame")
    }    else if (species == "Cow") {
        endf <- GenomicFeatures::genes(EnsDb.Btaurus.v103, return.type = "DataFrame")
    }    else if (species == "Horse") {
        endf <- GenomicFeatures::genes(EnsDb.Ecaballus.v103, return.type = "DataFrame")
    }    else if (species == "Zebrafish") {
        endf <- GenomicFeatures::genes(EnsDb.Drerio.v103, return.type = "DataFrame")
    }    else if (species == "Pig") {
        endf <- GenomicFeatures::genes(EnsDb.Sscrofa.v103, return.type = "DataFrame")
    }    else {
        endf <- GenomicFeatures::genes(EnsDb.Mmusculus.v103, return.type = "DataFrame")
    }
    
    en2gene <- as.data.frame(endf[, c("gene_id", "gene_name")])
    resultsTable <- merge(resultsTable, en2gene, by.x = "ID", 
        by.y = "gene_id")
    resultsTable <- resultsTable[order(resultsTable$comparisonNumber, 
        rev(abs(resultsTable$chrDir))), ]
    
    write.table(resultsTable, file = chrDirTable, col.names = TRUE, 
        row.names = FALSE, sep = "\t", quote = FALSE)
}

opt <- opt[names(opt) != "help"]
do.call(characteristicDirectionRNASeq, opt)
