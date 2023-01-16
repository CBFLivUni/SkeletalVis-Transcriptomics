#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("genefilter"))
suppressPackageStartupMessages(library("EnsDb.Hsapiens.v103"))
suppressPackageStartupMessages(library("EnsDb.Mmusculus.v103"))
suppressPackageStartupMessages(library("EnsDb.Rnorvegicus.v103"))
suppressPackageStartupMessages(library("EnsDb.Ecaballus.v103"))
suppressPackageStartupMessages(library("EnsDb.Btaurus.v103"))
suppressPackageStartupMessages(library("EnsDb.Sscrofa.v103"))
suppressPackageStartupMessages(library("EnsDb.Drerio.v103"))
suppressPackageStartupMessages(library("sva"))

option_list <- list()

option_list$txiData <- make_option('--txiData', type='character',help="path to the tximport output")
option_list$sampleTable <- make_option('--sampleTable', type='character',help="path to sample table detailing phenotype info for each sample")
option_list$comparisonsTable <- make_option('--comparisonsTable', type='character',help="path to comparisons table with contrasts to make")
option_list$species <- make_option('--species', type='character',help="species under study")
option_list$foldChangeOnly <- make_option('--foldChangeOnly', type='logical',help="does the experiment have replicates")
option_list$technicalReplicates <- make_option('--technicalReplicates', type='logical',help="are there technical replicates to merge")
option_list$batchCorrect <- make_option('--batchCorrect', type='logical',help="should batch correction be used")
option_list$supervised <- make_option('--supervised', type='logical',help="should supervision svaseq batch correct be used")
option_list$foldChangeTable <- make_option('--foldChangeTable', type='character',help="path to output fold change table of results")

opt <- parse_args(OptionParser(option_list=option_list,description = "differential expression analysis with DESeq2"))


getFoldChangeDataFrame <- function(rld, sampleTable, condition, 
                                   numerator, denominator) {
    
    numerator.ind <- which(sampleTable[, condition] %in% 
                               numerator)
    denominator.ind <- which(sampleTable[, condition] %in% 
                                 denominator)
    rld <- assay(rld)
    if (length(numerator.ind) > 1) {
        numerator.dat <- rowMeans(rld[, numerator.ind])
    }       else {
        numerator.dat <- rld[, numerator.ind]
    }
    if (length(denominator.ind) > 1) {
        denominator.dat <- rowMeans(rld[, denominator.ind])
    }        else {
        denominator.dat <- rld[, denominator.ind]
    }
    data <- data.frame(log2FC = numerator.dat - denominator.dat)
    colnames(data) <- paste(numerator, denominator, sep = "vs")
    return(data)
}

getResultsDataFrame <- function(dds, condition, numerator, 
                                denominator) {
    data <- as.data.frame(lfcShrink(dds, contrast = c(condition, 
                                                      numerator, denominator),parallel=TRUE,BPPARAM = MulticoreParam(4),type="ashr"))
    data <- data[, c("log2FoldChange", "padj")]
    colnames(data) <- paste(paste(numerator, denominator, sep = "vs"), colnames(data), sep = "_")
    return(data)
}

DESeq2FoldChange <- function (txiData, 
    sampleTable, 
    comparisonsTable, 
    species, foldChangeOnly=TRUE, 
    technicalReplicates=FALSE , 
    batchCorrect = FALSE, supervised = FALSE, 
    foldChangeTable = "foldChangeTable.txt") {

    txi <- readRDS(txiData)
    sampleTable <- read.delim(sampleTable)
    sampleTable <- sampleTable[, !grepl("File",colnames(sampleTable))]
    comparisonsTable <- read.delim(comparisonsTable, stringsAsFactors = F)
    txi$counts <- txi$counts[, match(as.character(sampleTable[, 1]), colnames(txi$counts))]
    txi$abundance <- txi$abundance[, match(as.character(sampleTable[,1]), colnames(txi$abundance))]
    txi$length <- txi$length[, match(as.character(sampleTable[, 1]), colnames(txi$length))]
    factors <- colnames(sampleTable)[-1]
    
    designFormula <- as.formula(paste("~", paste(factors, 
            collapse = " + ")))
    
    if (technicalReplicates) {
        dds <- DESeqDataSetFromTximport(txi = txi, colData = sampleTable, 
            design = designFormula, ignoreRank = T)
        dds <- collapseReplicates(dds, paste0(sampleTable[, factors[1]], 
            sampleTable[, factors[2]]))
        design(dds) <- as.formula(paste("~", paste(factors[-1], collapse = " + ")))
    }    else {
        dds <- DESeqDataSetFromTximport(txi = txi, colData = sampleTable, 
            design = designFormula)
    }
    
    if (foldChangeOnly == TRUE) {
        rld <- rlogTransformation(dds)
        resultsTable <- lapply(1:nrow(comparisonsTable), function(x) getFoldChangeDataFrame(rld, 
            sampleTable, comparisonsTable[x, 1], comparisonsTable[x,2], comparisonsTable[x, 3]))
    }    else {
         if (batchCorrect) {
                idx <- rowMeans(counts(dds)) > 1
                dat <- counts(dds)[idx, ]
                mod <- model.matrix(designFormula, colData(dds))
                mod0 <- model.matrix(~1, colData(dds))
                if (supervised) {
                  dds <- DESeq(dds[idx, ])
                  resultsTable <- as.data.frame(results(dds, 
                    contrast = c(comparisonsTable[1, 1], comparisonsTable[1,2], comparisonsTable[1, 3])))
                  resultsTable <- resultsTable[order(abs(resultsTable$pvalue)), ]
                  controls <- !rownames(dds) %in% rownames(resultsTable)[1:5000]
                  svseq <- svaseq(dat, mod, mod0, controls = controls)
                }
                else {
                  svseq <- svaseq(dat, mod, mod0)
                }
                colnames(svseq$sv) <- paste0("batch", 1:svseq$n.sv)
                sampleTable <- cbind(svseq$sv, sampleTable)
                designFormula <- as.formula(paste("~", colnames(svseq$sv), 
                  "+", paste(factors, collapse = " + ")))
                dds <- DESeqDataSetFromTximport(txi = txi, colData = sampleTable, 
                  design = designFormula)
            }
            dds <- DESeq(dds,parallel=T,BPPARAM = MulticoreParam(4))
            resultsTable <- lapply(1:nrow(comparisonsTable), 
                function(x) getResultsDataFrame(dds, comparisonsTable[x, 
                  1], comparisonsTable[x, 2], comparisonsTable[x, 
                  3]))
        
    }
    resultsTable <- do.call(cbind, resultsTable)
    if (species == "Human") {
        endf <- GenomicFeatures::genes(EnsDb.Hsapiens.v103, return.type = "DataFrame")
    }
    else if (species == "Rat") {
        endf <- GenomicFeatures::genes(EnsDb.Rnorvegicus.v103, 
            return.type = "DataFrame")
    }
    else if (species == "Cow") {
        endf <- GenomicFeatures::genes(EnsDb.Btaurus.v103, return.type = "DataFrame")
    }
    else if (species == "Horse") {
        endf <- GenomicFeatures::genes(EnsDb.Ecaballus.v103, return.type = "DataFrame")
    }
    else if (species == "Zebrafish") {
        endf <- GenomicFeatures::genes(EnsDb.Drerio.v103, return.type = "DataFrame")
    }
    else if (species == "Pig") {
        endf <- GenomicFeatures::genes(EnsDb.Sscrofa.v103, return.type = "DataFrame")
    }
    else {
        endf <- GenomicFeatures::genes(EnsDb.Mmusculus.v103, return.type = "DataFrame")
    }
    en2gene <- as.data.frame(endf[, c("gene_id", "gene_name")])
    resultsTable <- merge(resultsTable, en2gene, by.x = "row.names", 
        by.y = "gene_id")
    resultsTable <- resultsTable[, c(ncol(resultsTable), 3:ncol(resultsTable) - 
        1, 1)]
    colnames(resultsTable)[ncol(resultsTable)] <- "ID"
    write.table(resultsTable, file = foldChangeTable, col.names = TRUE, 
        row.names = FALSE, sep = "\t", quote = FALSE)
}


opt <- opt[names(opt) != "help"]
do.call(DESeq2FoldChange, opt)

