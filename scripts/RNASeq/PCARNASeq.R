#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("DESeq2"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("cowplot"))
suppressPackageStartupMessages(library("genefilter"))
suppressPackageStartupMessages(library("EnsDb.Hsapiens.v103"))
suppressPackageStartupMessages(library("EnsDb.Mmusculus.v103"))
suppressPackageStartupMessages(library("EnsDb.Rnorvegicus.v103"))
suppressPackageStartupMessages(library("EnsDb.Btaurus.v103"))
suppressPackageStartupMessages(library("EnsDb.Sscrofa.v103"))
suppressPackageStartupMessages(library("EnsDb.Drerio.v103"))
suppressPackageStartupMessages(library("sva"))

option_list <- list()

option_list$txiData <- make_option('--txiData', type='character')
option_list$sampleTable <- make_option('--sampleTable', type='character')
option_list$comparisonsTable <- make_option('--comparisonsTable', type='character')
option_list$species <- make_option('--species', type='character')
option_list$technicalReplicates <- make_option('--technicalReplicates', type='logical')
option_list$batchCorrect <- make_option('--batchCorrect', type='logical')
option_list$supervised <- make_option('--supervised', type='logical')
option_list$pcaPlot <- make_option('--pcaPlot', type='character')
option_list$geneInfluence <- make_option('--geneInfluence', type='character')
option_list$normalisedExp <- make_option('--normalisedExp', type='character')

opt <- parse_args(OptionParser(option_list=option_list))

PCA <- function(object, intgroup, ntop = 30000, correctedMatrix) {
  if (is.null(correctedMatrix)) {
    rv <- rowVars(assay(object))
    select <- order(rv, decreasing = TRUE)[seq_len(min(ntop,length(rv)))]
    pca <- prcomp(t(assay(object)[select, ]))
  }  else {
    pca <- prcomp(t(correctedMatrix))
  }
  percentVar <- pca$sdev^2/sum(pca$sdev^2)
  if (!all(intgroup %in% names(colData(object)))) {
    stop("the argument 'intgroup' should specify columns of colData(dds)")
  }
  
  intgroup.df <- as.data.frame(colData(object)[, intgroup, 
                                               drop = FALSE])
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], intgroup.df, name = colnames(object))

  if (length(intgroup) > 1) {
    
    if (length(unique(levels(intgroup.df[,2])) > 1)) {
      g <- ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = intgroup.df[, 2], shape = intgroup.df[,1])) +
        geom_point(size = 3) +
        xlab(paste0("PC1: ",round(percentVar[1] * 100), "% variance")) + 
        ylab(paste0("PC2: ", round(percentVar[2] * 100),"% variance")) + 
        coord_fixed() + scale_colour_discrete(name = colnames(intgroup.df)[2]) + 
        scale_shape_discrete(name = colnames(intgroup.df)[1])
    } else {
      g <- ggplot(data = d, aes_string(x = "PC1", y = "PC2", color = intgroup.df[, 1], shape = intgroup.df[,2])) +
        geom_point(size = 3) + xlab(paste0("PC1: ",round(percentVar[1] * 100), "% variance")) + 
        ylab(paste0("PC2: ", round(percentVar[2] * 100),"% variance")) +
        coord_fixed() + scale_colour_discrete(name = colnames(intgroup.df)[1]) + 
        scale_shape_discrete(name = colnames(intgroup.df)[2])
    }
  }  else {
    g <- ggplot(data = d, aes_string(x = "PC1", y = "PC2",color = intgroup.df[, 1])) +
      geom_point(size = 3) + xlab(paste0("PC1: ", round(percentVar[1] * 100), "% variance")) +
      ylab(paste0("PC2: ", round(percentVar[2] * 100), "% variance")) + 
      coord_fixed() + scale_colour_discrete(name = colnames(intgroup.df)[1])
  }
  g <- g + cowplot::theme_cowplot()
  loadings <- as.data.frame(pca$rotation[, 1:2])
  loadings$ID <- rownames(loadings)
  intgroup.df <- apply(intgroup.df, 1, function(x) paste(gsub('`',"",x),collapse = "_"))
  return(list(plot = g, loadings = loadings,sampleAnno = intgroup.df))
}

regressSVs <- function(y, mod, svaobj, P = ncol(mod)) {
  X <- cbind(mod, svaobj$sv)
  Hat <- solve(t(X) %*% X) %*% t(X)
  beta <- (Hat %*% t(y))
  cleany <- y - t(as.matrix(X[, -c(1:P)]) %*% beta[-c(1:P), ])
  return(cleany)
}

tidyExprs <- function(exprs,sampleAnno) {
  
  colnames(exprs) <- sampleAnno
  exprs <- as.data.frame(exprs)
  exprs$ID <- rownames(exprs)
  exprs <- exprs[,c(ncol(exprs),1:(ncol(exprs)-1))]
  
  return(exprs)
}

PCARNASeq <- function (txiData, 
                       sampleTable, 
                       comparisonsTable , 
                       species = c("Human", "Mouse", "Rat", "Pig","Cow", "Horse", "Zebrafish"),
                       technicalReplicates=FALSE, 
                       batchCorrect=FALSE,
                       supervised=FALSE, 
                       pcaPlot="PCA.png",
                       geneInfluence="geneInfluence.txt",
                       normalisedExp = "normalisedExpression.txt" )
{
  
  txi <- readRDS(txiData)
  sampleTable <- read.delim(sampleTable)
  comparisonsTable <- read.delim(comparisonsTable, stringsAsFactors = F)
  sampleTable <- sampleTable[, !grepl("File",colnames(sampleTable))]
  txi$counts <- txi$counts[, match(as.character(sampleTable[,1]), colnames(txi$counts))]
  txi$abundance <- txi$abundance[, match(as.character(sampleTable[,1]), colnames(txi$abundance))]
  txi$length <- txi$length[, match(as.character(sampleTable[,1]), colnames(txi$length))]
  factors <- colnames(sampleTable)[-1]
  designFormula <- as.formula(paste("~", paste(factors, collapse = " + ")))
  
  if (technicalReplicates) {
    dds <- DESeqDataSetFromTximport(txi = txi, colData = sampleTable, 
                                    design = designFormula, ignoreRank = T)
    dds <- collapseReplicates(dds, paste0(sampleTable[, factors[1]], 
                                          sampleTable[, factors[2]]))
    design(dds) <- as.formula(paste("~", paste(factors[-1], 
                                               collapse = " + ")))
  }  else {
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
      resultsTable <- as.data.frame(results(dds, contrast = c(comparisonsTable[1,1],
                                                              comparisonsTable[1,2],
                                                              comparisonsTable[1,3])))
      resultsTable <- resultsTable[order(abs(resultsTable$pvalue)), ]
      controls <- !rownames(dds) %in% rownames(resultsTable)[1:5000]
      svseq <- svaseq(dat, mod, mod0, controls = controls)
    }   else {
      svseq <- svaseq(dat, mod, mod0)
    }
    
    correctedMatrix <- regressSVs(dat, mod, svseq)
    correctedMatrix <- correctedMatrix[!apply(correctedMatrix,1, function(x) any(x < 0)), ]
    dds_new <- DESeqDataSetFromMatrix(round(correctedMatrix), colData = sampleTable, design = designFormula)
    correctedMatrix <- as.matrix(assay(varianceStabilizingTransformation(dds_new)))
    pca <- PCA(varianceStabilizingTransformation(dds_new), intgroup = factors, 
               correctedMatrix = correctedMatrix)
  }  else {
    pca <- PCA(varianceStabilizingTransformation(dds), intgroup = factors, 
               correctedMatrix = NULL)
  }
  
  ggsave(filename = pcaPlot, plot = pca$plot, device = "png",bg="white")
  loadings <- pca$loadings
  
  if (species == "Human") {
    endf <- GenomicFeatures::genes(EnsDb.Hsapiens.v103, return.type = "DataFrame")
  }  else if (species == "Rat") {
    endf <- GenomicFeatures::genes(EnsDb.Rnorvegicus.v103, 
                                   return.type = "DataFrame") 
  }  else if (species == "Cow") {
    endf <- GenomicFeatures::genes(EnsDb.Btaurus.v103, return.type = "DataFrame")
  }  else if (species == "Horse") {
    endf <- GenomicFeatures::genes(EnsDb.Ecaballus.v103, return.type = "DataFrame")
  }  else if (species == "Zebrafish") {
    endf <- GenomicFeatures::genes(EnsDb.Drerio.v103, return.type = "DataFrame")
  }  else if (species == "Pig") {
    endf <- GenomicFeatures::genes(EnsDb.Sscrofa.v103, return.type = "DataFrame")
  }  else {
    endf <- GenomicFeatures::genes(EnsDb.Mmusculus.v103, return.type = "DataFrame")
  }
  
  en2gene <- as.data.frame(endf[, c("gene_id", "gene_name")])
  loadings <- merge(loadings, en2gene, by.x = "ID", by.y = "gene_id")
  loadings <- loadings[, c(1, 4, 2, 3)]
  loadings <- loadings[order(abs(loadings$PC1), decreasing = TRUE),]
  write.table(loadings, file = geneInfluence, col.names = TRUE, 
              row.names = FALSE, sep = "\t", quote = FALSE)
  
  sampleAnno <- pca$sampleAnno
  
  exprs <- tidyExprs(assay(varianceStabilizingTransformation(dds)),sampleAnno)
  exprs <- merge(exprs, en2gene, by.x = "ID", by.y = "gene_id")
  exprs <- exprs[,c(1,ncol(exprs),2:(ncol(exprs)-1))]
  
  write.table(exprs, file = normalisedExp, col.names = TRUE, 
              row.names = FALSE, sep = "\t", quote = FALSE)
  
}

opt <- opt[names(opt) != "help"]
do.call(PCARNASeq, opt)

