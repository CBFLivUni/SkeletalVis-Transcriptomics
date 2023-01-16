#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(oligo))
suppressPackageStartupMessages(library(affy))
suppressPackageStartupMessages(library(gcrma))
suppressPackageStartupMessages(library(affyPLM))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(vsn))
suppressPackageStartupMessages(library(sva))
suppressPackageStartupMessages(library(GeoDE))

                
option_list <- list()

option_list$inputfile <- make_option('--inputfile', type='character')
option_list$comparisonsTable <- make_option('--comparisonsTable', type='character')
option_list$platform <- make_option('--platform', type='character')
option_list$annotationFile <- make_option('--annotationFile', type='character')
option_list$foldChangeOnly <- make_option('--foldChangeOnly', type='logical')
option_list$offset <- make_option('--offset', type='integer')
option_list$chrDirTable <- make_option('--chrDirTable', type='character')


opt <- parse_args(OptionParser(option_list=option_list))


regressSVs <- function(y, mod, svaobj, P = ncol(mod)) {
  X <- cbind(mod, svaobj$sv)
  Hat <- solve(t(X) %*% X) %*% t(X)
  beta <- (Hat %*% t(y))
  cleany <- y - t(as.matrix(X[, -c(1:P)]) %*% beta[-c(1:P),])
  return(cleany)
}

getResultsDataFrame <- function(comparisonNum, expTable, eset, 
                                comparisonsTable) {
  contrast <- comparisonsTable[comparisonNum, 1]
  numerator <- comparisonsTable[comparisonNum, 2]
  denominator <- comparisonsTable[comparisonNum, 3]
  keep <- which(pData(eset)[[contrast]] %in% c(numerator, denominator))
  expTable <- expTable[, keep]
  eset <- eset[, keep]
  classes <- as.factor(ifelse(pData(eset)[[contrast]] %in%  numerator, 1, 2))
  expTable <- as.data.frame(expTable)
  expTable <- na.omit(cbind(rownames(expTable), expTable))
  chdir <- chdirAnalysis(expTable, classes, 1, CalculateSig = F)
  results <- stack(chdir$chdirprops$chdir[[1]][,1])[, 2:1]
  results <- results[order(abs(results[,2]), decreasing = T), ]
  results$comparison <- paste(numerator, denominator, sep = "vs")
  results$comparisonNumber <- comparisonNum
  colnames(results)[1:2] <- c("ID", "chrDir")
  
  return(results)
}

annotateProbes <- function(expTable, annotationFile) {
  annotationFile <- gsub(".db", replacement = "", annotationFile)
  geneIDs <- na.omit(stack(mget(as.character(rownames(expTable)), 
                                get(paste(annotationFile, "SYMBOL", sep = "")), ifnotfound = NA)))
  expTable <- merge(expTable, geneIDs, by.x = "row.names", by.y = "ind")
  expTable$Row.names <-  NULL
  expTable <- aggregate(expTable[, -ncol(expTable)], by = list(expTable$values), 
                     FUN = median, na.rm = TRUE)
  rownames(expTable) <- expTable$Group.1
  expTable <-expTable[, -1]
  expTable <- as.matrix(expTable)
  return(expTable)
}

characteristicDirectionFromMicroarray <- function (inputfile, 
    comparisonsTable, 
    platform =c("Affy", "Affy-ST", "1C-Agilent", 
        "2C-Agilent", "Illumina", "Illumina_Exp"), 
    annotationFile = c("bovine.db", "hgu133a.db", 
        "hgu133plus2.db", "hgug4112a.db", "HsAgilentDesign026652.db", 
        "HsAgilentDesign026652.db", "hugene10sttranscriptcluster.db", 
        "hugene20sttranscriptcluster.db", "lumiHumanIDMapping", 
        "lumiMouseIDMapping", "lumiRatIDMapping", "mgug4122a.db", 
        "MmAgilentDesign026655.db", "moe430a.db", "mogene10sttranscriptcluster.db", 
        "mogene20sttranscriptcluster.db", "mouse4302.db", "mouse430a2.db", 
        "porcine.db", "ragene10sttranscriptcluster.db", "ragene20sttranscriptcluster.db", 
        "rat2302.db", "xlaevis.db", "hgu95av2.db", "hgfocus.db", 
        "mgu74av2.db", "u133x3p.db", "primeview.db", "AgilentMouse014868.db", 
        "AgilentMouse028005.db", "AgilentRat014879.db", "AgilentRat028279.db", 
        "AgilentHuman039494.db", "ArrayXHuman.db", "AgilentMouse074809.db", 
        "mta10transcriptcluster.db", "htmg430pm.db", "huex10sttranscriptcluster.db", 
        "mgu74a.db", "hugene21sttranscriptcluster.db", "hta20transcriptcluster.db", 
        "AgilentMouse079303.db"), foldChangeOnly = TRUE, 
    offset = 0, chrDirTable = "chrDirTable.txt") 
{

    suppressPackageStartupMessages(library(annotationFile, character.only = TRUE))

    

    
    inputfile <- readRDS(inputfile)
    
    comparisonsTable <- read.delim(comparisonsTable, header = T, 
        stringsAsFactors = F)
print("here")
    
    if (foldChangeOnly == TRUE) {
        return("need replicates to perform characteristic direction analysis")
    }
    if (platform == "Affy") {
affy::cdfFromBioC(affy::cdfName(inputfile),lib=".")
        eset <- affy::rma(inputfile)
        exprs <- exprs(eset)
        exprs <- annotateProbes(exprs, annotationFile)
    }   else if (platform == "Affy-ST") {
        suppressPackageStartupMessages(library(pd.hugene.1.0.st.v1))
        suppressPackageStartupMessages(library(pd.hugene.2.0.st))
        suppressPackageStartupMessages(library(pd.mogene.1.0.st.v1))
        suppressPackageStartupMessages(library(pd.mogene.2.0.st))
        suppressPackageStartupMessages(library(pd.ragene.1.0.st.v1))
        suppressPackageStartupMessages(library(pd.ragene.2.0.st))
        suppressPackageStartupMessages(library(pd.mta.1.0))
        suppressPackageStartupMessages(library(pd.huex.1.0.st.v2))
        suppressPackageStartupMessages(library(pd.hugene.2.1.st))
        suppressPackageStartupMessages(library(pd.hta.2.0))
        eset <- oligo::rma(inputfile)
        exprs <- exprs(eset)
        exprs <- annotateProbes(exprs, annotationFile)
    }   else if (platform == "Illumina") {
        exprs(inputfile) <- exprs(inputfile) + offset
        eset <- normalize(inputfile, transfn = "log")
        exprs <- log2(exprs(eset))
        if (annotationFile == "lumiMouseIDMapping") {
            rownames(exprs) <- as.data.frame(IlluminaID2nuID(rownames(exprs), 
                lib.mapping = "lumiMouseIDMapping", species = "Mouse"))$Symbol
        }    else {
            rownames(exprs) <- as.data.frame(IlluminaID2nuID(rownames(exprs), 
                lib.mapping = "lumiHumanIDMapping", species = "Human"))$Symbol
        }
        exprs <- aggregate(exprs, by = list(rownames(exprs)), 
            FUN = median, na.rm = TRUE)
        rownames(exprs) <- exprs$Group.1
        exprs[, 1] <- NULL
        exprs <- as.matrix(exprs)
    }  else if (platform == "1C-Agilent") {
        eset <- limma::normalizeBetweenArrays(inputfile)
        rownames(eset$E) <- eset$genes$ProbeName
        rownames(eset$targets) <- colnames(eset$E)
        eset$E <- aggregate(eset$E, by = list(eset$genes$ProbeName), 
            FUN = median, na.rm = TRUE)
        rownames(eset$E) <- eset$E$Group.1
        eset$E$Group.1 = NULL
        eset <- ExpressionSet(assayData = as.matrix(eset$E), 
            phenoData = new("AnnotatedDataFrame", eset$targets))
        exprs <- exprs(eset)
        exprs <- annotateProbes(exprs, annotationFile)
    } else if (platform == "2C-Agilent") {
        eset <- limma::backgroundCorrect(inputfile, method = "normexp", 
            offset = 50)
        eset <- limma::normalizeWithinArrays(eset, method = "loess")
        eset <- limma::normalizeBetweenArrays(eset, method = "Aquantile")
        annotationFile <- gsub(".db", replacement = "", annotationFile)
        geneIDs <- na.omit(stack(mget(as.character(eset$genes$ProbeName), 
            get(paste(annotationFile, "SYMBOL", sep = "")), ifnotfound = NA)))
        eset <- eset[which(eset$genes$ProbeName %in% geneIDs$ind), 
            ]
        eset$genes$GeneName <- geneIDs$values
        eset <- avereps(eset, eset$genes$GeneName)
        Cy3 <- comparisonsTable[1, "Cy3"]
        Cy5 <- comparisonsTable[1, "Cy5"]
        eset$targets <- data.frame(FileName = colnames(eset$M), 
            Cy5 = gsub(" +", " ", eset$targets[, Cy5]), Cy3 = gsub(" +", 
                " ", eset$targets[, Cy3]))
        exprs <- exprs.MA(eset)
        targets <- data.frame(as.vector(t(cbind(as.character(eset$targets$Cy5), 
            as.character(eset$targets$Cy3)))))
        colnames(targets) <- comparisonsTable$Factor[1]
        rownames(exprs) <- rownames(eset$M)
        eset <- ExpressionSet(assayData = as.matrix(exprs), phenoData = new("AnnotatedDataFrame", 
            as.data.frame(targets)))
    }
    
    
    if (platform == "2C-Agilent") {
        comparisonsTable <- comparisonsTable[, c(5, 4, 3)]
        resultsTable <- lapply(1:nrow(comparisonsTable), function(x) getResultsDataFrame(x, 
            exprs, eset, comparisonsTable))
        resultsTable <- do.call(rbind, resultsTable)
        write.table(resultsTable, file = chrDirTable, col.names = T, 
            row.names = F, sep = "\t", quote = F)
        return()
    }   else if ("Factor2" %in% colnames(comparisonsTable)) {
        Factor1 <- unique(as.character(comparisonsTable[, 1]))
        Factor2 <- unique(as.character(comparisonsTable[, 2]))
        factors <- paste(Factor1, Factor2, sep = ".")
        pData(eset)[factors] <- paste(pData(eset)[, Factor1],pData(eset)[, Factor2], sep = ".")
        factors <- gsub("/| ", ".", factors)
        colnames(pData(eset)) <- gsub("/| ", ".", colnames(pData(eset)))
        factorValues <- pData(eset)[factors]
        if (any(grepl("+/+|-/-", factorValues[, 1]))) {
            factorValues[, 1] <- gsub(x = factorValues[, 1], 
                pattern = "+/+", replacement = "WT", fixed = T)
            factorValues[, 1] <- gsub(x = factorValues[, 1], 
                pattern = "-/-", replacement = "KO", fixed = T)
            factorValues[, 1] <- gsub(x = factorValues[, 1], 
                pattern = "+/-", replacement = "HET", fixed = T)
        }
        factorValues <- as.data.frame(apply(factorValues, 2,as.factor))
        comparisonsTable[, 1] <- paste(comparisonsTable[, 1], comparisonsTable[, 2], sep = ".")
        comparisonsTable[, 1] <- gsub("/| ", ".", comparisonsTable[,1])
        comparisonsTable <- comparisonsTable[, -2]
        designFormula <- as.formula(paste("~0 +", paste(factors,collapse = " + ")))
        design <- model.matrix(designFormula, data = factorValues)
        colnames(design) <- make.names(colnames(design))
        
    }   else if ("Paired2" %in% colnames(comparisonsTable)) {
        Factor1 <- unique(as.character(comparisonsTable[, 1]))
        Factor2 <- unique(as.character(comparisonsTable[, 2]))
        factors <- c(Factor1, Factor2)
        factors <- gsub("/| ", ".", factors)
        colnames(pData(eset)) <- gsub("/| ,:", ".", colnames(pData(eset)))
        factorValues <- pData(eset)[factors]
        factorValues <- as.data.frame(apply(factorValues, 2, as.factor))
        comparisonsTable <- comparisonsTable[, -1]
        comparisonsTable[, 2] <- gsub("/| ,:", ".", comparisonsTable[,2])
        designFormula <- as.formula(paste("~0+", paste(factors[2:1], collapse = " + ")))
        design <- model.matrix(designFormula, data = factorValues)
        colnames(design) <- make.names(colnames(design))
        comparisonsTable$Paired2 <- make.names(comparisonsTable$Paired2)
        pData(eset)[, Factor2] <- make.names(pData(eset)[, Factor2])
    }   else {
        factors <- as.character(unique(comparisonsTable[, 1]))
        print(factors)
        factors <- gsub("/| ", ".", factors)
        colnames(pData(eset)) <- gsub("/| ,:", ".", colnames(pData(eset)))
        comparisonsTable[, 1] <- gsub("/| ,:", ".", comparisonsTable[,1])
        factorValues <- pData(eset)[factors]
        if (any(grepl("+/+|-/-|+", factorValues[, 1]))) {
            factorValues[, 1] <- gsub(x = factorValues[, 1], 
                pattern = "+/+", replacement = "WT", fixed = T)
            factorValues[, 1] <- gsub(x = factorValues[, 1], 
                pattern = "-/-", replacement = "KO", fixed = T)
            factorValues[, 1] <- gsub(x = factorValues[, 1], 
                pattern = "+/-", replacement = "HET", fixed = T)
            factorValues[, 1] <- gsub(x = factorValues[, 1], 
                pattern = "+", replacement = ".", fixed = T)
        }
        if (any(grepl("^[[:digit:]]", factorValues[, 1]))) {
            factorValues <- as.data.frame(apply(factorValues, 
                2, make.names))
        }
        factorValues <- as.data.frame(apply(factorValues, 2, 
            as.factor))
        designFormula <- as.formula(paste("~0 +", paste(factors, 
            collapse = " + ")))
        design <- model.matrix(designFormula, data = factorValues)
        colnames(design) <- make.names(colnames(design))
    }
    if (!"Paired2" %in% colnames(comparisonsTable)) {
        mod0 = model.matrix(~1, data = factorValues)
        svafit <- try(sva(exprs, mod = design, mod0 = mod0))
        if(inherits(svafit,"try-error")) {
          n <- num.sv(exprs, mod = design) - 1
          svafit <- sva(exprs, mod = design, mod0 = mod0,n.sv = n)
        }
        if (!is.null(dim(svafit$sv))) {
            exprs <- regressSVs(exprs, design, svafit)
        }
        if (any(grepl("+/+|-/-|+|-", pData(eset)[[factors]]))) {
            pData(eset)[[factors]] <- gsub(x = pData(eset)[[factors]], 
                pattern = "+/+", replacement = "WT", fixed = T)
            pData(eset)[[factors]] <- gsub(x = pData(eset)[[factors]], 
                pattern = "-/-", replacement = "KO", fixed = T)
            pData(eset)[[factors]] <- gsub(x = pData(eset)[[factors]], 
                pattern = "+/-", replacement = "HET", fixed = T)
            pData(eset)[[factors]] <- gsub(x = pData(eset)[[factors]], pattern = "+", replacement = ".", fixed = T)
            pData(eset)[[factors]] <- gsub(x = pData(eset)[[factors]], pattern = "-", replacement = ".", fixed = T)
            comparisonsTable[, 2] <- gsub(x = comparisonsTable[, 
                2], pattern = "-", replacement = ".", fixed = T)
            comparisonsTable[, 3] <- gsub(x = comparisonsTable[, 
                3], pattern = "-", replacement = ".", fixed = T)
            comparisonsTable[, 2] <- gsub(x = comparisonsTable[, 
                2], pattern = "+", replacement = ".", fixed = T)
            comparisonsTable[, 3] <- gsub(x = comparisonsTable[, 
                3], pattern = "+", replacement = ".", fixed = T)
        }
        if (any(grepl("^[[:digit:]]", pData(eset)[[factors]]))) {
            pData(eset)[[factors]] <- make.names(pData(eset)[[factors]])
        }
    }
    
    
    resultsTable <- lapply(1:nrow(comparisonsTable), function(x) getResultsDataFrame(x, exprs, eset, comparisonsTable))
    resultsTable <- do.call(rbind, resultsTable)
    write.table(resultsTable, file = chrDirTable, col.names = TRUE, row.names = FALSE, sep = "\t", quote = F)
}

opt <- opt[names(opt)!="help"]
do.call(characteristicDirectionFromMicroarray, opt)

