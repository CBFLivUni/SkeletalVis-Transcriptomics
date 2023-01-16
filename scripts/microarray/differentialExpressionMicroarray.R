#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(oligo))
suppressPackageStartupMessages(library(affy))
suppressPackageStartupMessages(library(gcrma))
suppressPackageStartupMessages(library(affyPLM))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(vsn))
suppressPackageStartupMessages(library(sva))
suppressPackageStartupMessages(library(zingeR))
  
  .libPaths( c(".", .libPaths()) )                          



suppressPackageStartupMessages(library("optparse"))


option_list <- list()

option_list$inputfile <- make_option('--inputfile', type='character',help="path to expressionSet RDS")
option_list$comparisonsTable <- make_option('--comparisonsTable', type='character',help="Table specifying the Factor, Numerator and Denominators of interest")
option_list$platform <- make_option('--platform', type='character',help="Type of microarray platform")
option_list$annotationFile <- make_option('--annotationFile', type='character',help="Annotation database")
option_list$foldChangeOnly <- make_option('--foldChangeOnly', type='logical',help="Fold change only?")
option_list$offset <- make_option('--offset', type='integer',help="offset for Illumina data")
option_list$foldChangeTable <- make_option('--foldChangeTable', type='character',help="path for output differential expression table")


opt <- parse_args(OptionParser(option_list=option_list))


getFoldChangeDataFrame <- function(eset, exprs, condition, 
                                   numerator, denominator) {
  numerator.ind <- which(pData(eset)[, condition] %in% 
                           numerator)
  denominator.ind <- which(pData(eset)[, condition] %in% 
                             denominator)
  if (length(numerator.ind) > 1) {
    numerator.dat <- rowMeans(exprs[, numerator.ind])
  }    else {
    numerator.dat <- exprs[, numerator.ind]
  }
  if (length(denominator.ind) > 1) {
    denominator.dat <- rowMeans(exprs[, denominator.ind])
  }     else {
    denominator.dat <- exprs[, denominator.ind]
  }
  data <- data.frame(log2FC = numerator.dat - denominator.dat)
  colnames(data) <- paste(numerator, denominator, sep = "vs")
  return(data)
}
getResultsDataFrame <- function(fit2, contrast, numerator, 
                                denominator) {
  data <- topTable(fit2, coef = contrast, number = Inf, 
                   sort.by = "none")
  print(dim(data))
  data <- independentFiltering(data, filter = data$AveExpr, 
                               objectType = "limma")
  data <- data[, c("logFC", "padjFilter")]
  print(dim(data))
  colnames(data) <- paste(paste(numerator, denominator, 
                                sep = "vs"), colnames(data), sep = "_")
  return(data)
}

annotateProbes <- function(expMat, annotationFile) {
  annotationFile <- gsub(".db", replacement = "", annotationFile)
  geneIDs <- na.omit(stack(mget(as.character(rownames(expMat)), 
                                get(paste(annotationFile, "SYMBOL", sep = "")), ifnotfound = NA)))
  expMat <- merge(expMat, geneIDs, by.x = "row.names", by.y = "ind")
  expMat$Row.names = NULL
  expMat <- aggregate(expMat[, -ncol(expMat)], by = list(expMat$values), 
                     FUN = median, na.rm = TRUE)
  rownames(expMat) <- expMat$Group.1
  expMat[, 1] <- NULL
  expMat <- as.matrix(expMat)
  return(expMat)
}

make2CModelMatrix <- function(targets, comparisonsTable) {
  designMatrix <- apply(comparisonsTable, 1, function(row) {
    numerator <- row[3]
    denominator <- row[4]
    design <- ifelse(targets$Cy5 == numerator & targets$Cy3 == 
                       denominator, 1, 0)
    design[targets$Cy3 == numerator & targets$Cy5 == 
             denominator] <- -1
    design
  })
  return(designMatrix)
}


differentialExpressionFromMicroarray <- function (inputfile, 
    comparisonsTable, 
    platform = c("Affy", "Affy-ST", "1C-Agilent", 
        "2C-Agilent", "Illumina"), annotationFile = c("bovine.db", 
        "hgu133a.db", "hgu133plus2.db", "hgug4112a.db", "HsAgilentDesign026652.db", 
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
    offset = 0, foldChangeTable = "foldChangeTable.txt") 
{
  suppressPackageStartupMessages(library(annotationFile, character.only = TRUE))
  
    
    #load the RDS of expression data
    inputfile <- readRDS(inputfile)
    
    #load the comparison table
    comparisonsTable <- read.delim(comparisonsTable, header = T, 
        stringsAsFactors = F)
    print("here")
    #normalise and annotate the expression data based on the platform used
    if (platform == "Affy") {
	affy::cdfFromBioC(affy::cdfName(inputfile),lib=".")
        eset <- affy::rma(inputfile)
        exprs <- exprs(eset)
        exprs <- annotateProbes(exprs, annotationFile)
    }    else if (platform == "Affy-ST") {
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
    }    else if (platform == "Illumina") {
        exprs(inputfile) <- exprs(inputfile) + offset
        eset <- normalize(inputfile, transfn = "log")
        exprs <- log2(exprs(eset))
        if (annotationFile == "lumiMouseIDMapping") {
            rownames(exprs) <- as.data.frame(IlluminaID2nuID(rownames(exprs), 
                lib.mapping = "lumiMouseIDMapping", species = "Mouse"))$Symbol
        }        else {
            rownames(exprs) <- as.data.frame(IlluminaID2nuID(rownames(exprs), 
                lib.mapping = "lumiHumanIDMapping", species = "Human"))$Symbol
        }
        exprs <- as.matrix(exprs)
    }    else if (platform == "1C-Agilent") {
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
    }    else if (platform == "2C-Agilent") {
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
    }
    
    if (platform == "2C-Agilent") {
        factors <- as.character(unique(comparisonsTable[, 1]))
        design <- make2CModelMatrix(eset$targets, comparisonsTable)
        fit <- lmFit(eset, design)
        fit2 <- eBayes(fit)
        resultsTable <- lapply(1:nrow(comparisonsTable), function(x) getResultsDataFrame(fit2, 
            x, comparisonsTable[x, 3], comparisonsTable[x, 4]))
        exprs <- eset$M
    }    else if (foldChangeOnly == TRUE) {
        if ("Factor2" %in% colnames(comparisonsTable)) {
            Factor1 <- unique(as.character(comparisonsTable[, 
                1]))
            Factor2 <- unique(as.character(comparisonsTable[, 
                2]))
            factors <- paste(Factor1, Factor2, sep = ".")
            pData(eset)[factors] <- paste(pData(eset)[, Factor1], 
                pData(eset)[, Factor2], sep = ".")
            factors <- gsub("/| ", ".", factors)
            colnames(pData(eset)) <- gsub("/| ", ".", colnames(pData(eset)))
            factorValues <- pData(eset)[factors]
            factorValues <- as.data.frame(apply(factorValues, 
                2, as.factor))
            comparisonsTable[, 1] <- paste(comparisonsTable[, 
                1], comparisonsTable[, 2], sep = ".")
            comparisonsTable[, 1] <- gsub("/| ", ".", comparisonsTable[, 
                1])
            comparisonsTable <- comparisonsTable[, -2]
        }
        resultsTable <- lapply(1:nrow(comparisonsTable), function(x) getFoldChangeDataFrame(eset, 
            exprs, comparisonsTable[x, 1], comparisonsTable[x, 
                2], comparisonsTable[x, 3]))
    }    else {
        if ("Factor2" %in% colnames(comparisonsTable)) {
            Factor1 <- unique(as.character(comparisonsTable[, 
                1]))
            Factor2 <- unique(as.character(comparisonsTable[, 
                2]))
            factors <- paste(Factor1, Factor2, sep = ".")
            pData(eset)[factors] <- paste(pData(eset)[, Factor1], 
                pData(eset)[, Factor2], sep = ".")
            factors <- gsub("/| ", ".", factors)
            colnames(pData(eset)) <- gsub("/| ", ".", colnames(pData(eset)))
            factorValues <- pData(eset)[factors]
            if (any(grepl("+/+|-/-", factorValues[, 1]))) {
                factorValues[, 1] <- gsub(x = factorValues[, 
                  1], pattern = "+/+", replacement = "WT", fixed = T)
                factorValues[, 1] <- gsub(x = factorValues[, 
                  1], pattern = "-/-", replacement = "KO", fixed = T)
                factorValues[, 1] <- gsub(x = factorValues[, 
                  1], pattern = "+/-", replacement = "HET", fixed = T)
            }
            factorValues <- as.data.frame(apply(factorValues, 
                2, as.factor))
            comparisonsTable[, 1] <- paste(comparisonsTable[, 
                1], comparisonsTable[, 2], sep = ".")
            comparisonsTable[, 1] <- gsub("/| ", ".", comparisonsTable[, 
                1])
            comparisonsTable <- comparisonsTable[, -2]
            designFormula <- as.formula(paste("~0 +", paste(factors, 
                collapse = " + ")))
            design <- model.matrix(designFormula, data = factorValues)
            colnames(design) <- make.names(colnames(design))
        }        else if ("Paired2" %in% colnames(comparisonsTable)) {
            Factor1 <- unique(as.character(comparisonsTable[,1]))
            Factor2 <- unique(as.character(comparisonsTable[,2]))
            factors <- c(Factor1, Factor2)
            factors <- gsub("/| ", ".", factors)
            colnames(pData(eset)) <- gsub("/| ,:", ".", colnames(pData(eset)))
            factorValues <- pData(eset)[factors]
            factorValues <- as.data.frame(apply(factorValues,2, as.factor))
            comparisonsTable <- comparisonsTable[, -1]
            comparisonsTable[, 2] <- gsub("/| ,:", ".", comparisonsTable[,2])
            designFormula <- as.formula(paste("~0+", paste(factors[2:1],collapse = " + ")))
            design <- model.matrix(designFormula, data = factorValues)
            colnames(design) <- make.names(colnames(design))
        }        else {
            factors <- as.character(unique(comparisonsTable[,1]))
            print(factors)
            factors <- gsub("/| ", ".", factors)
            colnames(pData(eset)) <- gsub("/| ,:", ".", colnames(pData(eset)))
            comparisonsTable[, 1] <- gsub("/| ,:", ".", comparisonsTable[,1])
            factorValues <- pData(eset)[factors]
            if (any(grepl("+/+|-/-", factorValues[, 1]))) {
                factorValues[, 1] <- gsub(x = factorValues[, 
                  1], pattern = "+/+", replacement = "WT", fixed = T)
                factorValues[, 1] <- gsub(x = factorValues[,1], pattern = "-/-", replacement = "KO", fixed = T)
                factorValues[, 1] <- gsub(x = factorValues[,1], pattern = "+/-", replacement = "HET", fixed = T)
            }
            if (any(grepl("^[[:digit:]]", factorValues[, 1]))) {
                factorValues <- as.data.frame(apply(factorValues,2, make.names))
            }
            factorValues <- as.data.frame(apply(factorValues, 2, as.factor))
            designFormula <- as.formula(paste("~0 +", paste(factors, collapse = " + ")))
            design <- model.matrix(designFormula, data = factorValues)
            colnames(design) <- make.names(colnames(design))
        }
      
      
        #batch effect correct the data using sva
        if (!"Paired2" %in% colnames(comparisonsTable)) {
            mod0 = model.matrix(~1, data = factorValues)
            svafit <- try(sva(exprs, mod = design, mod0 = mod0))
    	    if(inherits(svafit,"try-error")) {
    		n <- num.sv(exprs, mod = design) - 1
    		svafit <- sva(exprs, mod = design, mod0 = mod0,n.sv = n)
    	    }
		
            if (is.null(dim(svafit$sv))) {
                svafit$sv <- as.matrix(svafit$sv)
                colnames(svafit$sv) <- paste0("batch", 1:svafit$n.sv)
            }   else {
                colnames(svafit$sv) <- paste0("batch", 1:svafit$n.sv)
            }
            design <- cbind(design, svafit$sv)
        }
        
      
        #perform differential expression analysis limma for each of the comparisons defined in the comparison table       
        fit <- lmFit(exprs, design)
        contrasts <- apply(comparisonsTable, 1, function(x) paste0(make.names(x[1]), 
            make.names(x[2]), "-", make.names(x[1]), make.names(x[3])))
        contrast.matrix <- makeContrasts(contrasts = contrasts, 
            levels = design)
        fit2 <- contrasts.fit(fit, contrast.matrix)
        fit2 <- eBayes(fit2)
        resultsTable <- lapply(seq_along(contrasts), function(x) getResultsDataFrame(fit2, 
            x, comparisonsTable[x, 2], comparisonsTable[x, 3]))
    }
    
    #combine all the differential expression results
    resultsTable <- do.call(cbind, resultsTable)
    resultsTable$ID <- rownames(exprs)
    resultsTable <- resultsTable[, c(ncol(resultsTable), 1:ncol(resultsTable)-1)]
    write.table(resultsTable, file = foldChangeTable, col.names = TRUE,row.names = FALSE, sep = "\t", quote = F)
}


opt <- opt[names(opt)!="help"]
do.call(differentialExpressionFromMicroarray, opt)


