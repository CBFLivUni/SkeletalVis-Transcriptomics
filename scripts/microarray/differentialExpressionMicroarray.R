suppressPackageStartupMessages(library(oligo))
suppressPackageStartupMessages(library(affy))
suppressPackageStartupMessages(library(gcrma))
suppressPackageStartupMessages(library(affyPLM))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(vsn))
suppressPackageStartupMessages(library(sva))
suppressPackageStartupMessages(library(zingeR))

suppressPackageStartupMessages(library("optparse"))

option_list <- list()

option_list$inputfile <- make_option('--inputfile', type='character',help="path to expressionSet RDS")
option_list$comparisonsTable <- make_option('--comparisonsTable', type='character',help="Table specifying the Factor, Numerator and Denominators of interest")
option_list$platform <- make_option('--platform', type='character',help="Type of microarray platform")
option_list$annotationFile <- make_option('--annotationFile', type='character',help="Annotation database")
option_list$foldChangeOnly <- make_option('--foldChangeOnly', type='logical',help="Fold change only?")
option_list$offset <- make_option('--offset', type='integer',help="offset for Illumina data")
option_list$foldChangeTable <- make_option('--foldChangeTable', type='character',help="path for output differential expression table")
option_list$sampleSizeFile <- make_option('--sampleSizeFile', type='character', default='sampleSizes.txt', help="path for output sample sizes table")
option_list$covariates <- make_option('--covariates', type='character', default=NULL, help="Comma-separated list of covariate column names (e.g., 'patient,sex,age')")
option_list$covariate_types <- make_option('--covariate_types', type='character', default=NULL, help="Comma-separated list of covariate types: 'factor' or 'continuous' (e.g., 'factor,factor,continuous')")

opt <- parse_args(OptionParser(option_list=option_list))

getFoldChangeDataFrame <- function(eset, exprs, condition, 
                                   numerator, denominator) {
  numerator.ind <- which(pData(eset)[, condition] %in% 
                           numerator)
  denominator.ind <- which(pData(eset)[, condition] %in% 
                             denominator)
  if (length(numerator.ind) > 1) {
    numerator.dat <- rowMeans(exprs[, numerator.ind])
  } else {
    numerator.dat <- exprs[, numerator.ind]
  }
  if (length(denominator.ind) > 1) {
    denominator.dat <- rowMeans(exprs[, denominator.ind])
  } else {
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

  data <- independentFiltering(data, filter = data$AveExpr, 
                               objectType = "limma")
  data <- data[, c("logFC", "padjFilter")]

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


# Fixed getSampleSizes that handles 2C and 1C
getSampleSizes <- function(eset, condition = NULL, numerator, denominator, covariates = NULL) {
  if (inherits(eset, "ExpressionSet")) {
    sample_data <- pData(eset)
  } else if (!is.null(eset$targets)) {
    sample_data <- as.data.frame(eset$targets, stringsAsFactors = FALSE)
  } else {
    stop("Unable to find sample metadata in 'eset'.")
  }
  
  covariate_str <- if (!is.null(covariates) && length(covariates) > 0) {
    paste(covariates, collapse = ", ")
  } else {
    "None"
  }
  
  # Two-colour handling: count across Cy5 and Cy3
  if (all(c("Cy5", "Cy3") %in% colnames(sample_data))) {
    numerator_N <- sum(sample_data$Cy5 %in% numerator, na.rm = TRUE) + sum(sample_data$Cy3 %in% numerator, na.rm = TRUE)
    denominator_N <- sum(sample_data$Cy5 %in% denominator, na.rm = TRUE) + sum(sample_data$Cy3 %in% denominator, na.rm = TRUE)
    numerator_or_den_index <- which(sample_data$Cy5 %in% numerator | sample_data$Cy3 %in% numerator |
                                      sample_data$Cy5 %in% denominator | sample_data$Cy3 %in% denominator)
    total_N <- length(unique(numerator_or_den_index))
    factor_label <- "Cy5/Cy3"
  } else {
    if (is.null(condition) || !(condition %in% colnames(sample_data))) {
      condition_mn <- make.names(condition)
      if (condition_mn %in% colnames(sample_data)) {
        condition <- condition_mn
      } else {
        stop(paste("Condition column", condition, "not found in sample metadata."))
      }
    }
    numerator.ind <- which(make.names(as.character(sample_data[[condition]])) %in% make.names(as.character(numerator)))
    denominator.ind <- which(make.names(as.character(sample_data[[condition]])) %in% make.names(as.character(denominator)))
    numerator_N <- length(numerator.ind)
    denominator_N <- length(denominator.ind)
    total_N <- numerator_N + denominator_N
    factor_label <- condition
  }
  
  return(data.frame(
    Comparison = paste(numerator, "vs", denominator),
    Factor = factor_label,
    Numerator = numerator,
    Numerator_N = as.integer(numerator_N),
    Denominator = denominator,
    Denominator_N = as.integer(denominator_N),
    Total_N = as.integer(total_N),
    Covariates = covariate_str,
    stringsAsFactors = FALSE
  ))
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
                                                  offset = 0, foldChangeTable = "foldChangeTable.txt", sampleSizeFile = "sampleSizes.txt", covariates = NULL, 
                                                  covariate_types = NULL) 
{
  suppressPackageStartupMessages(library(annotationFile, character.only = TRUE))
  
  # Parse covariates if provided
  covariate_list <- NULL
  covariate_type_list <- NULL
  if (!is.null(covariates) && covariates != "") {
    covariate_list <- strsplit(covariates, ",")[[1]]
    covariate_list <- trimws(covariate_list)
    
    # Parse covariate types
    if (!is.null(covariate_types) && covariate_types != "") {
      covariate_type_list <- strsplit(covariate_types, ",")[[1]]
      covariate_type_list <- trimws(covariate_type_list)
      
      # Validate that lengths match
      if (length(covariate_list) != length(covariate_type_list)) {
        stop("Number of covariates must match number of covariate types")
      }
      
      # Validate types
      valid_types <- c("factor", "continuous")
      if (!all(covariate_type_list %in% valid_types)) {
        stop("Covariate types must be either 'factor' or 'continuous'")
      }
    } else {
      # Default all to factor if not specified
      covariate_type_list <- rep("factor", length(covariate_list))
      print("Warning: covariate_types not specified, defaulting all covariates to 'factor'")
    }
    
    print(paste("Using covariates:", paste(covariate_list, collapse = ", ")))
    print(paste("Covariate types:", paste(covariate_type_list, collapse = ", ")))
  }
  
  # Load the RDS of expression data
  inputfile <- readRDS(inputfile)
  
  # Load the comparison table
  comparisonsTable <- read.delim(comparisonsTable, header = T, 
                                 stringsAsFactors = F)

  # Normalize and annotate the expression data based on the platform used
  if (platform == "Affy") {
    affy::cdfFromBioC(affy::cdfName(inputfile),lib=".")
    eset <- affy::rma(inputfile)
    exprs <- exprs(eset)
    exprs <- annotateProbes(exprs, annotationFile)
  } else if (platform == "Affy-ST") {
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
  } else if (platform == "Illumina") {
    exprs(inputfile) <- exprs(inputfile) + offset
    eset <- normalize(inputfile, transfn = "log")
    exprs <- log2(exprs(eset))
    if (annotationFile == "lumiMouseIDMapping") {
      rownames(exprs) <- as.data.frame(IlluminaID2nuID(rownames(exprs), 
                                                       lib.mapping = "lumiMouseIDMapping", species = "Mouse"))$Symbol
    } else {
      rownames(exprs) <- as.data.frame(IlluminaID2nuID(rownames(exprs), 
                                                       lib.mapping = "lumiHumanIDMapping", species = "Human"))$Symbol
    }
    exprs <- as.matrix(exprs)
  } else if (platform == "1C-Agilent") {
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
    eset <- eset[which(eset$genes$ProbeName %in% geneIDs$ind), ]
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
  } else if (foldChangeOnly == TRUE) {
    if ("Factor2" %in% colnames(comparisonsTable)) {
      Factor1 <- unique(as.character(comparisonsTable[, 1]))
      Factor2 <- unique(as.character(comparisonsTable[, 2]))
      factors <- paste(Factor1, Factor2, sep = ".")
      pData(eset)[factors] <- paste(pData(eset)[, Factor1], 
                                    pData(eset)[, Factor2], sep = ".")
      factors <- gsub("/| ", ".", factors)
      colnames(pData(eset)) <- gsub("/| ", ".", colnames(pData(eset)))
      factorValues <- pData(eset)[factors]
      factorValues <- as.data.frame(apply(factorValues, 
                                          2, as.factor))
      comparisonsTable[, 1] <- paste(comparisonsTable[, 1], 
                                     comparisonsTable[, 2], sep = ".")
      comparisonsTable[, 1] <- gsub("/| ", ".", comparisonsTable[, 1])
      comparisonsTable <- comparisonsTable[, -2]
    }
    resultsTable <- lapply(1:nrow(comparisonsTable), function(x) getFoldChangeDataFrame(eset, 
                                                                                        exprs, comparisonsTable[x, 1], comparisonsTable[x, 2], comparisonsTable[x, 3]))
  } else {
    # Determine if we have Factor2 column (two-factor design)
    has_factor2 <- "Factor2" %in% colnames(comparisonsTable)
    
    # Clean up column names in pData
    colnames(pData(eset)) <- gsub("/| ,:", ".", colnames(pData(eset)))
    
    if (has_factor2) {
      # Two-factor design: concatenate Factor1 and Factor2
      Factor1 <- unique(as.character(comparisonsTable[, 1]))
      Factor2 <- unique(as.character(comparisonsTable[, 2]))
      factors <- paste(Factor1, Factor2, sep = ".")
      pData(eset)[factors] <- paste(pData(eset)[, Factor1], 
                                    pData(eset)[, Factor2], sep = ".")
      factors <- gsub("/| ", ".", factors)
      factorValues <- pData(eset)[factors]
      
      # Handle genotype notation
      if (any(grepl("+/+|-/-", factorValues[, 1]))) {
        factorValues[, 1] <- gsub(x = factorValues[, 1], pattern = "+/+", 
                                  replacement = "WT", fixed = T)
        factorValues[, 1] <- gsub(x = factorValues[, 1], pattern = "-/-", 
                                  replacement = "KO", fixed = T)
        factorValues[, 1] <- gsub(x = factorValues[, 1], pattern = "+/-", 
                                  replacement = "HET", fixed = T)
      }
      factorValues <- as.data.frame(apply(factorValues, 2, as.factor))
      
      # Update comparisons table
      comparisonsTable[, 1] <- paste(comparisonsTable[, 1], 
                                     comparisonsTable[, 2], sep = ".")
      comparisonsTable[, 1] <- gsub("/| ", ".", comparisonsTable[, 1])
      comparisonsTable <- comparisonsTable[, -2]
      
    } else {
      # Single-factor design
      factors <- as.character(unique(comparisonsTable[, 1]))
      print(paste("Main factor:", factors))
      factors <- gsub("/| ", ".", factors)
      comparisonsTable[, 1] <- gsub("/| ,:", ".", comparisonsTable[, 1])
      factorValues <- pData(eset)[factors]
      
      # Handle genotype notation
      if (any(grepl("+/+|-/-", factorValues[, 1]))) {
        factorValues[, 1] <- gsub(x = factorValues[, 1], pattern = "+/+", 
                                  replacement = "WT", fixed = T)
        factorValues[, 1] <- gsub(x = factorValues[, 1], pattern = "-/-", 
                                  replacement = "KO", fixed = T)
        factorValues[, 1] <- gsub(x = factorValues[, 1], pattern = "+/-", 
                                  replacement = "HET", fixed = T)
      }
      
      # Handle numeric starting values
      if (any(grepl("^[[:digit:]]", factorValues[, 1]))) {
        factorValues <- as.data.frame(apply(factorValues, 2, make.names))
      }
      factorValues <- as.data.frame(apply(factorValues, 2, as.factor))
    }
    
    # Build design formula with covariates
    design_terms <- factors
    if (!is.null(covariate_list)) {
      # Clean covariate names
      covariate_list_clean <- gsub("/| ,:", ".", covariate_list)
      
      # Add covariates to the design (covariates come first)
      design_terms <- c( factors,covariate_list_clean)
      
      # Extract covariate data and apply appropriate type conversion
      covariate_data <- pData(eset)[covariate_list_clean]
      for (i in seq_along(covariate_list_clean)) {
        if (covariate_type_list[i] == "factor") {
          covariate_data[[i]] <- as.factor(covariate_data[[i]])
        } else if (covariate_type_list[i] == "continuous") {
          # Extract numeric values from strings (e.g., "30y" -> 30, "5.5kg" -> 5.5)
          covariate_data[[i]] <- as.numeric(gsub("[^0-9.]", "", as.character(covariate_data[[i]])))
          # Check for any NA values introduced by coercion
          if (any(is.na(covariate_data[[i]]))) {
            warning(paste("Some values in covariate", covariate_list_clean[i], 
                         "could not be converted to numeric and are set to NA"))
          }
        }
      }
      covariate_data <- as.data.frame(covariate_data)
      factorValues <- cbind(factorValues,covariate_data)
    }
    # Create design formula: ~0 + covariate1 + ... + main_factor
    designFormula <- as.formula(paste("~0 +", paste(design_terms, collapse = " + ")))
    print(paste("Design formula:", deparse(designFormula)))
    
    design <- model.matrix(designFormula, data = factorValues)
    colnames(design) <- make.names(colnames(design))
    
    # Batch effect correction using SVA
    # Build null model: includes covariates but not main factor
    if (!is.null(covariate_list)) {
      covariate_list_clean <- gsub("/| ,:", ".", covariate_list)
      mod0_formula <- as.formula(paste("~", paste(covariate_list_clean, collapse = " + ")))
      print(paste("Null model formula:", deparse(mod0_formula)))
      mod0 = model.matrix(mod0_formula, data = factorValues)
    } else {
      mod0 = model.matrix(~1, data = factorValues)
    }
    svafit <- try(sva(exprs, mod = design, mod0 = mod0))
    if(inherits(svafit,"try-error")) {
      n <- num.sv(exprs, mod = design) - 1
      svafit <- sva(exprs, mod = design, mod0 = mod0, n.sv = n)
    }
    
    if (is.null(dim(svafit$sv))) {
      svafit$sv <- as.matrix(svafit$sv)
      colnames(svafit$sv) <- paste0("batch", 1:svafit$n.sv)
    } else {
      colnames(svafit$sv) <- paste0("batch", 1:svafit$n.sv)
    }
    design <- cbind(design, svafit$sv)
    
    print(paste("Design matrix columns:", paste(colnames(design), collapse = ", ")))
    
    # Perform differential expression analysis with limma
    fit <- lmFit(exprs, design)
    
    # Build contrasts based on main factor (ignoring covariates in contrast)
    # Find which columns in design correspond to the main factor
    main_factor_cols <- grep(paste0("^", make.names(factors), ".*"), colnames(design), value = TRUE)
    print(paste("Main factor columns:", paste(main_factor_cols, collapse = ", ")))
    
    contrasts <- apply(comparisonsTable, 1, function(x) paste0(make.names(x[1]), 
                                                               make.names(x[2]), "-", make.names(x[1]), make.names(x[3])))

    print(contrasts)
    print(design)
    contrast.matrix <- makeContrasts(contrasts = contrasts, levels = design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    resultsTable <- lapply(seq_along(contrasts), function(x) getResultsDataFrame(fit2, 
                                                                                 x, comparisonsTable[x, 2], comparisonsTable[x, 3]))
  }
  
  # Combine all the differential expression results
  resultsTable <- do.call(cbind, resultsTable)
  resultsTable$ID <- rownames(exprs)
  resultsTable <- resultsTable[, c(ncol(resultsTable), 1:ncol(resultsTable)-1)]
  write.table(resultsTable, file = foldChangeTable, col.names = TRUE,
              row.names = FALSE, sep = "\t", quote = F)
  
  
  
  # -------------------------
  # Build sample size table (safe for both 1C and 2C)
  # -------------------------
  # Find sample metadata (prefer eset$targets for 2C, otherwise pData)
  sample_metadata <- NULL
  if (!is.null(eset$targets)) {
    sample_metadata <- eset$targets
  } else if ("ExpressionSet" %in% class(eset)) {
    sample_metadata <- pData(eset)
  } else if (!is.null(inputfile$targets)) {
    # fallback: sometimes original input object has targets
    sample_metadata <- inputfile$targets
  } else {
    stop("Could not find sample metadata in 'eset' (neither targets nor pData).")
  }

  has_cy <- all(c("Cy5", "Cy3") %in% colnames(sample_metadata))

  if (has_cy) {
    num_col <- if ("Numerator" %in% colnames(comparisonsTable)) "Numerator" else 3
    den_col <- if ("Denominator" %in% colnames(comparisonsTable)) "Denominator" else 4
    sampleSizeTable <- lapply(1:nrow(comparisonsTable), function(x) {
      getSampleSizes(eset, condition = NULL,
                     numerator = comparisonsTable[x, num_col],
                     denominator = comparisonsTable[x, den_col],
                     covariates = covariate_list)
    })
  } else {
    # determine which column in comparisonsTable is the factor name (fallback to 1)
    col_candidates <- c("Factor", "Condition", "factor", "condition", "Variable")
    factor_col <- col_candidates[col_candidates %in% colnames(comparisonsTable)][1]
    if (is.na(factor_col)) factor_col <- 1
    sampleSizeTable <- lapply(1:nrow(comparisonsTable), function(x) {
      getSampleSizes(eset,
                     condition = comparisonsTable[x, factor_col],
                     numerator = comparisonsTable[x, 2],
                     denominator = comparisonsTable[x, 3],
                     covariates = covariate_list)
    })
  }

  sampleSizeTable <- do.call(rbind, sampleSizeTable)
  write.table(sampleSizeTable, file = sampleSizeFile, col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)

  
}

opt <- opt[names(opt)!="help"]
do.call(differentialExpressionFromMicroarray, opt)
