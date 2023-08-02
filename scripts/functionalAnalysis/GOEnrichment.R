#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("goseq"))
suppressPackageStartupMessages(library("GO.db"))
suppressPackageStartupMessages(library("GOSemSim"))
suppressPackageStartupMessages(library("tidyr"))
suppressPackageStartupMessages(library("plotly"))
suppressPackageStartupMessages(library("htmlwidgets"))
suppressPackageStartupMessages(library("dplyr"))

option_list <- list()

option_list$differentialExpression <- make_option('--differentialExpression', type='character', help="path to differential expression table")
option_list$geneLengths <- make_option('--geneLengths', type='character', help="path to gene length table")
option_list$foldChangeOnly <- make_option('--foldChangeOnly', type='logical', help="does the experiment have replicates")
option_list$species <- make_option('--species', type='character', help="species under study")
option_list$foldchange <- make_option('--foldchange', type='numeric', help="threshold for the fold change to define differentially expressed genes")
option_list$padj <- make_option('--padj', type='numeric', help="threshold for the pvalue to define differentially expressed genes")
option_list$enrichedTerms <- make_option('--enrichedTerms', type='character', help ="path for output enrichment results")
option_list$enrichedTermsReduced <- make_option('--enrichedTermsReduced', type='character', help ="path for output simplified enrichment results")
option_list$mdsPlot <- make_option('--mdsPlot', type='character', help= "path for output mds plot")


opt <- parse_args(OptionParser(option_list=option_list,description = "find enriched gene ontology terms in differential expression results"))

getEnrichedGOTerms <- function(diffExp, sigDiffExp, geneLengths, 
                               species, threshold = 0.05) {
  if (species == "Human") {
    suppressPackageStartupMessages(library("org.Hs.eg.db"))
    gene2GO <- AnnotationDbi::select(org.Hs.eg.db, keys(org.Hs.egGO2EG), 
                                     c("ENTREZID", "SYMBOL"), "GOALL")
  }        else if (species == "Mouse") {
    suppressPackageStartupMessages(library("org.Mm.eg.db"))
    gene2GO <- AnnotationDbi::select(org.Mm.eg.db, keys(org.Mm.egGO2EG), 
                                     c("ENTREZID", "SYMBOL"), "GOALL")
  }        else if (species == "Cow") {
    suppressPackageStartupMessages(library("org.Bt.eg.db"))
    gene2GO <- AnnotationDbi::select(org.Bt.eg.db, keys(org.Bt.egGO2EG), 
                                     c("ENTREZID", "SYMBOL"), "GOALL")
  }        else if (species == "Horse") {
    suppressPackageStartupMessages(library("org.Ecaballus.eg.db"))
    gene2GO <- AnnotationDbi::select(org.Ecaballus.eg.db, 
                                     keys(org.Ecaballus.egGO2EG), c("ENTREZID", "SYMBOL"), 
                                     "GOALL")
  }        else if (species == "Zebrafish") {
    
    suppressPackageStartupMessages(library("org.Dr.eg.db"))
    gene2GO <- AnnotationDbi::select(org.Dr.eg.db, keys(org.Dr.egGO2EG), 
                                     c("ENTREZID", "SYMBOL"), "GOALL")
  }        else if (species == "Pig") {
    suppressPackageStartupMessages(library("org.Ss.eg.db"))
    gene2GO <- AnnotationDbi::select(org.Ss.eg.db, keys(org.Ss.egGO2EG), 
                                     c("ENTREZID", "SYMBOL"), "GOALL")
  }        else {
    suppressPackageStartupMessages(library("org.Rn.eg.db"))
    gene2GO <- AnnotationDbi::select(org.Rn.eg.db, keys(org.Rn.egGO2EG), 
                                     c("ENTREZID", "SYMBOL"), "GOALL")
  }
  gene2GO <- unstack(gene2GO[, c(1, 5)])
  genes <- ifelse(diffExp[, 1] %in% sigDiffExp[, 1], 1, 0)
  names(genes) <- diffExp[, 1]
  geneLengths <- geneLengths[match(diffExp[,1], geneLengths$gene_name), ]
  x <- nullp(genes, bias.data = geneLengths$length, plot.fit = FALSE)
  GOTerms <- goseq(x, gene2cat = gene2GO)
  GOTerms <- GOTerms[GOTerms$ontology == "BP", ]
  GOTerms$padj <- p.adjust(GOTerms$over_represented_pvalue, method = "BH")
  GOTerms.sig <- GOTerms[GOTerms$padj <= threshold, ]
  if (nrow(GOTerms.sig) == 0) {
    print("no significant GO terms!")
    quit(status=3)
  }
  GOTerms.sig$enrichment <- GOTerms.sig$numDEInCat/GOTerms.sig$numInCat
  GOResults = list()
  
  for (i in 1:nrow(GOTerms.sig)) {
    GOTerm <- GOTerms.sig$category[i]
    index <- sapply(gene2GO, function(x) GOTerm %in% x)
    termIDs <- names(index[index == "TRUE"])
    sig <- sigDiffExp[sigDiffExp[, 1] %in% termIDs, ]
    GOResults[[GOTerm]] = sig[, 1]
  }
  names(GOResults) = GOTerms.sig$term
  GOResults <- lapply(GOResults, function(x) paste(x, sep = "", collapse = " "))
  GOResults <- data.frame(Term = names(GOResults), ID = GOTerms.sig$category, 
                          Genes = unlist(GOResults), Adj.pvalue = GOTerms.sig$padj, 
                          Enrichment = GOTerms.sig$enrichment)
  if(nrow(GOResults)==1) return(list(GOResults = GOResults,GOResults.reduced=NULL))

  GOResults.reduced <- try(simplify(GOResults, gene2GO,species))
  return(list(GOResults = GOResults, GOResults.reduced = GOResults.reduced))
}
simplify <- function(GORes, gene2GO,species) {
  if (species == "Human") {
    semData <- godata("org.Hs.eg.db", "SYMBOL", "BP")
  }        else if (species == "Mouse") {
    semData <- godata("org.Mm.eg.db", "SYMBOL", "BP")
  }        else if (species == "Zebrafish") {
    semData <- godata("org.Dr.eg.db", "SYMBOL", "BP")
  }        else if (species == "Horse") {
    semData <- godata("org.Mm.eg.db", "SYMBOL", "BP")
  }        else if (species == "Cow") {
    semData <- godata("org.Bt.eg.db", "SYMBOL", "BP")
  }        else if (species == "Pig") {
    semData <- godata("org.Ss.eg.db", "SYMBOL", "BP")
  }        else {
    semData <- godata("org.Rn.eg.db", "SYMBOL", "BP")
  }
  sim <- mgoSim(GORes$ID, GORes$ID, semData = semData, 
                measure = "Rel", combine = NULL)
  sim[is.na(sim)] <- 0
  go1 <- go2 <- similarity <- NULL
  sim.df <- as.data.frame(sim)
  sim.df$go1 <- row.names(sim.df)
  sim.df <- gather(sim.df, go2, similarity, -go1)
  sim.df <- sim.df[!is.na(sim.df$similarity), ]
  sim.df <- sim.df[order(sim.df$similarity, decreasing = T), 
  ]
  sim.df <- sim.df[sim.df$go1 != sim.df$go2, ]
  sim.df <- sim.df[sim.df$similarity > 0.4, ]
  GO2Gene <- unstack(stack(gene2GO)[2:1])
  freq <- sapply(GO2Gene, length)
  freqCutOff <- length(gene2GO) * 0.05
  highFreqTerms <- names(freq[freq > freqCutOff])
  
  sim.df$remove <- apply(sim.df, 1, function(x) {
    if (x[1] %in% highFreqTerms) {
      return(x[1])
    }
    if (x[2] %in% highFreqTerms) {
      return(x[2])
    }            else {
      return(NA)
    }
  })
  
  remove <- na.omit(sim.df$remove)
  sim.df <- sim.df[is.na(sim.df$remove), ]
  sim.df$go1.pval <- GORes$Adj.pvalue[match(sim.df$go1, GORes$ID)]
  sim.df$go2.pval <- GORes$Adj.pvalue[match(sim.df$go2, GORes$ID)]
  childTerms <- as.list(GOBPCHILDREN)
  
  for (i in 1:nrow(sim.df)) {
    if (sim.df[i, "go1"] %in% remove) {
      next
    }
    if (sim.df[i, "go2"] %in% remove) {
      next
    }
    go1.pval <- sim.df[i, "go1.pval"]
    go2.pval <- sim.df[i, "go2.pval"]
    if (go1.pval == go2.pval) {
      go1 <- sim.df[i, "go1"]
      go2 <- sim.df[i, "go2"]
      if (go2 %in% childTerms[[go1]]) {
        remove <- c(remove, go2)
        next
      }                else if (go1 %in% childTerms[[go2]]) 
        remove <- c(remove, go1)
      next
    }
    remove <- c(remove, sim.df[i, which.max(c(go1.pval, 
                                              go2.pval))])
  }
  GORes.filt <- GORes[!GORes$ID %in% remove, ]
  sim.filt <- sim[as.character(GORes.filt$ID), as.character(GORes.filt$ID)]
  fit <- cmdscale(1 - sim.filt, eig = TRUE, k = 2)
  x <- fit$points[, 1]
  y <- fit$points[, 2]
  GORes.filt.plot <- GORes.filt
  GORes.filt.plot$x <- x
  GORes.filt.plot$y <- y
  GORes.filt.plot$log10Adjpvalue <- -log10(GORes.filt.plot$Adj.pvalue)
  GO.MDS <- plot_ly(GORes.filt.plot, x = GORes.filt.plot$x, 
                    y = GORes.filt.plot$y, mode = "markers", type = "scatter", 
                    color = GORes.filt.plot$log10Adjpvalue, size = log2(GORes.filt.plot$Enrichment), 
                    text = paste("Term: ", GORes.filt$Term, ""), marker = list(sizeref = 0.05)) %>% 
    colorbar(title = "-log10 Adj.pvalue")
  return(list(GORes.filt, GO.MDS))
}

GOEnrichment <- function (differentialExpression, geneLengths, foldChangeOnly = TRUE, 
    species = c("Human", "Mouse", "Rat", "Horse","Zebrafish", "Cow", "Pig"), foldchange = 1.5, padj = 0.05,
    enrichedTerms = "enrichedTerms.txt", enrichedTermsReduced = "enrichedTerms.reduced.txt", mdsPlot = "GO.MDS.html") {

    differentialExpression <- na.omit(read.delim(differentialExpression))
    geneLengths <- read.delim(geneLengths)
    colnames(differentialExpression)[1:2] <- c("GeneSymbol","log2FC")

    if (foldChangeOnly == TRUE) {
        differentialExpression <- differentialExpression %>% 
            group_by(GeneSymbol) %>% slice(which.max(abs(log2FC))) %>% 
            as.data.frame
        differentialExpression.sig <- na.omit(differentialExpression[abs(differentialExpression[, 2]) >= log2(foldchange), ])
    }    else {
      colnames(differentialExpression)[3] <- "padj"
        differentialExpression <- differentialExpression %>% 
            group_by(GeneSymbol) %>% slice(which.min(padj)) %>% 
            as.data.frame
        differentialExpression.sig <- na.omit(differentialExpression[abs(differentialExpression[, 
            2]) >= log2(foldchange) & differentialExpression[, 3] <= padj, ])
    }
    GOResultsList <- getEnrichedGOTerms(differentialExpression, 
        differentialExpression.sig, geneLengths, species)
    GOResults <- GOResultsList[[1]]
    write.table(GOResults, file = enrichedTerms, col.names = T, 
        row.names = F, sep = "\t", quote = F)
    
    if (!inherits(GOResultsList[[2]], "try-error") & !is.null(GOResultsList[[2]])) {
        write.table(GOResultsList[[2]][[1]], file = enrichedTermsReduced, 
            col.names = T, row.names = F, sep = "\t", quote = F)
        GO.MDS <- GOResultsList[[2]][[2]]
        GO.MDS <- saveWidget(GO.MDS, file = mdsPlot)
    }
}

opt <- opt[names(opt) != "help"]
do.call(GOEnrichment, opt)
