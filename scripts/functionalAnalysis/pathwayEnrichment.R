#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("goseq"))
suppressPackageStartupMessages(library("dplyr"))


option_list <- list()

option_list$differentialExpression <- make_option('--differentialExpression', type='character',help="path to differential expression table")
option_list$reactomePathways <- make_option('--reactomePathways', type='character', help ="list of pathway gene sets")
option_list$geneExonLengths <- make_option('--geneExonLengths', type='character', help = "path to effective exon length data for genes")
option_list$foldChangeOnly <- make_option('--foldChangeOnly', type='logical', help= "does the data have replicate samples?" )
option_list$foldchange <- make_option('--foldchange', type='numeric', help ="fold change threshold for defining differentially expressed genes")
option_list$padj <- make_option('--padj', type='numeric', help= "p-value threshold for defining differentially expressed genes")
option_list$enrichedPathways <- make_option('--enrichedPathways', type='character', help = "path for enriched pathway results")



opt <- parse_args(OptionParser(option_list=option_list,description = "Find pathways enriched in differentially expressed genes"))


getEnrichedPathways <- function(data, sigData, pathway,geneExonLengths, 
                                species, threshold = 0.05) {
  data <- as.data.frame(data)
  sigData <- as.data.frame(sigData)
  genes <- ifelse(as.data.frame(data)[, 1] %in% as.data.frame(sigData)[, 1], 1, 0)
  names(genes) <- data[, 1]
  
  geneExonLengths <- geneExonLengths[match(data[, 1], geneExonLengths$gene_name),]
  x <- nullp(genes, bias.data = geneExonLengths$length, plot.fit = T)
  PATHWAY = goseq(x, gene2cat = pathway)
  PATHWAY$padj = p.adjust(PATHWAY$over_represented_pvalue, method = "BH")
  PATHWAY$Coverage = PATHWAY$numDEInCat/PATHWAY$numInCat * 100
  PATHWAY$Adjpvaluelog = -log10(PATHWAY$padj)
  PATHWAY$Term <- unlist(PATHWAY$Term)
  PATHWAY.sig = PATHWAY[PATHWAY$padj <= threshold, ]
  if (nrow(PATHWAY.sig) == 0) {
    print("no significant pathways")
    quit(status=3)
  }
  pathwayResults = list()
  for (i in 1:nrow(PATHWAY.sig)) {
    pathwayTerm = PATHWAY.sig$category[i]
    termIDs = pathway[[pathwayTerm]]
    sig = sigData[sigData[, 1] %in% termIDs, ]
    pathwayResults[[pathwayTerm]] = sig[, 1]
  }
  names(pathwayResults) = PATHWAY.sig$category
  pathwayResults = lapply(pathwayResults, function(x) paste(x, 
                                                            sep = "", collapse = " "))
  pathwayResults = data.frame(Term = names(pathwayResults), 
                              Genes = unlist(pathwayResults), Adj.pvalue = PATHWAY.sig$padj, 
                              Coverage = PATHWAY.sig$Coverage, Adjpvaluelog = PATHWAY.sig$Adjpvaluelog)
  return(pathwayResults)
}

pathwayEnrichment <- function (differentialExpression, reactomePathways, geneExonLengths, foldChangeOnly=FALSE, 
    foldchange = 1.5, padj = 0.05, enrichedPathways = "enrichedPathways.txt") {    

    differentialExpression <- na.omit(read.delim(differentialExpression))
    geneExonLengths <- read.delim(geneExonLengths)
    
    colnames(differentialExpression)[1:2] <- c("GeneSymbol","log2FC")
    if (foldChangeOnly) {
        differentialExpression <- differentialExpression %>% 
            group_by(GeneSymbol) %>% slice(which.max(abs(log2FC))) %>% 
            as.data.frame
        differentialExpression.sig <- na.omit(differentialExpression[abs(differentialExpression[,2]) >= log2(foldchange), ])
    }    else {
      colnames(differentialExpression)[3] <- c("padj")
        differentialExpression <- differentialExpression %>% 
            group_by(GeneSymbol) %>% slice(which.min(padj)) %>% 
            as.data.frame
        differentialExpression.sig <- na.omit(differentialExpression[abs(differentialExpression[,2]) >= log2(foldchange) & differentialExpression[,3] <= padj, ])
    }
    
    reactomePathways <- readRDS(reactomePathways)
    
    pathwayResults <- getEnrichedPathways(differentialExpression, 
        differentialExpression.sig, reactomePathways,geneExonLengths)
    write.table(pathwayResults, file = enrichedPathways, 
        col.names = TRUE, row.names = FALSE, sep = "\t", quote = FALSE)
}

opt <- opt[names(opt) != "help"]
do.call(pathwayEnrichment, opt)

