#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("ggplot2"))
suppressPackageStartupMessages(library("dplyr"))

option_list <- list()

option_list$differentialExpression <- make_option('--differentialExpression', type='character', help="path to differential expression table")
option_list$homology <- make_option('--homology', type='character', help="path to species to human gene mapping table")
option_list$foldChangeOnly <- make_option('--foldChangeOnly', type='logical', help="does the experiment have replicates")
option_list$species <- make_option('--species', type='character',help="species under study")
option_list$TFs <- make_option('--TFs', type='character',help ="path to list of top genes per transcription factor - from CHEA3")
option_list$TFUp <- make_option('--TFUp', type='character', help="path for output enriched TFs using upregulated genes")
option_list$TFDown <- make_option('--TFDown', type='character',help="path for output enriched TFs using downregulated genes")



opt <- parse_args(OptionParser(option_list=option_list,description = "Find TFs potentially regulating differential expression using CHEA3 approach"))

#get the genes that are coexpressing with a TF
getRegulatedGenes <- function(TF,TFset,datasets){
  
  TFset2<-TFset[[TF]]
  TFset2 <- sort(TFset2[TFset2%in% datasets])
  TFset2<- sort(unique(unlist(TFset2)))
  TFset2 <- paste(TFset2,collapse=" ")
  
  return(TFset2)
}

summariseGenes <- function(genes){
  
  genes <- lapply(genes,function(x) unlist(strsplit(x," ")))
  genes <- sort(unique(unlist(genes)))
  genes <- paste(genes,collapse=" ")
  genes
  
}

getTFRanks <- function(genes,TFsets){
  
  #for each gene set collection get fishers exact test
  ranks <- lapply(TFsets,getSetRanks,genes)
  ranks <- do.call(rbind,ranks)
  ranks <- ranks %>% group_by(TF) %>% summarise(meanRank=mean(values),combinedOverlap=sum(length),Targets=summariseGenes(targets))
  ranks
}

getSetRanks <- function(TFset,genes){
  TF <- sapply(strsplit(names(TFset),"_"),"[[",1)
  scores <- lapply(TFset,getHyper,genes)
  targets <- sapply(scores, "[[",2)
  lengths <- sapply(scores, "[[",3)
  scores <- sapply(scores, "[[",1)
  
  #calculate scaled ranks
  scores <- stack(scores)
  scores$TF <- TF
  scores$length <- lengths
  scores$targets <- sapply(targets,function(x) paste(sort(x),collapse=" "))
  scores <- scores %>% group_by(TF) %>% slice(which.min(values))
  scores$values <- rank(scores$values)
  return(scores)
}

getHyper <- function(query,geneSet){
  
  a <- length(query)
  b<- length(geneSet)
  n <- 20000
  overlap <- intersect(query,geneSet)
  t <- length(overlap)
  
  return(list(sum(dhyper(t:b, a, n - a, b)),overlap,t)) 
  
}

chea3 <- function (differentialExpression, homology = NULL , foldChangeOnly = FALSE,
                   species = c("Human", "Mouse", "Rat", "Horse", "Zebrafish", "Cow", "Pig"), topSig=300,
                   TFs, TFUp = "TFUp.txt", TFDown = "TFDown.txt") {

  
  differentialExpression <- na.omit(read.delim(differentialExpression))
  
  #convert the gene names to human if needed
  if (species != "Human") {
    homology <- read.delim(homology, stringsAsFactors = F)
    differentialExpression <- differentialExpression[differentialExpression[, 1] %in% homology[, 3], ]
    homology <- homology[homology[, 3] %in% differentialExpression[, 1], -1:-2]
    colnames(homology)[2] <- "gene_name.human"
    differentialExpression <- merge(differentialExpression, homology, by.x = 1, by.y = 1)
    differentialExpression[, 1] <- differentialExpression[, "gene_name.human"]
  }
  
  colnames(differentialExpression)[1:2] <- c("GeneSymbol","log2FC")
  
  if (foldChangeOnly == TRUE) {
    differentialExpression <- differentialExpression %>% 
      group_by(GeneSymbol) %>% dplyr::slice(which.max(abs(log2FC))) %>% 
      as.data.frame
    
    upRegulated <- differentialExpression[order(-differentialExpression[,2])[1:topSig] ,1]
    downRegulated <-  differentialExpression[order(differentialExpression[,2])[1:topSig] ,1]
    
  }  else {
    colnames(differentialExpression)[3] <- c("padj")
    differentialExpression <- differentialExpression %>% 
      group_by(GeneSymbol) %>% slice(which.min(padj)) %>% 
      as.data.frame %>% na.omit()
    upRegulated <- differentialExpression[order(-differentialExpression[,2])[1:topSig] ,1]
    downRegulated <-  differentialExpression[order(differentialExpression[,2])[1:topSig] ,1]
    
  }
  
  #read in the TF
  TFs <- readRDS(TFs)
  
  TFRanks.up <- getTFRanks(upRegulated,TFs)
  TFRanks.down <- getTFRanks(downRegulated,TFs)
  
  write.table(TFRanks.up, file=TFUp, sep = "\t",row.names=FALSE,quote = FALSE,col.names = TRUE)
  write.table(TFRanks.down, file=TFDown, sep = "\t",row.names=FALSE,quote = FALSE,col.names = TRUE)
  
  
}


opt <- opt[names(opt) != "help"]
do.call(chea3, opt)
