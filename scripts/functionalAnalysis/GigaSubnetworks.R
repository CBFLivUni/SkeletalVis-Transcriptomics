#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("igraph"))
suppressPackageStartupMessages(library("squash"))
suppressPackageStartupMessages(library("intergraph"))
suppressPackageStartupMessages(library("visNetwork"))
suppressPackageStartupMessages(library("goseq"))
suppressPackageStartupMessages(library("dplyr"))

option_list <- list()

option_list$differentialExpression <- make_option('--differentialExpression', type='character')
option_list$interactome <- make_option('--interactome', type='character')
option_list$foldChangeOnly <- make_option('--foldChangeOnly', type='logical')
option_list$species <- make_option('--species', type='character')
option_list$visNetworks <- make_option('--visNetworks', type='character')
option_list$summaryTable <- make_option('--summaryTable', type='character')


opt <- parse_args(OptionParser(option_list=option_list))


getTopGOTerms <- function(geneLists, geneLengths, differentialExpression, 
                          species) {
  if (species == "Human") {
    suppressPackageStartupMessages(library("org.Hs.eg.db"))
    gene2GO <- AnnotationDbi::select(org.Hs.eg.db, keys(org.Hs.egGO2EG), c("ENTREZID", "SYMBOL"), "GOALL")
  }      else if (species == "Mouse") {
    suppressPackageStartupMessages(library("org.Mm.eg.db"))
    gene2GO <- AnnotationDbi::select(org.Mm.eg.db, keys(org.Mm.egGO2EG), c("ENTREZID", "SYMBOL"), "GOALL")
  }      else if (species == "Cow") {
    suppressPackageStartupMessages(library("org.Bt.eg.db"))
    gene2GO <- AnnotationDbi::select(org.Bt.eg.db, keys(org.Bt.egGO2EG), c("ENTREZID", "SYMBOL"), "GOALL")
  }      else if (species == "Horse") {
    suppressPackageStartupMessages(library("org.Ecaballus.eg.db"))
    gene2GO <- AnnotationDbi::select(org.Ecaballus.eg.db,keys(org.Ecaballus.egGO2EG), c("ENTREZID", "SYMBOL"), "GOALL")
  }      else if (species == "Zebrafish") {
    suppressPackageStartupMessages(library("org.Dr.eg.db"))
    gene2GO <- AnnotationDbi::select(org.Dr.eg.db, keys(org.Dr.egGO2EG),c("ENTREZID", "SYMBOL"), "GOALL")
  }      else if (species == "Pig") {
    suppressPackageStartupMessages(library("org.Ss.eg.db"))
    gene2GO <- AnnotationDbi::select(org.Ss.eg.db, keys(org.Ss.egGO2EG),c("ENTREZID", "SYMBOL"), "GOALL")
  }      else {
    suppressPackageStartupMessages(library("org.Rn.eg.db"))
    gene2GO <- AnnotationDbi::select(org.Rn.eg.db, keys(org.Rn.egGO2EG),c("ENTREZID", "SYMBOL"), "GOALL")
  }
  
  gene2GO <- unstack(gene2GO[, c(1, 5)])
  TopGOTerms <- sapply(geneLists, getTopGOTerm, geneLengths, 
                       differentialExpression, gene2GO)
  return(TopGOTerms)
}

getTopGOTerm <- function(geneList, geneLengths, differentialExpression, 
                         gene2GO) {
  genes <- ifelse(differentialExpression[, 1] %in% geneList, 1, 0)
  genes <- as.data.frame(genes)
  colnames(genes) <- "DEgenes"
  genes$bias.data <- 1
  genes$pwf <- 1
  rownames(genes) <- differentialExpression[, 1]
  goTerms <- goseq(genes, gene2cat = gene2GO, method = "Hypergeometric")
  return(goTerms[1, "term"])
}

getVisSubnetworks <- function(subnetTable, nodeList, network) {
  maxVal <- max(abs(subnetTable[, 3])) + 0.05
  minVal <- -maxVal
  breaksList = seq(minVal, maxVal, by = 0.05)
  colfunc <- colorRampPalette(c("green", "white", "red"))
  map <- makecmap(subnetTable[, 3], n = 3, breaks = breaksList, symm = TRUE, colFn = colfunc)
  subnetworks <- lapply(nodeList, function(x) induced.subgraph(network, x))
  visNetworks <- lapply(subnetworks, getVisSubnetwork, map, subnetTable)
  return(visNetworks)
}

getVisSubnetwork <- function(subnetwork, colours, subnetTable) {
  FCs <- subnetTable[match(V(subnetwork)$name, subnetTable[, 1]), 3]
  colours <- cmap(FCs, map = colours)
  nodes <- data.frame(id = V(subnetwork)$name, title = as.character(FCs), 
                      label = V(subnetwork)$name, color = colours, size = 16, 
                      font.size = 22)
  edges <- as.data.frame(get.edgelist(subnetwork))
  edges <- data.frame(from = edges[, 1], to = edges[, 2], color = "black")
  network <- visNetwork(nodes, edges)
  return(network)
}

runGIGA <- function(expressionData, inputGraph, max_number = 20) {
  
  Scores <- data.frame(name = expressionData[, 1], geneRank = rank(-expressionData$Pi), 1:nrow(expressionData))
  adj <- get.adjacency(inputGraph, sparse = TRUE)
  neighboursAdj <- lapply(1:nrow(adj), getNeighboursAdj, adj)
  neighboursAdjSelf <- mapply(c, 1:nrow(adj), neighboursAdj)
  Scores <- Scores[order(Scores[, 3]), ]
  seedCluster <- initaliseCluster(Scores, inputGraph, neighboursAdj, neighboursAdjSelf)
  expandedCluster <- clusterExpand(Scores, seedCluster, completedCluster = c(), inputGraph, max_number, neighboursAdjSelf)
  
  while (length(expandedCluster[["cluster"]]) > 0) {
    seedCluster <- addNewMin(Scores, expandedCluster, 
                             inputGraph, neighboursAdjSelf)
    expandedCluster <- clusterExpand(Scores, seedCluster, 
                                     seedCluster[["completedCluster"]], inputGraph, 
                                     max_number, neighboursAdjSelf)
  }
  
  finalClusters <- getFinalClusters(Scores, expandedCluster, inputGraph)
  clusters <- finalClusters[[1]]
  pvalues <- finalClusters[[2]]
  clusters <- lapply(clusters, function(x) V(inputGraph)$name[x])
  return(list(clusters = clusters, pvalues = pvalues))
}

getNeighboursAdj <- function(i, adj) {
  neighbours <- which(adj[i, ] > 0)
  return(neighbours)
}

initaliseCluster = function(Scores, inputGraph, neighboursAdj, neighboursAdjSelf) {
  
  size <- length(neighboursAdjSelf)
  NodeDataRankNamed <- Scores[, 3]
  names(NodeDataRankNamed) <- Scores$geneRank
  neighbours <- lapply(neighboursAdjSelf, function(x) NodeDataRankNamed[x])
  localMinIndex <- lapply(neighbours, function(x) which.min(as.numeric(names(x))))
  localMin <- which(localMinIndex == 1)
  cluster <- as.list(localMin)
  names(cluster) <- V(inputGraph)$name[localMin]
  oldcluster <- cluster
  cluster <- lapply(cluster, function(x) NodeDataRankNamed[x])
  oldpvalue <- rep(1, length(cluster))
  names(oldpvalue) <- names(cluster)
  clustermax <- localMin
  neighbours <- neighboursAdj[localMin]
  neighbours <- sapply(neighbours, function(x) NodeDataRankNamed[x])
  neighbours <- sapply(neighbours, function(x) x[sort.list(as.numeric(names(x)))])
  neighbours <- lapply(neighbours, "[", 1)
  cluster <- mapply(c, cluster, neighbours, SIMPLIFY = FALSE)
  clustermax <- sapply(cluster, function(x) names(x[2]))
  return(list(clustermax = clustermax, oldpvalue = oldpvalue, 
              cluster = cluster, oldcluster = oldcluster))
}

clusterExpand = function(Scores, seedCluster, completedCluster, 
                         inputGraph, max_number, neighboursAdjSelf) {
  
  NodeDataRankNamed = Scores[, 3]
  names(NodeDataRankNamed) <- Scores$geneRank
  clustermax <- seedCluster[["clustermax"]]
  oldpvalue <- seedCluster[["oldpvalue"]]
  cluster <- seedCluster[["cluster"]]
  oldcluster <- seedCluster[["oldcluster"]]
  cluster_new <- cluster
  tobePValued <- c()
  size <- vcount(inputGraph)
  
  while (length(cluster_new) > 0) {
    neigbours <- getNeigbours(cluster_new, neighboursAdjSelf)
    lastexpansion <- cluster_new
    neigbours <- lapply(neigbours, function(x) NodeDataRankNamed[x])
    clustermax <- as.numeric(clustermax)
    expansion <- lapply(seq_along(neigbours), function(x) neigbours[[x]][(as.numeric(names(neigbours[[x]])) <= clustermax[x])])
    names(expansion) <- names(cluster_new)
    expansion_tooBig <- lapply(expansion, function(x) length(x) <  max_number)
    cluster_new <- expansion[unlist(expansion_tooBig)]
    completedCluster <- c(completedCluster, oldcluster[!unlist(expansion_tooBig)])
    oldcluster <- oldcluster[unlist(expansion_tooBig)]
    clustermax <- clustermax[unlist(expansion_tooBig)]
    clusterlength <- lapply(cluster_new, length)
    oldclusterlength <- lapply(lastexpansion, length)
    continueCluster <- lapply(names(cluster_new), function(x) clusterlength[[x]] > oldclusterlength[[x]])
    if (length(continueCluster > 1)) {
      tobePValued <- c(tobePValued, cluster_new[!unlist(continueCluster)])
    }
    cluster_new <- cluster_new[unlist(continueCluster)]
    clustermax <- clustermax[unlist(continueCluster)]
    oldcluster <- oldcluster[unlist(continueCluster)]
  }
  if (length(tobePValued) > 0) {
    oldcluster <- seedCluster[["oldcluster"]]
    oldcluster <- oldcluster[match(names(tobePValued), 
                                   names(oldcluster))]
    oldpvalue <- oldpvalue[names(oldpvalue) %in% names(tobePValued)]
    currentPValue <- lapply(tobePValued, getPvalue, size)
    expansion_improved <- lapply(names(currentPValue), 
                                 function(x) currentPValue[[x]] < oldpvalue[[x]])
    names(expansion_improved) <- names(currentPValue)
    cluster <- tobePValued[unlist(expansion_improved)]
    worseClusters <- oldcluster[!unlist(expansion_improved)]
    completedCluster <- c(completedCluster, worseClusters)
    oldpvalue <- currentPValue[unlist(expansion_improved)]
  }
  else {
    cluster <- NULL
  }
  return(list(oldpvalue = oldpvalue, cluster = cluster, 
              completedCluster = completedCluster))
}

addNewMin <- function(Scores, expandedCluster, inputGraph, neighboursAdjSelf) {
  
  NodeDataRankNamed <- Scores[, 3]
  names(NodeDataRankNamed) <- Scores$geneRank
  cluster <- expandedCluster[["cluster"]]
  oldpvalue <- expandedCluster[["oldpvalue"]]
  completedCluster <- expandedCluster[["completedCluster"]]
  oldcluster <- cluster
  newNeighbours <- getNeigbours(cluster, neighboursAdjSelf)
  clustermax <- c()
  newMin <- lapply(newNeighbours, function(x) NodeDataRankNamed[x])
  names(newMin) <- names(cluster)
  newMin <- lapply(newMin, function(x) x[sort.list(as.numeric(names(x)))])
  addition <- lapply(seq_along(newMin), function(x) newMin[[x]][!(newMin[[x]] %in% cluster[[x]])][1])
  cluster <- mapply(c, cluster, addition, SIMPLIFY = FALSE)
  clustermax <- sapply(cluster, function(x) max(as.numeric(names(x))))
  return(list(clustermax = clustermax, oldpvalue = oldpvalue, 
              cluster = cluster, oldcluster = oldcluster, completedCluster = completedCluster))
}

getPvalue <- function(x, size) {
  x <- as.numeric(names(x))
  maxRank <- max(x)
  pvalue <- maxRank/size
  if (length(x) == 1) 
    return(pvalue)
  for (k in 1:(length(x) - 1)) {
    pvalue <- pvalue * ((maxRank - k)/(size - k))
  }
  return(pvalue)
}

getNeigbours = function(cluster, neighboursAdj) {
  newNeighbours <- lapply(cluster, function(x) unlist(neighboursAdj[unlist(x)]))
  newNeighbours <- lapply(newNeighbours, function(x) x[!duplicated(x)])
  return(newNeighbours)
}

getFinalClusters = function(Scores, expandedCluster, inputGraph) {
  size <- vcount(inputGraph)
  NodeDataRankNamed <- Scores[, 3]
  names(NodeDataRankNamed) <- Scores$geneRank
  completedCluster <- expandedCluster[["completedCluster"]]
  completedCluster <- lapply(completedCluster, function(x) NodeDataRankNamed[x])
  completedClusterpValue <- lapply(completedCluster, getPvalue, size)
  completedCluster <- completedCluster[completedClusterpValue <  1/size]
  completedClusterpValue <- lapply(completedCluster, getPvalue, size)
  completedClusterpValue <- data.frame(name = names(completedCluster), 
                                       pvalue = unlist(completedClusterpValue))
  completedClusterpValue <- completedClusterpValue[order(completedClusterpValue$pvalue),]
  completedClusterpValue$localmin <- match(completedClusterpValue$name, V(inputGraph)$name)
  localmin <- completedClusterpValue[, 3]
  finalCluster <- c()
  localminused <- c()
  alreadydone <- c()
  finalPvalues <- c()
  for (i in 1:nrow(completedClusterpValue)) {
    if ((completedClusterpValue[i, 3] %in% localminused) == FALSE) {
      candidateCluster <- completedCluster[[as.character(completedClusterpValue[i, 1])]]
      localminpresent <- candidateCluster[candidateCluster %in% localmin]
      localminused <- c(localminpresent, localminused)
      finalCluster <- c(finalCluster, list(candidateCluster))
      finalPvalues <- c(finalPvalues, completedClusterpValue[i, "pvalue"])
    }
  }
  return(list(finalCluster, finalPvalues))
}


GigaSubnetworks <- function (differentialExpression, 
                             interactome, foldChangeOnly, 
                             species = c("Human", "Mouse", "Rat", "Horse","Zebrafish", "Cow", "Pig"), 
                             visNetworks = "visNetworks.RDS", summaryTable = "summaryTable.txt") {

  differentialExpression <- read.delim(differentialExpression)
  interactome <- read.table(interactome, sep = "\t", as.is = TRUE)
  interactome <- graph.data.frame(interactome, directed = F)
  colnames(differentialExpression)[1:2] <- c("GeneSymbol","log2FC")
  differentialExpression$absFC <- abs(differentialExpression[, 
                                                             2])
  if (foldChangeOnly == T) {
    differentialExpression <- differentialExpression %>% 
      group_by(GeneSymbol) %>% dplyr::slice(which.max(abs(log2FC))) %>% 
      as.data.frame
    differentialExpression$Pi <- differentialExpression$absFC
  }  else {
    colnames(differentialExpression)[3] <- "padj"
    differentialExpression <- differentialExpression %>% 
      group_by(GeneSymbol) %>% dplyr::slice(which.min(padj)) %>% 
      as.data.frame
    differentialExpression$logAdjPval <- -log10(differentialExpression[,3])
    differentialExpression$Pi <- (differentialExpression$absFC)
  }
  
  differentialExpression <- na.omit(differentialExpression)
  presentList <- na.omit(match(differentialExpression[, 1], V(interactome)$name))
  interactome <- induced.subgraph(interactome, presentList)
  interactome <- decompose.graph(interactome)
  interactome <- interactome[[which.max(sapply(interactome, vcount))]]
  presentList <- na.omit(match(V(interactome)$name, differentialExpression[, 1]))
  differentialExpression <- differentialExpression[presentList, ]
  results <- runGIGA(expressionData = differentialExpression, inputGraph = interactome, max_number = 20)
  pvalues <- results$pvalues
  output <- data.frame(subnetworks = rep(1:length(results$clusters), 
                                         sapply(results$clusters, length)), gene_name = unlist(results$clusters))
  
  output <- merge(output, differentialExpression, by.x = "gene_name", by.y = 1, all.x = TRUE)
  output <- output[order(output$subnetworks, output$gene_name), ]
  topGOTerms <- getTopGOTerms(results$clusters, geneLengths, differentialExpression, species)
  resultsSummary <- data.frame(Network = 1:length(results$clusters), 
                               Size = sapply(results$clusters, length), pvalue = results$pvalues, 
                               topGOBPTerm = topGOTerms)
  
  write.table(resultsSummary, summaryTable, sep = "\t", row.names = FALSE, quote = FALSE)
  networks <- getVisSubnetworks(output, results$clusters, interactome)
  saveRDS(networks, file = visNetworks)
  
}

opt <- opt[names(opt) != "help"]
do.call(GigaSubnetworks, opt)