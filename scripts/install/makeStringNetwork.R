#make a protein-protein interaction network using stringDB
library(readr)
library(biomaRt)
library(igraph)
httr::set_config(httr::config(ssl_verifypeer = FALSE))

#specify the stringdb network version
version <- 11.5
url <- sprintf("https://stringdb-static.org/download/protein.links.detailed.v%s/9606.protein.links.detailed.v%s.txt.gz",version,version)

#download and parse the network
tmp <- tempfile()
download.file(url,tmp)

network <- read_delim(gzfile(tmp), 
                      " ", escape_double = FALSE, trim_ws = TRUE)

network <- as.data.frame(apply(network,2,function(x) gsub('\\s+', '',x)))

#keep interactions with experimental or database support
network <- network[ apply(network,1,function(x) sum(as.numeric(x[7:8]))>0),]

#map the ensemblIDs to gene names
network$protein1 <- gsub(pattern = "9606.",x=network$protein1,replacement = "")
network$protein2 <- gsub(pattern = "9606.",x=network$protein2,replacement = "")

mart <- useDataset(dataset = "hsapiens_gene_ensembl",mart = useMart("ensembl"))

genes <- unique(c(network$protein1,network$protein2))

G_list <- getBM(filters= "ensembl_peptide_id", attributes= c("ensembl_peptide_id","hgnc_symbol"),values=genes,mart= mart)
network <- merge(network,G_list,by.x="protein1",by.y="ensembl_peptide_id")
network <- merge(network,G_list,by.x="protein2",by.y="ensembl_peptide_id")


network <- network[,c(1,2,7,8,11,12)]
colnames(network)[5:6] <- c("gene1","gene2")

network$experimental <- as.numeric(network$experimental)
network$database <- as.numeric(network$database)
network$score <- network$experimental + network$database 

#make a igraph network
network <- graph.data.frame(network[,5:7],directed = F)
network.simple <- simplify(network,edge.attr.comb = list(score="max"))
network.strict <- delete.edges(network.simple,E(network.simple)[E(network.simple)$score<400])
network.strict <- decompose.graph(network.strict)[[which.max(sapply(decompose.graph(network.strict),vcount))]]

#write out as table for future use
humanStringTable <- as.data.frame(get.edgelist(network.strict))
humanStringTable$score <- E(network.strict)$score
write.table(humanStringTable, file="networks/humanStringNetwork.txt", sep = "\t",row.names=F)



#map the human genes to other species
mapGenes <- function(genes,homology){
  genes <- homology[ na.omit(match(genes,homology[,4])),3]
  return(genes)
}

#convert the human network to other species using gene homology
makeNonHumanNework <- function(species,network){
  
  homologyFile <- sprintf("homology/Homology.%s.txt",species)
  homology <- read.delim(homologyFile,stringsAsFactors = F)
  nonHumanNodes <- sapply(V(network)$name,mapGenes,homology)
  nonHumanNodes[sapply(nonHumanNodes, function(x) length(x)==0)]  <-  NA
  V(network)$mapped <- unlist(nonHumanNodes)
  
  network.nonHuman <- delete.vertices(network,V(network)$name[is.na(V(network)$mapped)])
  V(network.nonHuman)$name <- V(network.nonHuman)$mapped  
  network.nonHuman <- decompose.graph(network.nonHuman)[[1]]
  nonHumanStringNetwork <- network.nonHuman
  
  nonHumanStringTable <- as.data.frame(get.edgelist(nonHumanStringNetwork))
  nonHumanStringTable$score <- E(nonHumanStringNetwork)$score
  write.table(nonHumanStringTable, file=sprintf("networks/%sStringNetwork.txt",species), sep = "\t",row.names=F)
  
}

speciesList <- c("mouse","rat","cow","horse","zebrafish","pig")
lapply(speciesList,makeNonHumanNework,network.strict)


