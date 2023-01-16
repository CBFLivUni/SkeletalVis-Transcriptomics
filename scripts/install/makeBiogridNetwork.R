#make a protein-protein interaction network using biogrid
library(readr)
library(biomaRt)
library(igraph)
library(org.Hs.eg.db)

#set options to avoid timeout when download large files
httr::set_config(httr::config(ssl_verifypeer = FALSE))
options(timeout=100000)

#specify the stringdb network version to download
version <- "4.4.200"
url <- sprintf("https://downloads.thebiogrid.org/Download/BioGRID/Release-Archive/BIOGRID-%s/BIOGRID-ALL-%s.mitab.zip",version,version)

tmp <- tempfile(fileext = ".zip")
print(tmp)
download.file(url,tmp)


#read in the stringdb network
network <- read_delim(unzip(tmp),delim = "\t", escape_double = FALSE, trim_ws = TRUE)

#keep only the human data
network <- network[ network$`Taxid Interactor A` == "taxid:9606" & network$`Taxid Interactor B` == "taxid:9606", ]

#get just the entrez IDs
network <- na.omit(network[,1:2])
network$`ID Interactor B` <- sapply(strsplit(network$`ID Interactor B`,split = ":"),"[[",2)
network$`#ID Interactor A` <- sapply(strsplit(network$`#ID Interactor A`,split = ":"),"[[",2)

#convert the entrez ids to gene symbols
network$symbol1 <- stack(mget(as.character(network$`#ID Interactor A`), org.Hs.egSYMBOL, ifnotfound = NA))[,1]
network$symbol2 <- stack(mget(as.character(network$`ID Interactor B`), org.Hs.egSYMBOL, ifnotfound = NA))[,1]

#make the igraph network
network<-graph.data.frame(na.omit(network[,3:4]),directed = F)
network<-decompose.graph(network)[[which.max(sapply(decompose.graph(network),vcount))]]
network<-simplify(network)

#write to file
humanBioGridTable<-as.data.frame(get.edgelist(network))
write.table(humanBioGridTable, file="networks/humanBioGridNetwork.txt", sep = "\t",row.names=F)


#map the human genes to other species
mapGenes <- function(genes,homology){
  genes <- homology[ na.omit(match(genes,homology[,4])),3]
  return(genes)
}

#convert the human network to other species using gene homology
makeNonHumanNework <- function(species,network){
  
  homology<-read.delim(sprintf("homology/Homology.%s.txt",species),stringsAsFactors = F)
  nonHumanNodes<-sapply(V(network)$name,mapGenes,homology)
  nonHumanNodes[sapply(nonHumanNodes, function(x) length(x)==0)] <- NA
  V(network)$nonHuman<-unlist(nonHumanNodes)
  
  #keep the largest connect component
  network.nonHuman<-delete.vertices(network,V(network)$name[is.na(V(network)$nonHuman)])
  V(network.nonHuman)$name<-V(network.nonHuman)$nonHuman  
  network.nonHuman<-decompose.graph(network.nonHuman)[[1]]
  nonHumanBioGridNetwork<-network.nonHuman
  
  nonHumanBioGridTable<-as.data.frame(get.edgelist(nonHumanBioGridNetwork))
  write.table(nonHumanBioGridTable, file=sprintf("networks/%sBioGridNetwork.txt",species), sep = "\t",row.names=F)
  
}

speciesList <- c("mouse","rat","cow","horse","zebrafish","pig")
lapply(speciesList,makeNonHumanNework,network)