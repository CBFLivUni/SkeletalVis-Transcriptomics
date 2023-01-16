#make maps of pathways to genes for each species

library("reactome.db")
library("org.Hs.eg.db")
library("dplyr")
library("tidyr")
library("org.Mm.eg.db")
library("org.Bt.eg.db")
library("org.Dr.eg.db")
library("org.Ss.eg.db")
library("org.Rn.eg.db")



getReactomePathways <- function(species){
  
  
  #the reactome pathways are only provided with entrez IDs
  #so we need to convert them to gene symbols to match the data we have
  
  symbol2eg <- dplyr::case_when(
   species == "Human" ~ "org.Hs.egSYMBOL",
   species == "Mouse" ~ "org.Mm.egSYMBOL",
   species == "Pig" ~ "org.Ss.egSYMBOL",
   species == "Cow" ~ "org.Bt.egSYMBOL",
   species == "Rat" ~ "org.Rn.egSYMBOL",
  )
  
  symbol2eg <- as.list(get(symbol2eg))

  #get eg to reactome pathway
  reactome2eg <- as.list(reactomePATHID2EXTID)
  
  speciesID <- dplyr::case_when(
    species == "Human" ~ "R-HSA",
    species == "Mouse" ~ "R-MMU",
    species == "Pig" ~ "R-SSC",
    species == "Cow" ~ "R-BTA",
    species == "Rat" ~"R-RNO",
  )
  
  #filter to just species pathways
  reactome2eg <- reactome2eg[grep(speciesID,names(reactome2eg))]
  
  #function to search 
  grepREACTOME <- function(id,mapkeys){
    unique(unlist(mapkeys[id],use.names=FALSE))
  }
  
  #convert the entrez ids to gene symbols for each pathway
  reactome <- lapply(reactome2eg,grepREACTOME,symbol2eg)
  
  
  #get the pathway names rather than the ids
  reactome2name <- as.list(reactomePATHID2NAME)
  reactomeNames <- sapply(names(reactome),grepREACTOME,reactome2name)
  names(reactome) <- reactomeNames
  
  saveRDS(reactome,file=sprintf("enrichment/reactomePathways%s.RDS",species))
  
  
}

speciesList <- c("Human","Mouse","Rat","Cow","Pig")
lapply(speciesList,getReactomePathways)

