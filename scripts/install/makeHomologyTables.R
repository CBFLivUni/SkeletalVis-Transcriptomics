library("biomaRt")

httr::set_config(httr::config(ssl_verifypeer = FALSE))
options(timeout = 120)
#use ensembl to create a human to other species gene name map
makeHomologyTable <- function(martName,speciesName,version=103){

  mart1 = useEnsembl("ensembl", dataset="hsapiens_gene_ensembl",version = version)
  mart2 = useEnsembl("ensembl", dataset=paste0(martName,"_gene_ensembl"),version = version) 
  
  mapping <- getLDS(attributes=c("ensembl_gene_id","external_gene_name"), mart=mart1,
         attributesL=c("ensembl_gene_id","external_gene_name"), martL=mart2)
  
  mapping <- mapping[,c(3,1,4,2)]
  colnames(mapping)<-c(paste0(speciesName,"ID"),"HumanID",paste0(speciesName,"GeneSymbol"),"HumanGeneSymbol")
  mapping[ mapping[,3]=="",3] <- NA
  write.table(mapping, file=sprintf("homology/Homology.%s.txt",speciesName),sep = "\t",row.names=F,quote = F)

}

#species to perform mapping for
martNames<-c("mmusculus","rnorvegicus","ecaballus","btaurus","sscrofa","drerio")
species<-c("mouse","rat","horse","cow","pig","zebrafish")

mapply(makeHomologyTable,martNames,species)
