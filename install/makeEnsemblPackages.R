library("ensembldb")
library("AnnotationHub")
library("AnnotationForge")


#create the ensembldb packages
# Load the annotation resource.
ah <- AnnotationHub()

makeEnsemblDBPackage <- function(speciesName,version=103){
    
  ahDb <- query(ah, pattern = c(speciesName, "EnsDb",version))[[1]]
  
  package <- makeEnsembldbPackage(ensdb = dbfile(dbconn(ahDb)), version = "0.0.1",
                       maintainer = "Jamie Soul <jamie.soul@newcastle.ac.uk>",
                       author = "Jamie Soul")
  return(package)
}

#make a ensembldb for each of the designated species
species <- c("Homo sapiens","Mus musculus","Rattus norvegicus","Equus caballus","Bos taurus","Sus scrofa","Danio rerio")
packages <- sapply(species,makeEnsemblDBPackage)
install.packages(list.files(pattern = "EnsDb*"),repos = NULL)

#tidy up
removeCache(ah,ask=FALSE)
unlink("EnsDb*", recursive = TRUE)


