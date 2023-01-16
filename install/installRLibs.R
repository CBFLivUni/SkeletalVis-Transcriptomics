options(Ncpus = 7)

#function to install libararies from bioconductor
installPackage <- function(libName) {
  if(libName %in% rownames(installed.packages()) == FALSE){
    BiocManager::install(libName,ask = FALSE)
  }}

#Install R libraries need for pipeline
install.packages(c("curl", "httr"))
install.packages("devtools")
install.packages("intergraph")

require(devtools)

if (!requireNamespace("BiocManager"))
  install.packages("BiocManager")
BiocManager::install()
#avoid openblas phread22 errors
BiocManager::install("preprocessCore", configure.args="--disable-threading",force=TRUE)

install_github("cstubben/ENAbrowseR")
install_github("briatte/ggnet")
install_github("statOmics/zingeR")
install_github("cran/GeoDE")




#read the libraries needed
packagesToInstall <- read.delim("install/packagesToInstall.txt",header=F,stringsAsFactors = F)

#install all the libraries
sapply(packagesToInstall[,1],installPackage)

install_github("grimbough/arrayQualityMetrics")


#make and install the databases and files needed for the pipeline
source("install/makeMicroarrayAnnotations.R")
source("install/makeEnsemblPackages.R")