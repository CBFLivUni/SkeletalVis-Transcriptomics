#download the chea3 transcription factor libraries and format as a list
libs <- c("ARCHS4_Coexpression.gmt",
          "ENCODE_ChIP-seq.gmt",
          "Enrichr_Queries.gmt",
          "GTEx_Coexpression.gmt",
          "Literature_ChIP-seq.gmt",
          "ReMap_ChIP-seq.gmt")

downloadLib <- function(lib){
  url <- sprintf("https://maayanlab.cloud/chea3/assets/tflibs/%s",lib)
  download.file(url,lib)
  lib <- readLines(lib)
  lib <- strsplit(lib,"\t")
  names(lib) <- sapply(lib,"[[",1)
  lib <-  lapply(lib,"[",-1)
}

libsTFs <- lapply(libs,downloadLib)
names(libsTFs) <- libs

saveRDS(libsTFs,"enrichment/TFs.RDS")