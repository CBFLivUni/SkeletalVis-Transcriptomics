#!/bin/sh

#make the subdirectories
mkdir transcriptomes
mkdir homology
mkdir networks
mkdir enrichment

#download the fasta and gtf files needed
bash install/downloadReferences.sh

#make all the supporting files needed for the pipeline
Rscript install/makeHomologyTables.R
Rscript install/makeBiogridNetwork.R
Rscript install/makeStringNetwork.R
Rscript install/makeReactomePathways.R
Rscript install/makeTFLibrary.R
Rscript install/getFastqScreenGenomes.sh

