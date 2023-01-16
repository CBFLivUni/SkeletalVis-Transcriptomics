#!/bin/sh

#download the ensembl fasta and gtf files for the specified species and the release
#concatenates the cdna and the ncrna 
RELEASE=103
SPECIES_NAME=(human mouse rat horse cow zebrafish pig)
SPECIES_LIST=(homo_sapiens mus_musculus rattus_norvegicus equus_caballus bos_taurus danio_rerio sus_scrofa)


for ((i=0;i<${#SPECIES_LIST[@]};++i)); do

#get the cnda
wget "ftp://ftp.ensembl.org/pub/release-$RELEASE/fasta/"${SPECIES_LIST[i]}"/cdna/${SPECIES_LIST[i]^}*.all.*"

#get the ncrna
wget "ftp://ftp.ensembl.org/pub/release-$RELEASE/fasta/"${SPECIES_LIST[i]}"/ncrna/${SPECIES_LIST[i]^}*"

#concentate the two
cat `ls ${SPECIES_LIST[i]^}*.fa.gz` > ${SPECIES_NAME[i]}_release${RELEASE}.fa.gz
rm ${SPECIES_LIST[i]^}*


#get the correct gtf file
wget "ftp://ftp.ensembl.org/pub/release-$RELEASE/gtf/"${SPECIES_LIST[i]}"/${SPECIES_LIST[i]^}*${RELEASE}.gtf.gz"
mv `ls ${SPECIES_LIST[i]^}*gtf* | grep -v "abinitio" | grep -v "chr"` transcriptomes/${SPECIES_NAME[i]}_release${RELEASE}.gtf.gz

#ungzip
gzip -d < transcriptomes/${SPECIES_NAME[i]}_release${RELEASE}.gtf.gz > transcriptomes/${SPECIES_NAME[i]}_release${RELEASE}.gtf

done

zcat human_release103.fa.gz  | awk '/^>/ {P=index($0,"CHR_HS")==0} {if(P) print} ' | gzip > human_release103_nohaplo.fa.gz
