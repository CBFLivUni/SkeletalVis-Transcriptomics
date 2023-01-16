#create DB annotation packages for those microarray chips where annotation libraries are not available from Bioconductor
library(RSQLite)
library(AnnotationForge)
library(human.db0)
library(mouse.db0)
library(rat.db0)
library(GEOquery)

options(timeout = 100000)

#download primeview annotations
tmp <- "install/PrimeView.na36.annot.zip"
tmp <- unzip(tmp)[1]

#make and install package
makeDBPackage("HUMANCHIP_DB",affy=TRUE,prefix="primeview",fileName=tmp,baseMapType="eg",outputDir = ".",
              version="1.0.0",manufacturer = "Affymetrix",chipName = "Human Prime View Array",manufacturerUrl = "http://www.affymetrix.com")
install.packages("primeview.db",repos = NULL)


# download.file("http://www.affymetrix.com/Auth/analysis/downloads/na36/ivt/HT_MG-430_PM.na36.annot.csv.zip",tmp)
# tmp <- unzip(tmp)[1]
GPL11180 <- getGEO('GPL11180', destdir=".")
GPL11180 <- Table(GPL11180)[,c("ID","Entrez Gene")]
GPL11180 <- GPL11180[ !is.na(GPL11180$ID),]
write.table(GPL11180,file = "GPL11180.txt",row.names = F,quote = F,sep="\t")

#make and install package
makeDBPackage("MOUSECHIP_DB",affy=F,prefix="htmg430pm",fileName="GPL11180.txt",baseMapType="eg",outputDir = ".",
              version="1.0.0",manufacturer = "Affymetrix",chipName = "htmg430pm")
install.packages("htmg430pm.db",repos = NULL)

#retrieve and extract the probe ids for GPL7202
GPL7202 <- getGEO('GPL7202', destdir=".")
GPL7202 <- Table(GPL7202)[,c("SPOT_ID","GENE"),]
GPL7202 <- GPL7202[ !is.na(GPL7202$SPOT_ID),]
write.table(GPL7202,file = "GPL7202.txt",row.names = F,quote = F,sep="\t")

#make and install package
makeDBPackage("MOUSECHIP_DB",affy=F,prefix="AgilentMouse014868",fileName="GPL7202.txt",baseMapType="eg",outputDir = ".",
              version="1.0.0",manufacturer = "Agilent",chipName = "AgilentMouse014868")
install.packages("AgilentMouse014868.db",repos = NULL)


GPL10787 <- getGEO('GPL10787', destdir=".")
GPL10787 <- Table(GPL10787)[,c("SPOT_ID","GENE")]
GPL10787 <- GPL10787[ !is.na(GPL10787$SPOT_ID),]
write.table(GPL10787,file = "GPL10787.txt",row.names = F,quote = F,sep="\t")

makeDBPackage("MOUSECHIP_DB",affy=F,prefix="AgilentMouse028005",fileName="GPL10787.txt",baseMapType="eg",outputDir = ".",
              version="1.0.0",manufacturer = "Agilent",chipName = "AgilentMouse028005")
install.packages("AgilentMouse028005.db",repos = NULL)


GPL7294 <- getGEO('GPL7294', destdir=".")
GPL7294 <- Table(GPL7294)[,c("SPOT_ID","GENE")]
GPL7294 <- GPL7294[ !is.na(GPL7294$SPOT_ID),]
write.table(GPL7294,file = "GPL7294.txt",row.names = F,quote = F,sep="\t")
makeDBPackage("RATCHIP_DB",affy=F,prefix="AgilentRat014879",fileName="GPL7294.txt",baseMapType="eg",outputDir = ".",
              version="1.0.0",manufacturer = "Agilent",chipName = "AgilentRat014879")
install.packages("AgilentRat014879.db",repos = NULL)

  GPL14797 <- getGEO('GPL14797', destdir=".")
GPL14797 <- Table(GPL14797)[,c("SPOT_ID","REFSEQ")]
GPL14797 <- GPL14797[ !is.na(GPL14797$SPOT_ID),]
write.table(GPL14797,file = "GPL14797.txt",row.names = F,quote = F,sep="\t")

makeDBPackage("RATCHIP_DB",affy=F,prefix="AgilentRat028279",fileName="GPL14797.txt",baseMapType="refseq",outputDir = ".",
              version="1.0.0",manufacturer = "Agilent",chipName = "	AgilentRat028279")

install.packages("AgilentRat028279.db",repos = NULL)



GPL17077 <- getGEO('GPL17077', destdir=".")
GPL17077 <- Table(GPL17077)[,c("SPOT_ID","REFSEQ")]
GPL17077 <- GPL17077[ !is.na(GPL17077$SPOT_ID),]
write.table(GPL17077,file = "GPL17077.txt",row.names = F,quote = F,sep="\t")

makeDBPackage("HUMANCHIP_DB",affy=F,prefix="AgilentHuman039494",fileName="GPL17077.txt",baseMapType="refseq",outputDir = ".",
              version="1.0.0",manufacturer = "Agilent",chipName = "AgilentHuman039494")
install.packages("AgilentHuman039494.db",repos = NULL)


GPL21163 <- getGEO('GPL21163', destdir=".")
GPL21163 <- Table(GPL21163)[,c("ID","REFSEQ")]
GPL21163 <- GPL21163[ !is.na(GPL21163$ID),]
write.table(GPL21163,file = "GPL21163.txt",row.names = F,quote = F,sep="\t")

makeDBPackage("MOUSECHIP_DB",affy=F,prefix="AgilentMouse074809",fileName="GPL21163.txt",baseMapType="refseq",outputDir = ".",
              version="1.0.0",manufacturer = "Agilent",chipName = "AgilentMouse074809")
install.packages("AgilentMouse074809.db",repos = NULL)



GPL14550 <- getGEO('GPL14550', destdir=".")
GPL14550 <- Table(GPL14550)[,c("ID","REFSEQ")]
GPL14550 <- GPL14550[ !is.na(GPL14550$ID),]
write.table(GPL14550,file = "GPL14550.txt",row.names = F,quote = F,sep="\t")

makeDBPackage("HUMANCHIP_DB",affy=F,prefix="AgilentHuman028004",fileName="GPL14550.txt",baseMapType="refseq",outputDir = ".",
              version="1.0.0",manufacturer = "Agilent",chipName = "AgilentHuman028004")
install.packages("AgilentHuman028004.db",repos = NULL)


GPL6480 <- getGEO('GPL6480', destdir=".")
GPL6480 <- Table(GPL6480)[,c("SPOT_ID","REFSEQ")]
GPL6480 <- GPL6480[ !is.na(GPL6480$SPOT_ID),]
write.table(GPL6480,file = "GPL6480.txt",row.names = F,quote = F,sep="\t")

makeDBPackage("HUMANCHIP_DB",affy=F,prefix="AgilentHuman014850",fileName="GPL6480.txt",baseMapType="refseq",outputDir = ".",
              version="1.0.0",manufacturer = "Agilent",chipName = "AgilentHuman014850")
install.packages("AgilentHuman014850.db",repos = NULL)