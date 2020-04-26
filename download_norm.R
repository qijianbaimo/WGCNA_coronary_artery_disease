
###################
#
#  GSE20686_RAW and filelist.txt were downloaded from NCBI "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE20nnn/GSE20686/suppl/"
#
#  GSE20680 and GSE20681 were the SubSeries of GSE20686_RAW
#
###################

options(stringsAsFactors = F)

library(GEOquery)

## download filelist.txt
system('wget https://ftp.ncbi.nlm.nih.gov/geo/series/GSE20nnn/GSE20686/suppl/filelist.txt')

## download raw data
geo <- getGEOSuppFiles('GSE20686', fetch_files = TRUE, makeDirectory = F)

dir.create('GSE20686_RAW')
system('tar xvf GSE20686_RAW.tar -C ./GSE20686_RAW/')
system('gunzip ./GSE20686_RAW/*')

Group<-read.csv("Clinic.txt",sep="\t",header = T)
Group$Sample_name<-paste(Group$Sample,Group$Point2,sep=" - ")
file_name<-read.csv("filelist.txt",skip = 3,header = F,sep="\t")
file_name$Sample<-gsub("_.*","",file_name$V2)
For_nor<-merge(Group,file_name[,c(2,6)],by="Sample")
For_nor$V2<-paste("GSE20686_RAW/",For_nor$V2,sep="")
rownames(For_nor)<-For_nor$Sample
For_nor<-For_nor[,c(5,3,3)]
colnames(For_nor)<-c("FileName","Treatment","GErep")
For_nor$FileName<-gsub("[.]gz","",For_nor$FileName)

#install.packages('https://www.bioconductor.org/packages/2.13/bioc/src/contrib/Agi4x44PreProcess_1.22.0.tar.gz', repos = NULL)
#BiocManager::install('hgug4112a.db')

library(Agi4x44PreProcess)
library(hgug4112a.db)

targets <-For_nor[Group$Sample,]
dd=read.AgilentFE(targets,makePLOT=FALSE)

ddNORM=BGandNorm(dd,BGmethod="half",NORMmethod="quantile",
                 foreground="MeanSignal",background="BGMedianSignal",
                 offset=50,makePLOTpre=FALSE,makePLOTpost=FALSE)
ddFILT=filter.probes(ddNORM,
                     control=TRUE,
                     wellaboveBG=TRUE,
                     isfound=TRUE,
                     wellaboveNEG=TRUE,
                     sat=TRUE,
                     PopnOL=TRUE,
                     NonUnifOL=T,
                     nas=TRUE,
                     limWellAbove=75,
                     limISF=75,
                     limNEG=75,
                     limSAT=75,
                     limPopnOL=75,
                     limNonUnifOL=75,
                     limNAS=100,
                     makePLOT=F,annotation.package="hgug4112a.db",flag.counts=T,targets)

system('gzip RawDataNOCtrl.txt')

sessionInfo()
