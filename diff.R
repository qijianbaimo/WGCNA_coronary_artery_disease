
### Data Preprocessing

options(stringsAsFactors = FALSE);

nor_data=read.table("RawDataNOCtrl.txt.gz",sep="\t",header = T,stringsAsFactors = F)
#nor_data

my.median<-function(x,cat){
  tapply(x,cat,median)
}
nor_data2<-apply(nor_data[,-1:-7],2,my.median,cat=nor_data$SYMBOL)
nor_data3=data.frame(Sample=colnames(nor_data2),t(nor_data2))

a=read.table("RawDataNOCtrl.txt.gz",sep="\t",header = F,stringsAsFactors = F)[1,]

SaGo=t(a[-1:-7])
#SaGo
SG=t(as.data.frame(strsplit(SaGo,split = " - ")))

SG=data.frame(Sample=c("P1",SG[,1]), Group=c("P2",SG[,2]))

nor.data=cbind(SG[-1,],nor_data3)


#Group=read.table("../Clinic",sep="\t",header = T,row.names = 1)
#Group2=as.data.frame(t(Group[,-1]))
#colnames(Group2)=row.names(Group)
#For_ana=rbind(Group2,nor.data)
For_ana=nor.data[,-3]
n=ncol(For_ana)

k=data.frame(Average1=c(1),Average0=c(1),FC=c(1),FClog2=c(1),absFC=c(1),AUC=c(1),p_value=c(1))

library(pROC)
i=5891


for(i in 3:n){
  a=For_ana[For_ana[,2]=="1",i]
  Average1=log2(sum(2^a)/length(a))
  b=For_ana[For_ana[,2]=="0",i]
  Average2=log2(sum(2^b)/length(b))
  FC_log2=Average1-Average2
  FC=2^FC_log2
  if (FC_log2 < 0){
    absFC=-1/FC
  }else{
    absFC=FC
  }
#  pvalue=wilcox.test(a,b)$p.value####
  for_ROC=For_ana[,c(2,i)]
  pvalue=t.test(a,b)$'p.value'
  auc=auc(for_ROC$Group, for_ROC[,2])
  
  k[i-2,]<-c(Average1,Average2,FC,FC_log2,absFC,auc,pvalue)
  
}
#shapiro.test(a)
anno<- data.frame(t(nor_data[,1:7]))
anno<-nor_data[!duplicated(nor_data[,1]),1:6]


kk=data.frame(SYMBOL=colnames(For_ana[,-1:-2]),k)
final<-merge(anno,kk,by="SYMBOL")

write.table(final,"All_compare_gene.txt",sep="\t",quote = F ,row.names = F)

#write.csv(final)

write.table(For_ana,"for_WGCNA_gene.txt",sep="\t",quote = F ,row.names = F)

