rm(list = ls())
# 02 差异分析---------
setwd("/data/nas1/luchunlin/project/JNZK-218-8/")
if (! dir.exists("./06_PPI")){
  dir.create("./06_PPI")
}
setwd("./06_PPI")
MCC<-read.delim2('../06_PPI/PPI.tsv_MCC_top30 default node.csv',sep = ',')
MNC<-read.delim2('../06_PPI/PPI.tsv_MNC_top30 default node.csv',sep = ',')
Degree<-read.delim2('../06_PPI/PPI.tsv_Degree_top30 default node.csv',sep = ',')
Closeness<-read.delim2('../06_PPI/PPI.tsv_Closeness_top30 default node.csv',sep = ',')
Radiality<-read.delim2('../06_PPI/PPI.tsv_Radiality_top30 default node.csv',sep = ',')
Stress<-read.delim2('../06_PPI/PPI.tsv_Stress_top30 default node.csv',sep = ',')
EPC<-read.delim2('../06_PPI/PPI.tsv_EPC_top30 default node.csv',sep = ',')
hubgene<-MCC$name[which(MCC$name%in%MNC$name,)]%>%as.data.frame()
hubgene<-hubgene[hubgene$.%in%Degree$name,]%>%as.data.frame()
hubgene<-hubgene[hubgene$.%in%Closeness$name,]%>%as.data.frame()
hubgene<-hubgene[hubgene$.%in%Radiality$name,]%>%as.data.frame()
hubgene<-hubgene[hubgene$.%in%Stress$name,]%>%as.data.frame()
hubgene<-hubgene[hubgene$.%in%EPC$name,]%>%as.data.frame()
colnames(hubgene)<-'symbol'
write.table(hubgene,file = 'hubgene.xls',row.names = F,quote = F,sep = '\t')
