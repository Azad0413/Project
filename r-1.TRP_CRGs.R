rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/YQ444-8/")
if (! dir.exists("./01_TRP_CRGs")){
  dir.create("./01_TRP_CRGs")
}
setwd("./01_TRP_CRGs")
library(readxl)
msigdb<-read_xlsx('MsigDB.xlsx')
kegg<-read_xlsx('KEGG.xlsx')
genecards<-read.csv('GeneCards.csv')
genecards<-genecards[which(genecards$Relevance.score>2),]
genecards<-data.frame(Symbol=genecards$Gene.Symbol)
TRP<-rbind(msigdb,kegg)
TRP<-rbind(TRP,genecards)
TRP<-TRP[!duplicated(TRP$Symbol),]
##468
write.table(TRP,file = 'TRP_CRGs.xls',sep = '\t',row.names = F,quote = F)


