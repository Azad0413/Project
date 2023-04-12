rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/BJTC-327")
if (! dir.exists("./00_rawdata")){
  dir.create("./00_rawdata")
}
setwd("./00_rawdata")
## 使用客户提供的测序数据，对数据进行预处理
## 3个AIH样本、3个PBC样本、3个HBV样本。
library(lance)
library(tidyverse)
library(readxl)
dat.fpkm<-read_xlsx('dat.fpkm.xlsx')
dat.fpkm<-dat.fpkm%>%
  mutate(rowMean=rowMeans(.[grep('D',names(.))]))%>%
  arrange(desc(rowMean))%>%
  distinct(gene_name,.keep_all = T)%>%
  select(-rowMean)%>%
  tibble::column_to_rownames(colnames(.)[1])
colnames(dat.fpkm)<-gsub('_FPKM','',colnames(dat.fpkm),fixed = T)

dat.count<-read_xlsx('dat.counts.xlsx')
dat.count<-dat.count%>%
  mutate(rowMean=rowMeans(.[grep('D',names(.))]))%>%
  arrange(desc(rowMean))%>%
  distinct(gene_name,.keep_all = T)%>%
  select(-rowMean)%>%
  tibble::column_to_rownames(colnames(.)[1]) 
colnames(dat.count)<-gsub('_count','',colnames(dat.count),fixed = T)
## 分成AIH+PBC  AIH+HBV
AIH.sample<-c('D192100008','D192100018','D192100025')
PBC.sample<-c('D192100002','D192100020','D192100024')
HBV.sample<-c('D192100026','D192100027','D192100033')
## 
AIH.PBC.count<-cbind(dat.count[,AIH.sample],dat.count[,PBC.sample])
AIH.PBC.fpkm<-cbind(dat.fpkm[,AIH.sample],dat.fpkm[,PBC.sample])
AIH.HBV.count<-cbind(dat.count[,AIH.sample],dat.count[,HBV.sample])
AIH.HBV.fpkm<-cbind(dat.fpkm[,AIH.sample],dat.fpkm[,HBV.sample])
write.table(AIH.PBC.count,file = 'AIH.PBC.count.xls',sep = '\t',row.names = T,quote = F)
write.table(AIH.HBV.count,file = 'AIH.HBV.count.xls',sep = '\t',row.names = T,quote = F)
write.table(AIH.PBC.fpkm,file = 'AIH.PBC.fpkm.xls',sep = '\t',row.names = T,quote = F)
write.table(AIH.HBV.fpkm,file = 'AIH.HBV.fpkm.xls',sep = '\t',row.names = T,quote = F)
