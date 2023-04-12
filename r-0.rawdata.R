rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/BJTC-321")
if (! dir.exists("./00_rawdata")){
  dir.create("./00_rawdata")
}
setwd("./00_rawdata")
library(GEOquery)
library(Biobase)
library(tidyr)
library(AnnoProbe)
library(tidyverse)
library(readxl)
## 01-1 GSE67472------
gset<-getGEO("GSE67472",destdir = '.',GSEMatrix = T,getGPL = F)
a=gset[[1]]
pd<-pData(a)
expr<-as.data.frame(exprs(gset[[1]]))

gpl<-getGEO("GPL16311",destdir = '.')
gpl<-Table(gpl) 
colnames(gpl)
probe2entrizd<-gpl %>%
  select('ID','SPOT_ID')%>%
  filter('SPOT_ID'!='')%>%
  separate('SPOT_ID',c('ENTRIZD','drop'),sep = '///')%>%
  select(-drop)
#probe2symbol=probe2symbol[probe2symbol$symbol!='',]
colnames(probe2entrizd)<-c('ID','ENTRIZD')
dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2entrizd$ID<-as.character(probe2entrizd$ID)
dat<-dat %>%
  inner_join(probe2entrizd,by='ID')%>% 
  select(-ID)%>%     ## 去除多余信息
  select(ENTRIZD,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(ENTRIZD,.keep_all = T)%>%      ## symbol留下第一个
  select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
 
## entrizd转symbol
# 将Gene_ID更换为Symbol NCBI的Gene_id,就是Entrezid。
library(org.Hs.eg.db)
library(clusterProfiler)
geneid<-rownames(dat)
id2symbol<-bitr(geneid,fromType = "ENTREZID",
                     toType = "SYMBOL",
                     OrgDb = org.Hs.eg.db)
colnames(id2symbol)=c('Entrezid','Symbol')
dat=dat[id2symbol$Entrezid,]
rownames(dat)=id2symbol$Symbol

group<-data.frame(sample=pd$geo_accession,group=pd$`disease state:ch1`)
group$group<-ifelse(group$group=='healthy','control','asthma')
group<-group[order(group$group),]
dat.final<-dat[,group$sample]
write.table(dat.final,file = 'dat.final.xls',quote = F,sep = '\t',row.names = T)
write.table(group,file = 'group.xls',quote = F,sep = '\t',row.names = F)
clinical<-data.frame(sample=pd$geo_accession,group=pd$`disease state:ch1`,age=pd$`age:ch1`,gender=pd$`gender:ch1`)
clinical$group<-ifelse(clinical$group=='healthy','control','asthma')
clinical<-clinical[order(clinical$group),]
write.table(clinical,file = 'clinical.xls',sep = '\t',row.names = F,quote = F)
