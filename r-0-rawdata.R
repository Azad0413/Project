rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/BJTC-385-10/")
if (! dir.exists("./00_rawdata")){
  dir.create("./00_rawdata")
}
setwd("./00_rawdata")
library(GEOquery)
library(Biobase)
library(GEOquery)
library(Biobase)
## GSE110811----
gset<-getGEO("GSE110811",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL16686",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)    
## 将GB_ACC换成Symbol
library(org.Hs.eg.db)
library(clusterProfiler)
gene_transform<-bitr(gpl$GB_ACC,fromType = "REFSEQ",
                     toType = "SYMBOL",
                     OrgDb = org.Hs.eg.db)
probe2symbol <- data.frame(ID=gpl$ID,REFSEQ=gpl$GB_ACC)
probe2symbol <- merge(probe2symbol,gene_transform,by='REFSEQ')
probe2symbol <- probe2symbol[,-1]
colnames(probe2symbol) <- c('ID','symbol')

dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symbol$ID<-as.character(probe2symbol$ID)
library(tidyverse)
dat<-dat %>%
  inner_join(probe2symbol,by='ID')%>% 
  dplyr::select(-ID)%>%     ## 去除多余信息
  dplyr::select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
pd<-pData(a)
##31个RB样本和3个正常视网膜样本。
group<-data.frame(sample=pd$geo_accession,group=pd$`tissue:ch1`)
table(group$group)
group$group<-ifelse(group$group=='Normal retina','control','RB')
group<-group[order(group$group),]
dat<-dat[,group$sample]
write.table(dat,file = 'dat(GSE110811).xls',sep = '\t',row.names = T,quote = F)
write.table(group,file = 'group(GSE110811).xls',sep = '\t',row.names = F,quote = F)

phenotype

## GSE24673----
gset<-getGEO("GSE24673",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
a=gset[[1]]
library(AnnoProbe)
probe2symbol <- idmap(gpl = 'GPL6244',type = 'bioc')
colnames(probe2symbol) <- c('ID','symbol')
dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symbol$ID<-as.character(probe2symbol$ID)
dat<-dat %>%
  inner_join(probe2symbol,by='ID')%>% 
  dplyr::select(-ID)%>%     ## 去除多余信息
  dplyr::select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
pd<-pData(a)
## 敏感和耐药
group<-data.frame(sample=pd$geo_accession,group=pd$`cell type:ch1`)
table(group$group)
# group<-subset(group,group=='pCR'|group=='RD')
group$group<-ifelse(group$group=='primary cultured cells','control','RB')
group<-group[order(group$group),]
dat<-dat[,group$sample]
write.table(dat,file = 'dat(GSE24673).xls',sep = '\t',row.names = T,quote = F)
write.table(group,file = 'group(GSE24673).xls',sep = '\t',row.names = F,quote = F)











## GSE111168----
gset<-getGEO("GSE111168",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-read.delim2('GSE111168_gene_exp.txt.gz')
expr <- expr[,c(1,5:10)]
dat <- expr%>%column_to_rownames(var = 'gene')
write.table(dat,file = 'dat(GSE111168).xls',sep = '\t',row.names = T,quote = F)
group <- data.frame(sample=colnames(dat),group=c(rep('control',3),rep('RB',3)))
write.table(group,file = 'group(GSE111168).xls',sep = '\t',row.names = F,quote = F)
