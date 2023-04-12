rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/BJTC-420-1/")
if (! dir.exists("./00_rawdata")){
  dir.create("./00_rawdata")
}
setwd("./00_rawdata")
library(tidyverse)
library(lance)
##GEO------
library(GEOquery)
library(Biobase)
gset<-getGEO("GSE171110",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-read.delim2('GSE171110_Data_RNAseq_raw_counts_geo.txt.gz')

dat<-expr[!duplicated(expr$Sample_sheet),]%>%column_to_rownames(var = 'Sample_sheet')

pd<-pData(gset[[1]])
colnames(dat) <- pd$geo_accession
dat <- dat[-c(1:2),]

group<-data.frame(sample=pd$geo_accession,group=pd$`status:ch1`)
table(group$group)
# 93 Coronary Artery Disease patients and 48 healthy controls
group$group <- c(rep('control',10),rep('Covid_19',44))
group<-group[order(group$group),]
dat<-dat[,group$sample]
write.table(dat,file = 'dat(GSE171110).xls',sep = '\t',row.names = T,quote = F)
write.table(group,file = 'group(GSE171110).xls',sep = '\t',row.names = F,quote = F)


## validation 1---------
library(GEOquery)
library(Biobase)
library(org.Hs.eg.db)
library(clusterProfiler)
gset<-getGEO("GSE152418",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-read.delim2('GSE152418_p20047_Study1_RawCounts.txt.gz')
gene_transform <- bitr(expr$ENSEMBLID,
                       fromType = "ENSEMBL",
                       toType = c("SYMBOL",'ENSEMBL'),
                       OrgDb = "org.Hs.eg.db")

probe2symbol<-gene_transform
colnames(probe2symbol)<-c('ENSEMBLID','symbol')
dat<-expr
dat$ENSEMBLID<-as.character(dat$ENSEMBLID)
probe2symbol$ENSEMBLID<-as.character(probe2symbol$ENSEMBLID)
dat<-dat %>%
  inner_join(probe2symbol,by='ENSEMBLID')%>% 
  dplyr::select(-ENSEMBLID)%>%     ## 去除多余信息
  dplyr::select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('S',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
pd<-pData(gset[[1]])

group<-data.frame(sample=pd$geo_accession,group=pd$`severity:ch1`)
table(group$group)
group <- group[-1,]
group$group <- ifelse(group$group=='Healthy','control','COVID_19')
group<-group[order(group$group),]
colnames(dat) <- pd$geo_accession
dat<-dat[,group$sample]
write.table(dat,file = 'dat(GSE152418).xls',sep = '\t',row.names = T,quote = F)
write.table(group,file = 'group(GSE152418).xls',sep = '\t',row.names = F,quote = F)


# RA---------
library(GEOquery)
library(Biobase)
gset<-getGEO("GSE55457",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
#library(AnnoProbe)
gpl<-getGEO("GPL96",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)    
colnames(gpl)
probe2symbol<-gpl %>%
  dplyr::select('ID','Gene Symbol')%>%
  filter('Gene Symbol'!='')%>%
  separate('Gene Symbol',c('symbol','drop'),sep = '///')%>%
  dplyr::select(-drop)
probe2symbol=probe2symbol[probe2symbol$symbol!='',]
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

group<-data.frame(sample=pd$geo_accession,group=pd$`clinical status:ch1`)
table(group$group)
group <- subset(group,group=='normal control'|group=='rheumatoid arthritis')
group$group <- ifelse(group$group=='normal control','control','RA')
group<-group[order(group$group),]
dat<-dat[,group$sample]
write.table(dat,file = 'dat(GSE55457).xls',sep = '\t',row.names = T,quote = F)
write.table(group,file = 'group(GSE55457).xls',sep = '\t',row.names = F,quote = F)


# RA validation--------
library(GEOquery)
library(Biobase)
gset<-getGEO("GSE12021",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
#library(AnnoProbe)
gpl<-getGEO("GPL96",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)    
colnames(gpl)
probe2symbol<-gpl %>%
  dplyr::select('ID','Gene Symbol')%>%
  filter('Gene Symbol'!='')%>%
  separate('Gene Symbol',c('symbol','drop'),sep = '///')%>%
  dplyr::select(-drop)
probe2symbol=probe2symbol[probe2symbol$symbol!='',]
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

group<-data.frame(sample=pd$geo_accession,group=pd$`disease:ch1`)
table(group$group)
group<-group[order(group$group),]
group <- group[-c(1:9),]
group <- group[-17,]
group$group <- c(rep('RA',12),rep('control',9))
group<-group[order(group$group),]
dat<-dat[,group$sample]
write.table(dat,file = 'dat(GSE12021).xls',sep = '\t',row.names = T,quote = F)
write.table(group,file = 'group(GSE12021).xls',sep = '\t',row.names = F,quote = F)

