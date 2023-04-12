rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/XA-0214-1/")
if (! dir.exists("./00_rawdata")){
  dir.create("./00_rawdata")
}
setwd("./00_rawdata")
##GSE97537---------
library(GEOquery)
library(Biobase)
library(tidyverse)
gset<-getGEO("GSE97537",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)

expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL1355",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)    
colnames(gpl)
probe2symbol<-gpl %>%
  dplyr::select('ID','Gene Symbol')%>%
  filter('Gene Symbol'!='')%>%
  separate('Gene Symbol',c('symbol','drop'),sep = ' ')%>%
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

group<-data.frame(sample=pd$geo_accession,group=pd$`stress:ch1`)
table(group$group)
#group<-subset(group,group=='Normal'|group=='Tumor')
group$group<-ifelse(group$group=='Sham','control','CIRI')
group<-group[order(group$group),]
dat<-dat[,group$sample]
write.table(dat,file = 'dat(GSE97537).xls',sep = '\t',row.names = T,quote = F)
write.table(group,file = 'group(GSE97537).xls',sep = '\t',row.names = F,quote = F)


##GSE61616---------
library(GEOquery)
library(Biobase)
library(tidyverse)
gset<-getGEO("GSE61616",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL1355",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)    
colnames(gpl)
probe2symbol<-gpl %>%
  dplyr::select('ID','Gene Symbol')%>%
  filter('Gene Symbol'!='')%>%
  separate('Gene Symbol',c('symbol','drop'),sep = ' ')%>%
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

group<-data.frame(sample=pd$geo_accession,group=pd$title)
table(group$group)
group <- group[c(1:10),]
group$group <- c(rep('control',5),rep('CIRI',5))

#group<-subset(group,group=='Normal'|group=='Tumor')
# group$group<-ifelse(group$group=='Sham','control','CIRI')
group<-group[order(group$group),]
dat<-dat[,group$sample]
write.table(dat,file = 'dat(GSE61616).xls',sep = '\t',row.names = T,quote = F)
write.table(group,file = 'group(GSE61616).xls',sep = '\t',row.names = F,quote = F)


##GSE78731---------
library(GEOquery)
library(Biobase)
library(tidyverse)
gset<-getGEO("GSE78731",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)

expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL15084",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)    
colnames(gpl)
probe2symbol<-gpl %>%
  dplyr::select('ID','GENE_SYMBOL')%>%
  filter('GENE_SYMBOLl'!='')%>%
  separate('GENE_SYMBOL',c('symbol','drop'),sep = ' ')%>%
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


group<-data.frame(sample=pd$geo_accession,group=pd$title)
table(group$group)
group <- group[-c(5:10),]
group$group <- c(rep('control',4),rep('CIRI',6))

#group<-subset(group,group=='Normal'|group=='Tumor')
# group$group<-ifelse(group$group=='Sham','control','CIRI')
group<-group[order(group$group),]
dat<-dat[,group$sample]
write.table(dat,file = 'dat(GSE78731).xls',sep = '\t',row.names = T,quote = F)
write.table(group,file = 'group(GSE78731).xls',sep = '\t',row.names = F,quote = F)
