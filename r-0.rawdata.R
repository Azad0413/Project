rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/SJZZK-428-10/")
if (! dir.exists("./00_rawdata")){
  dir.create("./00_rawdata")
}
setwd("./00_rawdata")

## TCGA-OV 用BJTC-258整理的数据

##GEO
library(GEOquery)
library(Biobase)
library(tidyverse)
##1 GSE66957-----
gset<-getGEO("GSE66957",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL15048",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)    
colnames(gpl)
probe2symbol<-gpl %>%
  dplyr::select('ID','GeneSymbol')%>%
  filter('GeneSymbol'!='')%>%
  separate('GeneSymbol',c('symbol','drop'),sep = ' ')%>%
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
group<-data.frame(sample=pd$geo_accession,group=pd$`tissue:ch1`)
table(group$group)
#group<-subset(group,group=='Normal'|group=='Tumor')
group$group<-ifelse(group$group=='normal ovary','control','OV')
group<-group[order(group$group),]
dat<-dat[,group$sample]
write.table(dat,file = 'dat(GSE66957).xls',sep = '\t',row.names = T,quote = F)
write.table(group,file = 'group(GSE66957).xls',sep = '\t',row.names = F,quote = F)

##2 GSE18520 GSE40595 
gset<-getGEO("GSE18520",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL570",destdir = '.')
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
group<-data.frame(sample=pd$geo_accession,group=pd$title)
group$group <- c(rep('OV',53),rep('control',10))
table(group$group)
#group<-subset(group,group=='Normal'|group=='Tumor')
#group$group<-ifelse(group$group=='normal ovary','control','OV')
group<-group[order(group$group),]
dat<-dat[,group$sample]
write.table(dat,file = 'dat(GSE18520).xls',sep = '\t',row.names = T,quote = F)
write.table(group,file = 'group(GSE18520).xls',sep = '\t',row.names = F,quote = F)


##2 GSE40595 
gset<-getGEO("GSE40595",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL570",destdir = '.')
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
group<-data.frame(sample=pd$geo_accession,group=pd$`tissue:ch1`)
table(group$group)
group <- subset(group,group=='Microdissected normal ovarian stroma'|group=='Microdissected ovarian cancer stroma')
group$group<-ifelse(group$group=='Microdissected normal ovarian stroma','control','OV')
group<-group[order(group$group),]
dat<-dat[,group$sample]
write.table(dat,file = 'dat(GSE40595).xls',sep = '\t',row.names = T,quote = F)
write.table(group,file = 'group(GSE40595).xls',sep = '\t',row.names = F,quote = F)
