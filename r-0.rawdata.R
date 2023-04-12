rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/BJTC-149(modify)/")
if (! dir.exists("./00_rawdata")){
  dir.create("./00_rawdata")
}
setwd("./00_rawdata")

##GSE49541---------
library(GEOquery)
library(Biobase)
library(tidyverse)
gset<-getGEO("GSE49541",
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

# group<-data.frame(sample=pd$geo_accession,group=pd$`group:ch1`)
# table(group$group)
# #group<-subset(group,group=='Normal'|group=='Tumor')
# group$group<-ifelse(group$group=='normal weight','Normal','Obese')
# group<-group[order(group$group),]
# dat<-dat[,group$sample]
write.table(dat,file = 'dat(GSE49541).xls',sep = '\t',row.names = T,quote = F)
# write.table(group,file = 'group(GSE49541).xls',sep = '\t',row.names = F,quote = F)


##GSE151158---------
library(GEOquery)
library(Biobase)
library(tidyverse)
gset<-getGEO("GSE159088",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL21185",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)    
colnames(gpl)
probe2symbol<-gpl %>%
  dplyr::select('ID','GENE_SYMBOL')%>%
  filter('GENE_SYMBOL'!='')%>%
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

# group<-data.frame(sample=pd$geo_accession,group=pd$`group:ch1`)
# table(group$group)
# #group<-subset(group,group=='Normal'|group=='Tumor')
# group$group<-ifelse(group$group=='normal weight','Normal','Obese')
# group<-group[order(group$group),]
# dat<-dat[,group$sample]
write.table(dat,file = 'dat(GSE159088).xls',sep = '\t',row.names = T,quote = F)
# write.table(group,file = 'group(GSE151158).xls',sep = '\t',row.names = F,quote = F)


##GSE130970--------
library(GEOquery)
library(Biobase)
library(tidyverse)
gset<-getGEO("GSE130970",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-read_csv('GSE130970_all_sample_salmon_tximport_TPM_entrez_gene_ID.csv.gz')
id2symbol=AnnotationDbi::select(org.Hs.eg.db, as.character(expr$entrez_id), "SYMBOL", "ENTREZID")
expr$entrez_id <- id2symbol$SYMBOL
expr <- na.omit(expr)
expr <- expr[!duplicated(expr$entrez_id),]%>%as.data.frame()%>%column_to_rownames(var = 'entrez_id')

pd<-pData(gset[[1]])
dat <- expr
# group<-data.frame(sample=pd$geo_accession,group=pd$`group:ch1`)
# table(group$group)
# #group<-subset(group,group=='Normal'|group=='Tumor')
# group$group<-ifelse(group$group=='normal weight','Normal','Obese')
# group<-group[order(group$group),]
# dat<-dat[,group$sample]
write.table(dat,file = 'dat(GSE130970).xls',sep = '\t',row.names = T,quote = F)
# write.table(group,file = 'group(GSE151158).xls',sep = '\t',row.names = F,quote = F)


