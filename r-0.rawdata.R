rm(list = ls())
# 01 获取数据集--------------
setwd("/data/nas1/luchunlin/project/JNZK-218-8/")
if (! dir.exists("./00_rawdata")){
  dir.create("./00_rawdata")
}
setwd("./00_rawdata")
library(GEOquery)
library(Biobase)
library(limma)
library(tidyverse)
## GSE27276
gset<-getGEO("GSE27276",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
a=gset[[1]]
pd<-pData(a)
library(AnnoProbe)
gpl<-getGEO("GPL2507",destdir = '.')
gpl1<-Table(gpl)                  
colnames(Table(gpl)) 
## 没有Symbol名，先换成GB_ACC     
probe2GB<-gpl1[,c(1,3)]           
colnames(probe2GB)=c('probe_id','GB_ACC')
length(unique(probe2GB$GB_ACC))
# [1] 25360
ids=probe2GB[probe2GB$GB_ACC!='',]
#判断是否匹配
ids=probe2GB[probe2GB$probe_id%in%rownames(expr),]
dat=expr[ids$probe_id,]
ids$mean=apply(dat,1,mean)
ids=ids[order(ids$GB_ACC,ids$mean,decreasing = T),]
ids<-separate(ids,GB_ACC,into = c('GB_ACC'),sep = '\\.')
ids=ids[!duplicated(ids$GB_ACC),]
dat=dat[ids$probe_id,]
rownames(dat)=ids$GB_ACC
dat=dat[-25354,]
library(org.Hs.eg.db)
library(clusterProfiler)
geneid<-rownames(dat)
gene_transform<-bitr(geneid,fromType = "REFSEQ",
                     toType = "SYMBOL",
                     OrgDb = org.Hs.eg.db)
dat2<-dat
colnames(gene_transform)=c('REFSEQ','Symbol')
gene_transform=gene_transform[!duplicated(gene_transform$Symbol),]
dat2=dat2[gene_transform$REFSEQ,]
rownames(dat2)=gene_transform$Symbol
## 16351
write.table(dat2,file = "dat(GSE27276).xls",
            quote = F,
            sep = '\t',
            row.names = T)
group<-data.frame(sample=pd$geo_accession,group=pd$`disease state:ch1`)
table(group$group)
group$group<-ifelse(group$group=='Control','control','POAG')
write.table(group,file = 'group(GSE27276).xls',sep = '\t',row.names = F,quote = F)


## GSE100927
gset<-getGEO("GSE100927",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
a=gset[[1]]
pd<-pData(a)
gpl<-getGEO("GPL17077",destdir = '.')
gpl<-Table(gpl)    
colnames(gpl)
probe2symbol<-gpl %>%
  dplyr::select('ID','GENE_SYMBOL')%>%
  filter('GENE_SYMBOL'!='')%>%
  separate('GENE_SYMBOL',c('symbol','drop'),sep = '///')%>%
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

table(pd$source_name_ch1)
group<-data.frame(sample=pd$geo_accession,group=pd$source_name_ch1)
group<-subset(group,group=='Atherosclerotic carotid artery'|group=='Control carotid artery')
table(group$group)
group$group<-ifelse(group$group=='Atherosclerotic carotid artery','AS','control')
group<-group[order(group$group),]
dat<-dat[,group$sample]
write.table(dat,file = 'dat(GSE100927).xls',sep = '\t',row.names = T,quote = F)
write.table(group,file = 'group(GSE100927).xls',sep = '\t',row.names = F,quote = F)
