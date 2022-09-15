##具体分析点的图包括：表达差异分析（柱状图）；
## 临床分析（箱线图）；
## 一致性聚类（累积分布曲线，矩阵热图等）；
## 差异基因（火山图，热图）；GO和KEGG（气泡图，和弦图）；
## TOP5强相关基因的相关性分析（相关性热图，加上两个目的基因）。

rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/LLZK-505")
if (! dir.exists("./00_rawdata")){
  dir.create("./00_rawdata")
}
setwd("./00_rawdata")
##01-1 TCGA OV--------
library(TCGAbiolinks)
library(readr)
library(readxl)
library(tidyverse)
## 读取从xena下载的数据
tcga.expr<-read_tsv(file = 'TCGA-OV.htseq_counts.tsv')
tcga.expr<-as.data.frame(tcga.expr)
rownames(tcga.expr)<-tcga.expr[,1]
tcga.expr<-tcga.expr[,-1]
## xena下载的数据经过了log2+1转化，需要将其还原
tcga.expr<-2^tcga.expr-1
## 对数据进行id转化
genecode<-read.table(file = 'gencode.v22.annotation.gene.probeMap')
probe2symbol<-genecode[,(1:2)]
colnames(probe2symbol)<-c('ID','symbol')
probe2symbol<-probe2symbol[-1,]
dat.tcga<-tcga.expr
dat.tcga$ID <- rownames(dat.tcga)
dat.tcga$ID<-as.character(dat.tcga$ID)
probe2symbol$ID<-as.character(probe2symbol$ID)
dat.tcga<-dat.tcga %>%
  inner_join(probe2symbol,by='ID')%>% 
  select(-ID)%>%     ## 去除多余信息
  select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('TCGA',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
dim(dat.tcga)
keep<-rowSums(dat.tcga>0)>=floor(0.75*ncol(dat.tcga))
dat.final<-dat.tcga[keep,]
write.table(dat.final,file = 'dat.tcga.xls',sep = '\t',quote = F,row.names = T)

##fpkm
expr_fpkm<-read_tsv(file = 'TCGA-OV.htseq_fpkm.tsv')
expr_fpkm<-as.data.frame(expr_fpkm)
rownames(expr_fpkm)<-expr_fpkm[,1]
expr_fpkm<-expr_fpkm[,-1]
## xena下载的数据经过了log2+1转化，需要将其还原
expr_fpkm<-2^expr_fpkm-1
## 对数据进行id转化
genecode<-read.table(file = 'gencode.v22.annotation.gene.probeMap')
probe2symbol<-genecode[,(1:2)]
colnames(probe2symbol)<-c('ID','symbol')
probe2symbol<-probe2symbol[-1,]
dat_fpkm<-expr_fpkm
dat_fpkm$ID <- rownames(dat_fpkm)
dat_fpkm$ID<-as.character(dat_fpkm$ID)
probe2symbol$ID<-as.character(probe2symbol$ID)
dat_fpkm<-dat_fpkm %>%
  inner_join(probe2symbol,by='ID')%>% 
  select(-ID)%>%     ## 去除多余信息
  select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('TCGA',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
dim(dat_fpkm)
dat_fpkm<-dat_fpkm[rownames(dat_fpkm)%in%rownames(dat.final),]
write.table(dat_fpkm,file = 'dat.fpkm.xls',sep = '\t',row.names = T,quote = F)
# fpkm转TPM
FPKM2TPM <- function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

dat_tpm <- apply(dat_fpkm,2,FPKM2TPM)
write.table(dat_tpm,file = 'dat.tpm.xls',sep = '\t',row.names = T,quote = F)

## 01-2 GSE26712------
library(GEOquery)
library(Biobase)
gset<-getGEO("GSE26712",
              destdir = '.',
              GSEMatrix = T,
              getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL96",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)    
colnames(gpl)
probe2symobl<-gpl %>%
  select('ID','Gene Symbol')%>%
  filter('Gene Symbol'!='')%>%
  separate('Gene Symbol',c('symbol','drop'),sep = '///')%>%
  select(-drop)
probe2symobl=probe2symobl[probe2symobl$symbol!='',]
dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symobl$ID<-as.character(probe2symob$ID)
dat<-dat %>%
  inner_join(probe2symobl,by='ID')%>% 
  select(-ID)%>%     ## 去除多余信息
  select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
pd<-pData(a)

write.table(dat,file = 'dat.xls',sep = '\t',quote = F,row.names = T)
group<-data.frame(sample=pd$geo_accession,
                  group=pd$title)
group$group<-c(rep('control',10),rep('OV',185))
write.table(group,file = 'group.xls',sep = '\t',row.names = F,quote = F)
## 01-2 GSE14407------
library(GEOquery)
library(Biobase)
gset<-getGEO("GSE14407",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL570",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)    
colnames(gpl)
probe2symobl<-gpl %>%
  select('ID','Gene Symbol')%>%
  filter('Gene Symbol'!='')%>%
  separate('Gene Symbol',c('symbol','drop'),sep = '///')%>%
  select(-drop)
probe2symobl=probe2symobl[probe2symobl$symbol!='',]
dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symobl$ID<-as.character(probe2symobl$ID)
dat<-dat %>%
  inner_join(probe2symobl,by='ID')%>% 
  select(-ID)%>%     ## 去除多余信息
  select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
pd<-pData(a)
write.table(dat,file = 'dat(GSE14407).xls',sep = '\t',quote = F,row.names = T)
group<-data.frame(sample=pd$geo_accession,
                  group=pd$`disease state:ch1`)

group$group<-c(rep('control',12),rep('OV',12))
write.table(group,file = 'group(GSE14407).xls',sep = '\t',row.names = F,quote = F)

## 01-2 GSE105437------
library(GEOquery)
library(Biobase)
gset<-getGEO("GSE105437",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL570",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)    
colnames(gpl)
probe2symobl<-gpl %>%
  select('ID','Gene Symbol')%>%
  filter('Gene Symbol'!='')%>%
  separate('Gene Symbol',c('symbol','drop'),sep = '///')%>%
  select(-drop)
probe2symobl=probe2symobl[probe2symobl$symbol!='',]
dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symobl$ID<-as.character(probe2symobl$ID)
dat<-dat %>%
  inner_join(probe2symobl,by='ID')%>% 
  select(-ID)%>%     ## 去除多余信息
  select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
pd<-pData(a)
write.table(dat,file = 'dat(GSE105437).xls',sep = '\t',quote = F,row.names = T)
group<-data.frame(sample=pd$geo_accession,
                  group=pd$title)
group<-group[-c(6:12),]
group$group<-c(rep('control',5),rep('OV',10))
write.table(group,file = 'group(GSE105437).xls',sep = '\t',row.names = F,quote = F)

## 01-2 GSE14001------
library(GEOquery)
library(Biobase)
gset<-getGEO("GSE14001",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))

a=gset[[1]]
dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symobl$ID<-as.character(probe2symobl$ID)
dat<-dat %>%
  inner_join(probe2symobl,by='ID')%>% 
  select(-ID)%>%     ## 去除多余信息
  select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
pd<-pData(a)
write.table(dat,file = 'dat(GSE14001).xls',sep = '\t',quote = F,row.names = T)
group<-data.frame(sample=pd$geo_accession,
                  group=pd$title)
group$group<-c(rep('control',1),rep('OV',20),rep('control','2'))
write.table(group,file = 'group(GSE14001).xls',sep = '\t',row.names = F,quote = F)
