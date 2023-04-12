rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/BJTC-300-8/")
if (! dir.exists("./00_rawdata")){
  dir.create("./00_rawdata")
}
setwd("./00_rawdata")
library(GEOquery)
library(Biobase)
library(tidyverse)
##1 GSE16088-----
gset<-getGEO("GSE16088",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
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
group<-data.frame(sample=pd$geo_accession,group=pd$`disease state:ch1`)
table(group$group)
group<-subset(group,group=='Normal'|group=='Tumor')
group$group<-ifelse(group$group=='Normal','control','OS')
group<-group[order(group$group),]
dat<-dat[,group$sample]
write.table(dat,file = 'dat(GSE16088).xls',sep = '\t',row.names = T,quote = F)
write.table(group,file = 'group(GSE16088).xls',sep = '\t',row.names = F,quote = F)

##2 GSE99671------
gset<-getGEO("GSE99671",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
library(readxl)
library(lance)
dat2 <- read_xlsx('GSE99671_Osteosarcoma_filecounts.xlsx')%>%column_to_rownames(var = 'symbol')%>%lc.tableToNum()
pd2 <- pData(gset[[1]])
pd2 <- pd2[c(1:36),]
group2 <- data.frame(sample=pd2$geo_accession,group=pd2$`tumor:ch1`)
colnames(dat2) <- pd2$geo_accession
table(group2$group)
group2$group <- ifelse(group2$group=='NORMAL','control','OS')
group2 <- group2[order(group2$group),]
dat2 <- dat2[,group2$sample]
write.table(dat2,file = 'dat(GSE99671).xls',sep = '\t',row.names = T,quote = F)
write.table(group2,file = 'group(GSE99671).xls',sep = '\t',row.names = F,quote = F)

##3 GSE19276------
gset<-getGEO("GSE19276",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL6848",destdir = '.')
a=gset[[1]]
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
pd<-pData(a)
library(readxl)
# group<-data.frame(sample=pd$geo_accession,group=c(rep('OS',44),rep('control',5)))
group <- read_xlsx('pd.xlsx')
table(group$group)
colnames(group) <- c('sample','group')
group<-group[order(group$group),]
dat<-dat[,group$sample]
write.table(dat,file = 'dat(GSE19276).xls',sep = '\t',row.names = T,quote = F)
write.table(group,file = 'group(GSE19276).xls',sep = '\t',row.names = F,quote = F)

##TARGET-OS---------
library(TCGAbiolinks)
library(readr)
library(readxl)
library(tidyverse)
## 读取从xena下载的数据
tcga.expr<-read_tsv(file = '/data/nas1/luchunlin/TCGA.matrix/TARGET-OS.htseq_counts.tsv.gz')
tcga.expr<-as.data.frame(tcga.expr)
tcga.expr<-column_to_rownames(tcga.expr,var = 'Ensembl_ID')
## xena下载的数据经过了log2+1转化，需要将其还原
tcga.expr<-2^tcga.expr-1
## 对数据进行id转化
genecode<-read.table(file = '/data/nas1/luchunlin/pipeline/GENEANNO/gencode.v22.annotation.gene.probeMap')
probe2symbol<-genecode[,(1:2)]
colnames(probe2symbol)<-c('ID','symbol')
probe2symbol<-probe2symbol[-1,]
dat.tcga<-tcga.expr
dat.tcga$ID <- rownames(dat.tcga)
dat.tcga$ID<-as.character(dat.tcga$ID)
probe2symbol$ID<-as.character(probe2symbol$ID)
dat.tcga<-dat.tcga %>%
  inner_join(probe2symbol,by='ID')%>% 
  dplyr::select(-ID)%>%     ## 去除多余信息
  dplyr::select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('TARGET',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
dim(dat.tcga)
## 下载临床数据
# clinical<-GDCquery_clinic(project = "TARGET-OS",type = "clinical")
pcg<-read.delim2('/data/nas1/luchunlin/pipeline/PCG/PCG.xls(v22)')
dat.final<-dat.tcga[rownames(dat.tcga)%in%pcg$gene_name,]
write.table(dat.final,file = 'dat.tcga.xls',sep = '\t',quote = F,row.names = T)
##fpkm
expr_fpkm<-read_tsv(file = '/data/nas1/luchunlin/TCGA.matrix/TARGET-OS.htseq_fpkm.tsv.gz')%>%as.data.frame()%>%
  column_to_rownames(var = 'Ensembl_ID')
## xena下载的数据经过了log2+1转化，需要将其还原
expr_fpkm<-2^expr_fpkm-1
## 对数据进行id转化
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
# dat_fpkm<-dat_fpkm[mRNA$gene_name,]
dat_fpkm<-dat_fpkm[,colnames(dat.final)]
dat_fpkm<-dat_fpkm[rownames(dat.final),]
#dat_fpkm<-na.omit(dat_fpkm)

write.table(dat_fpkm,file = 'dat.fpkm.xls',sep = '\t',row.names = T,quote = F)
