rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/JNZK-255-2/")
if (! dir.exists("./00_rawdata")){
  dir.create("./00_rawdata")
}
setwd("./00_rawdata")
##01-1 TCGA--------
library(TCGAbiolinks)
library(readr)
library(readxl)
library(tidyverse)
## 读取从xena下载的数据
tcga.expr<-read_tsv(file = '/data/nas1/luchunlin/TCGA.matrix/TCGA-BLCA.htseq_counts.tsv.gz')%>%as.data.frame()%>%
  column_to_rownames(var = 'Ensembl_ID')
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
  select(-ID)%>%     ## 去除多余信息
  select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('TCGA',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
dim(dat.tcga)

## 筛选癌症组织，去掉癌旁组织。01-09为肿瘤，10-19为正常对照
mete=data.frame(colnames(dat.tcga))  # 取第一行样本id
for (i in 1:length(mete[,1])) {
  num=as.numeric(as.character(substring(mete[i,1],14,15)))
  if(num %in% seq(1,9)){mete[i,2]="T"}
  if(num %in% seq(10,29)){mete[i,2]="N"}
}
names(mete)=c("id","group")
table(mete$group)
mete$group=as.factor(mete$group)
mete=subset(mete,mete$group=="T")
exp_tumor<-dat.tcga[,which(colnames(dat.tcga)%in%mete$id)]
exp_tumor<-as.data.frame(exp_tumor)
# 411
## 保留有生存数据的
survival<-read.delim2('/data/nas1/luchunlin/TCGA_survival/TCGA-BLCA.survival.tsv')

exp_tumor<-exp_tumor[,colnames(exp_tumor)%in%survival$sample]
##406
exp_control<-dat.tcga[,which(!colnames(dat.tcga)%in%mete$id)]
exp_control<-as.data.frame(exp_control)
# 19
dat.final<-cbind(exp_control,exp_tumor)
##425
# pcg <- read.delim2('/data/nas1/luchunlin/pipeline/PCG/PCG.xls(v22)')
# dat.final <- dat.final[pcg$gene_name,]
# dat.final <- na.omit(dat.final)
write.table(dat.final,file = 'dat.tcga.xls',sep = '\t',quote = F,row.names = T)
##fpkm
expr_fpkm<-read_tsv(file = '/data/nas1/luchunlin/TCGA.matrix/TCGA-BLCA.htseq_fpkm.tsv.gz')%>%as.data.frame()%>%
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
# dat_fpkm<-dat_fpkm[pcg$gene_name,]
# dat_fpkm<-na.omit(dat_fpkm)

write.table(dat_fpkm,file = 'dat.fpkm.xls',sep = '\t',row.names = T,quote = F)

# fpkm转TPM
FPKM2TPM <- function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

dat_tpm <- apply(dat_fpkm,2,FPKM2TPM)
write.table(dat_tpm,file = 'dat.tpm.xls',sep = '\t',row.names = T,quote = F)

###验证集--------
## GSE48075
library(GEOquery)
library(Biobase)
library(GEOquery)
library(Biobase)

gset<-getGEO("GSE1827",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL1479",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)
colnames(gpl)
library(AnnoProbe)


probe2symobl<-gpl %>%
  dplyr::select('ID','Symbol')%>%
  filter('Symbol'!='')%>%
  separate('Symbol',c('symbol','drop'),sep = '///')%>%
  dplyr::select(-drop)
probe2symobl=probe2symobl[probe2symobl$symbol!='',]
dat.va<-expr
dat.va$ID<-rownames(dat.va)
dat.va$ID<-as.character(dat.va$ID)
probe2symobl$ID<-as.character(probe2symobl$ID)
dat.va<-dat.va %>%
  inner_join(probe2symobl,by='ID')%>%
  dplyr::select(-ID)%>%     ## 去除多余信息
  dplyr::select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除

pd<-pData(a)
#table(pd$source_name_ch1)
#pd <- subset(pd,source_name_ch1=='Primary bladder cancer')

survival.va <- data.frame(sample=pd$geo_accession,OS=pd$`os censor:ch1`,OS.time=pd$`survival (mo):ch1`)

table(survival.va$OS)
survival.va$OS <- ifelse(survival.va$OS=='uncensored',0,1)
survival.va$OS<-as.numeric(survival.va$OS)
survival.va$OS.time <- as.numeric(survival.va$OS.time)*30

dat.va <- dat.va[,survival.va$sample]
write.table(dat.va,file = 'dat.va(GSE48276).xls',sep = '\t',quote = F,row.names = T)
write.table(survival.va,file = 'survival.va(GSE48276).xls',sep = '\t',row.names = F,quote = F)

#---------
library(lance)
gset<-getGEO("GSE216037",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-read.csv('GSE216037_Processed_data.csv.gz',row.names = 1)%>%lc.tableToNum()
# dat.va<-expr[,-c(1,3)]
# dat.va <- dat.va[!duplicated(dat.va$gene_name),]%>%as.data.frame()
# rownames(dat.va) <- dat.va$gene_name
# dat.va <- dat.va[,-1]


pd<-pData(gset[[1]])
#table(pd$source_name_ch1)
#pd <- subset(pd,source_name_ch1=='Primary bladder cancer')
dat.va <- expr
colnames(dat.va) <- pd$geo_accession

survival.va <- data.frame(sample=pd$geo_accession,OS=pd$`survival state:ch1`,OS.time=pd$`overall survival(month):ch1`)

table(survival.va$OS)
# survival.va$OS <- ifelse(survival.va$OS=='survival',0,1)
survival.va$OS<-as.numeric(survival.va$OS)
survival.va$OS.time <- as.numeric(survival.va$OS.time)*30

dat.va <- dat.va[,survival.va$sample]
write.table(dat.va,file = 'dat.va(GSE216037).xls',sep = '\t',quote = F,row.names = T)
write.table(survival.va,file = 'survival.va(GSE216037).xls',sep = '\t',row.names = F,quote = F)

###3-----
gset<-getGEO("GSE32894",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL6947",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)
colnames(gpl)

probe2symobl<-gpl %>%
  dplyr::select('ID','Symbol')%>%
  filter('Symbol'!='')%>%
  separate('Symbol',c('symbol','drop'),sep = '///')%>%
  dplyr::select(-drop)
probe2symobl=probe2symobl[probe2symobl$symbol!='',]
dat.va<-expr
dat.va$ID<-rownames(dat.va)
dat.va$ID<-as.character(dat.va$ID)
probe2symobl$ID<-as.character(probe2symobl$ID)
dat.va<-dat.va %>%
  inner_join(probe2symobl,by='ID')%>%
  dplyr::select(-ID)%>%     ## 去除多余信息
  dplyr::select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除

pd<-pData(a)
#table(pd$source_name_ch1)
#pd <- subset(pd,source_name_ch1=='Primary bladder cancer')

survival.va <- data.frame(sample=pd$geo_accession,OS=pd$`dod_event_(yes/no):ch1`,OS.time=pd$`time_to_dod_(months):ch1`)

table(survival.va$OS)
survival.va$OS <- ifelse(survival.va$OS=='no',0,1)
survival.va$OS<-as.numeric(survival.va$OS)
survival.va$OS.time <- as.numeric(survival.va$OS.time)*30
survival.va <- na.omit(survival.va)
dat.va <- dat.va[,survival.va$sample]
write.table(dat.va,file = 'dat.va(GSE32894).xls',sep = '\t',quote = F,row.names = T)
write.table(survival.va,file = 'survival.va(GSE32894).xls',sep = '\t',row.names = F,quote = F)
