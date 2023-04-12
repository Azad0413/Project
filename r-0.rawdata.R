rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/HF-0103-1/")
if (! dir.exists("./00_rawdata")){
  dir.create("./00_rawdata")
}
setwd("./00_rawdata")
library(GEOquery)
library(tidyverse)
library(lance)
### GSE10846-------
gset<-getGEO("GSE10846",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL570",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)    
colnames(gpl)
probe2symobl<-gpl %>%
  dplyr::select('ID','Gene Symbol')%>%
  filter('Gene Symbol'!='')%>%
  separate('Gene Symbol',c('symbol','drop'),sep = '///')%>%
  dplyr::select(-drop)
probe2symobl=probe2symobl[probe2symobl$symbol!='',]
dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symobl$ID<-as.character(probe2symobl$ID)
dat<-dat %>%
  inner_join(probe2symobl,by='ID')%>% 
  dplyr::select(-ID)%>%     ## 去除多余信息
  dplyr::select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
pcg <- read.delim2('/data/nas1/luchunlin/pipeline/PCG/PCG.xls(v22)')
dat <- dat[pcg$gene_name,]
dat <- na.omit(dat)
pd <- pData(a)
pd <- pd[c(1:414),]
table(pd$characteristics_ch1.9)
survival <- data.frame(sample=pd$geo_accession,OS=pd$characteristics_ch1.7,OS.time=pd$characteristics_ch1.8)
table(survival$OS)
survival$OS <- ifelse(survival$OS=='Clinical info: Follow up status: ALIVE',0,1)
table(survival$OS.time)
survival$OS.time <- gsub('Clinical info: Follow up years: ','',survival$OS.time,fixed = T)
survival$OS.time <- as.numeric(survival$OS.time)*365
class(survival$OS.time)
##有2个样本时间为0，需要去除
survival <- survival[which(survival$OS.time>0),]

dat<-dat[,survival$sample]

write.table(dat,file = 'dat(GSE10846).xls',sep = '\t',row.names = T,quote = F)
write.table(survival,file = 'survival(GSE10846).xls',sep = '\t',row.names = F,quote = F)

##clinical----
phenotype <- data.frame(sample=pd$geo_accession,
                        Age=pd$`Age:ch1`,
                        Gender=pd$`Gender:ch1`,
                        Stage=pd$characteristics_ch1.11,
                        ECOG_status=pd$characteristics_ch1.10,
                        LDH_Ratio=pd$characteristics_ch1.12)
write.table(phenotype,file = 'phenotype.xls',sep = '\t',row.names = F,quote = F)


### GSE53786-------
gset<-getGEO("GSE53786",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL570",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)
colnames(gpl)
probe2symobl<-gpl %>%
  dplyr::select('ID','Gene Symbol')%>%
  filter('Gene Symbol'!='')%>%
  separate('Gene Symbol',c('symbol','drop'),sep = '///')%>%
  dplyr::select(-drop)
probe2symobl=probe2symobl[probe2symobl$symbol!='',]
dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symobl$ID<-as.character(probe2symobl$ID)
dat<-dat %>%
  inner_join(probe2symobl,by='ID')%>%
  dplyr::select(-ID)%>%     ## 去除多余信息
  dplyr::select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除

pd <- pData(a)

survival <- data.frame(sample=pd$geo_accession,OS=pd$characteristics_ch1.7,OS.time=pd$characteristics_ch1.8)
table(survival$OS)
survival$OS <- ifelse(survival$OS=='clinical info: Follow up status: Alive',0,ifelse(survival$OS=='Clinical info: Follow up status: ALIVE',0,1))
table(survival$OS.time)
survival$OS.time <- gsub('Clinical info: Follow up years: ','',survival$OS.time,fixed = T)
survival$OS.time <- gsub('clinical info: Follow up years: ','',survival$OS.time,fixed = T)
survival$OS.time <- as.numeric(survival$OS.time)*365
class(survival$OS.time)
##有2个样本时间为0，需要去除
# survival <- survival[which(survival$OS.time>0),]

dat<-dat[,survival$sample]

write.table(dat,file = 'dat(GSE53786).xls',sep = '\t',row.names = T,quote = F)
write.table(survival,file = 'survival(GSE53786).xls',sep = '\t',row.names = F,quote = F)


# 
# ##01-1 TCGA--------
# library(TCGAbiolinks)
# library(readr)
# library(readxl)
# library(tidyverse)
# ## 读取从xena下载的数据
# tcga.expr<-read_tsv(file = '/data/nas1/luchunlin/TCGA.matrix/TCGA-DLBC.htseq_counts.tsv.gz')%>%as.data.frame()%>%
#   column_to_rownames(var = 'Ensembl_ID')
# ## xena下载的数据经过了log2+1转化，需要将其还原
# tcga.expr<-2^tcga.expr-1
# ## 对数据进行id转化
# genecode<-read.table(file = '/data/nas1/luchunlin/pipeline/GENEANNO/gencode.v22.annotation.gene.probeMap')
# probe2symbol<-genecode[,(1:2)]
# colnames(probe2symbol)<-c('ID','symbol')
# probe2symbol<-probe2symbol[-1,]
# dat.tcga<-tcga.expr
# dat.tcga$ID <- rownames(dat.tcga)
# dat.tcga$ID<-as.character(dat.tcga$ID)
# probe2symbol$ID<-as.character(probe2symbol$ID)
# dat.tcga<-dat.tcga %>%
#   inner_join(probe2symbol,by='ID')%>% 
#   select(-ID)%>%     ## 去除多余信息
#   select(symbol,everything())%>%     ## 重新排列
#   mutate(rowMean=rowMeans(.[grep('TCGA',names(.))]))%>%    ## 求出平均数
#   arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
#   distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
#   select(-rowMean)%>%     ## 反向选择去除rowMean这一列
#   tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
# dim(dat.tcga)
# exp_tumor <- dat.tcga
# ## 保留有生存数据的
# survival<-read.delim2('/data/nas1/luchunlin/TCGA_survival/TCGA-DLBC.survival.tsv')
# exp_tumor<-exp_tumor[,colnames(exp_tumor)%in%survival$sample]
# ##47
# # pcg <- read.delim2('/data/nas1/luchunlin/pipeline/PCG/PCG.xls(v22)')
# # dat.final<-exp_tumor[pcg$gene_name,]
# # dat.final <- na.omit(dat.final)
# dat.final <- exp_tumor
# ##524
# write.table(dat.final,file = 'dat.tcga.xls',sep = '\t',quote = F,row.names = T)
# ##fpkm
# expr_fpkm<-read_tsv(file = '/data/nas1/luchunlin/TCGA.matrix/TCGA-DLBC.htseq_fpkm.tsv.gz')%>%as.data.frame()%>%
#   column_to_rownames(var = 'Ensembl_ID')
# ## xena下载的数据经过了log2+1转化，需要将其还原
# expr_fpkm<-2^expr_fpkm-1
# ## 对数据进行id转化
# dat_fpkm<-expr_fpkm
# dat_fpkm$ID <- rownames(dat_fpkm)
# dat_fpkm$ID<-as.character(dat_fpkm$ID)
# probe2symbol$ID<-as.character(probe2symbol$ID)
# dat_fpkm<-dat_fpkm %>%
#   inner_join(probe2symbol,by='ID')%>% 
#   select(-ID)%>%     ## 去除多余信息
#   select(symbol,everything())%>%     ## 重新排列
#   mutate(rowMean=rowMeans(.[grep('TCGA',names(.))]))%>%    ## 求出平均数
#   arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
#   distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
#   select(-rowMean)%>%     ## 反向选择去除rowMean这一列
#   tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
# dim(dat_fpkm)
# 
# dat_fpkm<-dat_fpkm[,colnames(dat.final)]
# dat_fpkm<-dat_fpkm[rownames(dat.final),]
# dat_fpkm <- na.omit(dat_fpkm)
# write.table(dat_fpkm,file = 'dat.fpkm.xls',sep = '\t',row.names = T,quote = F)
# # dat.final <- dat.final[rownames(dat_fpkm),]
# # write.table(dat.final,file = 'dat.tcga.xls',sep = '\t',quote = F,row.names = T)
# # # fpkm转TPM
# # FPKM2TPM <- function(fpkm){
# #   exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
# # }
# # 
# # dat_tpm <- apply(dat_fpkm,2,FPKM2TPM)
# # write.table(dat_tpm,file = 'dat.tpm.xls',sep = '\t',row.names = T,quote = F)