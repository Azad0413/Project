rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/NN-0118-2/")
if (! dir.exists("./00_rawdata")){
  dir.create("./00_rawdata")
}
setwd("./00_rawdata")
##01-1 TCGA -------
library(TCGAbiolinks)
library(readr)
library(readxl)
library(tidyverse)
library(lance)
## 读取从xena下载的数据
tcga.expr1<-read_tsv(file = '/data/nas1/luchunlin/TCGA.matrix/TCGA-LUAD.htseq_counts.tsv.gz')
tcga.expr2<-read_tsv(file = '/data/nas1/luchunlin/TCGA.matrix/TCGA-LUSC.htseq_counts.tsv.gz')
tcga.expr <- merge(tcga.expr1,tcga.expr2,by='Ensembl_ID')%>%as.data.frame()%>%
  column_to_rownames(var = 'Ensembl_ID')%>%lc.tableToNum()
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
  mutate(rowMean=rowMeans(.[grep('TCGA',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
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
# 1027
## 保留有生存数据的
survival1<-read.delim2('/data/nas1/luchunlin/TCGA_survival/TCGA-LUAD.survival.tsv')
survival2 <- read.delim2('/data/nas1/luchunlin/TCGA_survival/TCGA-LUSC.survival.tsv')
survival <- rbind(survival1,survival2)
write.table(survival,file = 'survival.xls',sep = '\t',row.names = F,quote = F)

exp_tumor<-exp_tumor[,colnames(exp_tumor)%in%survival$sample]
# ##1006
exp_control<-dat.tcga[,which(!colnames(dat.tcga)%in%mete$id)]
exp_control<-as.data.frame(exp_control)
# 108
dat.final<-cbind(exp_control,exp_tumor)

##1114
pcg <- read.delim2('/data/nas1/luchunlin/pipeline/PCG/PCG.xls(v22)')
dat.final <- dat.final[pcg$gene_name,]
dat.final <- na.omit(dat.final)
write.table(dat.final,file = 'dat.tcga.xls',sep = '\t',quote = F,row.names = T)
##fpkm
expr_fpkm1<-read_tsv(file = '/data/nas1/luchunlin/TCGA.matrix/TCGA-LUAD.htseq_fpkm.tsv.gz')
expr_fpkm2<-read_tsv(file = '/data/nas1/luchunlin/TCGA.matrix/TCGA-LUSC.htseq_fpkm.tsv.gz')
expr_fpkm <- merge(expr_fpkm1,expr_fpkm2,by='Ensembl_ID')%>%as.data.frame()%>%
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
  dplyr::select(-ID)%>%     ## 去除多余信息
  dplyr::select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('TCGA',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
dim(dat_fpkm)
# dat_fpkm<-dat_fpkm[mRNA$gene_name,]
dat_fpkm<-dat_fpkm[,colnames(dat.final)]
dat_fpkm<-dat_fpkm[pcg$gene_name,]
dat_fpkm<-na.omit(dat_fpkm)
write.table(dat_fpkm,file = 'dat.fpkm.xls',sep = '\t',row.names = T,quote = F)

# fpkm转TPM
FPKM2TPM <- function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

dat_tpm <- apply(dat_fpkm,2,FPKM2TPM)
write.table(dat_tpm,file = 'dat.tpm.xls',sep = '\t',row.names = T,quote = F)

library(GEOquery)
###GSE40275--------
gset<-getGEO("GSE40275",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL15974",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)    
colnames(gpl)
probe2symobl<-gpl %>%
  dplyr::select('ID','ENTREZ_GENE_ID')%>%
  filter('ENTREZ_GENE_ID'!='')

library(clusterProfiler)
library(org.Hs.eg.db)
gene_transform <- bitr(probe2symobl$ENTREZ_GENE_ID,
                       fromType = "ENTREZID",
                       toType = c("SYMBOL"),
                       OrgDb = "org.Hs.eg.db")
colnames(probe2symobl)[2] <- 'ENTREZID'
probe2symobl <- merge(probe2symobl,gene_transform,by='ENTREZID')
colnames(probe2symobl)[3] <- 'symbol'

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

table(pd$source_name_ch1)
group<-data.frame(sample=pd$geo_accession,title=pd$title,group=pd$source_name_ch1)
group <- subset(group,group=='normal lung'|group=='Adenocarcinoma of lung'|group=='Carcinoma of lung, large cell'|
               group=='Carcinoma of lung, neuroendocrine'|group=='Carcinoma of lung, non-small cell'|
               group=='Carcinoma of lung, squamous cell')

table(group$group)
group$group<-ifelse(group$group=='normal lung','Normal','Tumor')
group<-group[order(group$group),]%>%dplyr::select(-'title')
dat<-dat[,group$sample]

write.table(dat,file = 'dat(GSE40275).xls',sep = '\t',row.names = T,quote = F)
write.table(group,file = 'group(GSE40275).xls',sep = '\t',row.names = F,quote = F)


### GSE21933-------
gset<-getGEO("GSE21933",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL6254",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)    
colnames(gpl)
probe2symobl<-gpl %>%
  dplyr::select('ID','GENE_SYMBOL')%>%
  filter('GENE_SYMBOL'!='')%>%
  separate('GENE_SYMBOL',c('symbol','drop'),sep = '///')%>%
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
###TOP2A
top2a <- expr['PH_hs_0009437',]
rownames(top2a) <- 'TOP2A'
dat <- rbind(top2a,dat)
pd <- pData(a)
table(pd$`tissue:ch1`)
group<-data.frame(sample=pd$geo_accession,group=pd$`tissue:ch1`)
table(group$group)
group$group<-ifelse(group$group=='primary normal lung tissue','Normal','Tumor')
group<-group[order(group$group),]
dat<-dat[,group$sample]

write.table(dat,file = 'dat(GSE21933).xls',sep = '\t',row.names = T,quote = F)
write.table(group,file = 'group(GSE21933).xls',sep = '\t',row.names = F,quote = F)

###

### GSE18842-------
gset<-getGEO("GSE18842",
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
table(pd$`sample type:ch1`)
group<-data.frame(sample=pd$geo_accession,group=pd$`sample type:ch1`)
table(group$group)

group$group<-ifelse(group$group=='control','Normal','Tumor')
group<-group[order(group$group),]
dat<-dat[,group$sample]

write.table(dat,file = 'dat(GSE18842).xls',sep = '\t',row.names = T,quote = F)
write.table(group,file = 'group(GSE18842).xls',sep = '\t',row.names = F,quote = F)

### GSE19188-------
gset<-getGEO("GSE19188",
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
table(pd$`tissue type:ch1`)
group<-data.frame(sample=pd$geo_accession,group=pd$`tissue type:ch1`)
table(group$group)

group$group<-ifelse(group$group=='healthy','Normal','Tumor')
group<-group[order(group$group),]
dat<-dat[,group$sample]

write.table(dat,file = 'dat(GSE19188).xls',sep = '\t',row.names = T,quote = F)
write.table(group,file = 'group(GSE19188).xls',sep = '\t',row.names = F,quote = F)


