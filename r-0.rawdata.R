rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/YCZK-127")
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
tcga.expr<-read_tsv(file = 'TCGA-BRCA.htseq_counts.tsv')
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
#keep<-rowSums(dat.tcga>0)>=floor(0.75*ncol(dat.tcga))
#dat.final<-dat.tcga[keep,]


## 下载临床数据
clinical<-GDCquery_clinic(project = "TCGA-BRCA",type = "clinical")
## 筛选癌症组织，去掉癌旁组织。01-09为肿瘤，10-19为正常对照
mete=data.frame(colnames(dat.tcga))  # 取第一行样本id
for (i in 1:length(mete[,1])) {
  num=as.numeric(as.character(substring(mete[i,1],14,15)))
  if(num %in% seq(1,9)){mete[i,2]="T"}
  if(num %in% seq(10,29)){mete[i,2]="N"}
}
names(mete)=c("id","group")
mete$group=as.factor(mete$group)
mete=subset(mete,mete$group=="T")
exp_tumor<-dat.tcga[,which(colnames(dat.tcga)%in%mete$id)]
exp_tumor<-as.data.frame(exp_tumor)

# 1104
exp_control<-dat.tcga[,which(!colnames(dat.tcga)%in%mete$id)]
exp_control<-as.data.frame(exp_control)
# 113
dat.final<-cbind(exp_control,exp_tumor)

write.table(dat.final,file = 'dat.tcga.xls',sep = '\t',quote = F,row.names = T)

##fpkm
expr_fpkm<-read_tsv(file = 'TCGA-BRCA.htseq_fpkm.tsv')
expr_fpkm<-as.data.frame(expr_fpkm)
rownames(expr_fpkm)<-expr_fpkm[,1]
expr_fpkm<-expr_fpkm[,-1]
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
write.table(dat_fpkm,file = 'dat.fpkm.xls',sep = '\t',row.names = T,quote = F)
# fpkm转TPM
FPKM2TPM <- function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

dat_tpm <- apply(dat_fpkm,2,FPKM2TPM)
write.table(dat_tpm,file = 'dat.tpm.xls',sep = '\t',row.names = T,quote = F)
## 01-2 GSE42568------
library(GEOquery)
library(Biobase)
gset<-getGEO("GSE42568",
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
MIR7703<-expr['201762_s_at',]
rownames(MIR7703)<-'PSME2'
##201762_s_at
dat<-rbind(MIR7703,dat)
write.table(dat,file = 'dat.va.xls',sep = '\t',quote = F,row.names = T)
group<-data.frame(sample=pd$geo_accession,
                  group=pd$`tissue:ch1`)
survival<-data.frame(sample=pd$geo_accession,OS=pd$`overall survival event:ch1`,OS.time=pd$`overall survival time_days:ch1`)
write.table(survival,file = 'survival.va.xls',sep = '\t',row.names = F,quote = F)
#group$group<-c(rep('control',17),rep('Tumor',104))
#write.table(group,file = 'group.xls',sep = '\t',row.names = F,quote = F)
