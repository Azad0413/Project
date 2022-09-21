rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/BJTC-357")
if (! dir.exists("./00_rawdata")){
  dir.create("./00_rawdata")
}
setwd("./00_rawdata")
##01-1 TCGA ESCA--------
library(TCGAbiolinks)
library(readr)
library(readxl)
library(tidyverse)
## 读取从xena下载的数据
tcga.expr<-read_tsv(file = 'TCGA-ESCA.htseq_counts.tsv.gz')
tcga.expr<-as.data.frame(tcga.expr)
tcga.expr<-column_to_rownames(tcga.expr,var = "Ensembl_ID")
## xena下载的数据经过了log2+1转化，需要将其还原
tcga.expr<-2^tcga.expr-1
## 加载注释文件
library("rtracklayer")
  #gtf_data = import('/data/nas1/luchunlin/pipeline/GENEANNO/gencode.v22.annotation.gtf.gz') #gtf的路径
  gtf_data = import('gencode.v41.long_noncoding_RNAs.gtf.gz') #gtf的路径
  gtf_data = as.data.frame(gtf_data)
table(gtf_data$gene_type)
lncRNA=gtf_data%>%
  dplyr::filter(type=="gene",gene_type=="lncRNA")%>%
  dplyr::select(gene_id,gene_type,gene_name)
write.table(lncRNA,file = 'lncRNA.xls',sep = '\t',row.names = F,quote = F)
## 对数据进行id转化
genecode<-read.table(file = 'gencode.v22.annotation.gene.probeMap')
probe2symbol<-genecode[,(1:2)]
colnames(probe2symbol)<-c('ID','symbol')
probe2symbol<-probe2symbol[-1,]
dat.tcga<-tcga.expr
#dat.tcga<-tcga.expr[rownames(tcga.expr)%in%lncRNA$gene_id,]
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
#at.final<-dat.tcga[keep,]
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
## 162
exp_control<-dat.tcga[,which(!colnames(dat.tcga)%in%mete$id)]
exp_control<-as.data.frame(exp_control)
# 11

phenotype<-read_tsv(file = 'TCGA-ESCA.GDC_phenotype.tsv.gz')
table(phenotype$primary_diagnosis.diagnoses)
clinical_ESCC<-subset(phenotype,primary_diagnosis.diagnoses=='Squamous cell carcinoma, keratinizing, NOS'|primary_diagnosis.diagnoses=='Squamous cell carcinoma, NOS')
clinical_ESCC<-data.frame(sample=clinical_ESCC$submitter_id.samples,type=clinical_ESCC$primary_diagnosis.diagnoses)
exp_tumor<-exp_tumor[,colnames(exp_tumor)%in%clinical_ESCC$sample]

dat.final<-cbind(exp_control,exp_tumor)
## 92  8975
write.table(dat.final,file = 'dat.tcga.xls',sep = '\t',quote = F,row.names = T)
##fpkm
expr_fpkm<-read_tsv(file = 'TCGA-ESCA.htseq_fpkm.tsv.gz')
expr_fpkm<-as.data.frame(expr_fpkm)
expr_fpkm<-column_to_rownames(expr_fpkm,var = "Ensembl_ID")
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
dat_fpkm<-dat_fpkm[,colnames(dat.final)]
dat_fpkm<-dat_fpkm[rownames(dat.final),]
write.table(dat_fpkm,file = 'dat.fpkm.xls',sep = '\t',row.names = T,quote = F)
# fpkm转TPM
FPKM2TPM <- function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

dat_tpm <- apply(dat_fpkm,2,FPKM2TPM)
write.table(dat_tpm,file = 'dat.tpm.xls',sep = '\t',row.names = T,quote = F)

### TCGAall
tcga.all<-tcga.expr
## 8798
tcga.all$ID <- rownames(tcga.all)
tcga.all$ID<-as.character(tcga.all$ID)
probe2symbol$ID<-as.character(probe2symbol$ID)
tcga.all<-tcga.all %>%
  inner_join(probe2symbol,by='ID')%>% 
  select(-ID)%>%     ## 去除多余信息
  select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('TCGA',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
dim(tcga.all)
tcga.all<-tcga.all[,colnames(dat.final)]
write.table(tcga.all,file = 'dat.all(fpkm).xls',sep = '\t',row.names = T,quote = F)

library(GEOquery)
library(Biobase)
gset<-getGEO("GSE53625",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
a=gset[[1]]
library(idmap3)
probe2symobl<-idmap3::get_pipe_IDs('GPL18109')
colnames(probe2symobl)<-c('ID','symbol')
length(unique(probe2symobl$symbol))
#probe2symobl=probe2symobl[probe2symobl$symbol!='',]
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
#dat.lnc<-dat[rownames(dat)%in%rownames(dat.final),]
dat.lnc<-dat[rownames(dat)%in%lncRNA$gene_name,]
write.table(dat,file = 'dat.geo.xls',sep = '\t',row.names = T,quote = F)
write.table(dat.lnc,file = 'dat.geo(lnc).xls',sep = '\t',row.names = T,quote = F)

dat<-dat.lnc
pd<-pData(a)
survival<-data.frame(sample=pd$geo_accession,group=pd$title,OS=pd$characteristics_ch1.14,OS.time=pd$characteristics_ch1.15)
survival<-survival[order(survival$group),]
survival<-survival[c(1:179),]
table(survival$OS)
survival$OS<-ifelse(survival$OS=='death at fu: no',0,1)
table(survival$OS.time)
survival$OS.time<-gsub('survival time(months): ','',survival$OS.time,fixed = T)
survival$OS.time<-as.numeric(survival$OS.time)*30
survival<-survival[,-2]
write.table(survival,file = 'survival.xls',sep = '\t',row.names = F,quote = F)
group<-data.frame(sample=pd$geo_accession,group=pd$title)
group<-group[order(group$group),]
group$group<-c(rep('ESCC',179),rep('control',179))
## 把2个数据集拆开
group<-group[order(group$sample),]
group1<-group[-c(1:120),]
group2<-group[c(1:120),]
write.table(group,file = 'group(geo).xls',sep = '\t',quote = F,row.names = F)
write.table(group1,file = 'group1.xls',sep = '\t',quote = F,row.names = F)
write.table(group2,file = 'group2.xls',sep = '\t',quote = F,row.names = F)
dat1<-dat[,group1$sample]
dat2<-dat[,group2$sample]
write.table(dat1,file = 'dat1.xls',sep = '\t',row.names = T,quote = F)
write.table(dat2,file = 'dat2.xls',sep = '\t',row.names = T,quote = F)
phenotype<-data.frame(sample=pd$geo_accession,
                      gender=pd$characteristics_ch1.2,
                      T.stage=pd$characteristics_ch1.7,
                      N.stage=pd$characteristics_ch1.8,
                      TNM.stage=pd$characteristics_ch1.9)
#phenotype<-phenotype[phenotype$sample%in%colnames(dat2),]
write.table(phenotype,file = 'phenotype.xls',sep = '\t',row.names = F,quote = F)


