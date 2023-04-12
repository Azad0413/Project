rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/BJTC-386-10/")
if (! dir.exists("./00_rawdata")){
  dir.create("./00_rawdata")
}
setwd("./00_rawdata")
##01-1 TCGA KIRC-------
library(TCGAbiolinks)
library(readr)
library(readxl)
library(tidyverse)
## 读取从xena下载的数据
tcga.expr<-read_tsv(file = '/data/nas1/luchunlin/TCGA.matrix/TCGA-KIRC.htseq_counts.tsv.gz')
tcga.expr<-as.data.frame(tcga.expr)
tcga.expr<-column_to_rownames(tcga.expr,var = "Ensembl_ID")
## xena下载的数据经过了log2+1转化，需要将其还原
tcga.expr<-2^tcga.expr-1
## 加载注释文件
library("rtracklayer")
gtf_data = import('/data/nas1/luchunlin/pipeline/GENEANNO/gencode.v22.annotation.gtf.gz') #gtf的路径
gtf_data = as.data.frame(gtf_data)
table(gtf_data$gene_type)
lncRNA=gtf_data%>%
  dplyr::filter(type=="gene",gene_type=="lincRNA")%>%
  dplyr::select(gene_id,gene_type,gene_name)
write.table(lncRNA,file = 'lncRNA.xls',sep = '\t',row.names = F,quote = F)

## 对数据进行id转化
genecode<-read.table(file = '/data/nas1/luchunlin/pipeline/GENEANNO/gencode.v22.annotation.gene.probeMap')
probe2symbol<-genecode[,(1:2)]
colnames(probe2symbol)<-c('ID','symbol')
probe2symbol<-probe2symbol[-1,]
dat.tcga<-tcga.expr
dat.tcga<-tcga.expr[rownames(tcga.expr)%in%lncRNA$gene_id,]
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
#7648  607
keep<-rowSums(dat.tcga>0)>=floor(0.50*ncol(dat.tcga))
dat.tcga<-dat.tcga[keep,]
## 3536
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
## 535
exp_control<-dat.tcga[,which(!colnames(dat.tcga)%in%mete$id)]
exp_control<-as.data.frame(exp_control)
# 72
## 保留有生存数据的
survival<-read.delim2('/data/nas1/luchunlin/TCGA_survival/TCGA-KIRC.survival.tsv')
exp_tumor<-exp_tumor[,colnames(exp_tumor)%in%survival$sample]
##531
exp_control<-dat.tcga[,which(!colnames(dat.tcga)%in%mete$id)]
exp_control<-as.data.frame(exp_control)
# 72
dat.final<-cbind(exp_control,exp_tumor)
## 603 3536
write.table(dat.final,file = 'dat.tcga.xls',sep = '\t',quote = F,row.names = T)
##fpkm
expr_fpkm<-read_tsv(file = '/data/nas1/luchunlin/TCGA.matrix/TCGA-KIRC.htseq_fpkm.tsv.gz')
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

# ### TCGAall
tcga.all<-expr_fpkm
## 8798
tcga.all$ID <- rownames(tcga.all)
tcga.all$ID<-as.character(tcga.all$ID)
probe2symbol$ID<-as.character(probe2symbol$ID)
tcga.all<-tcga.all %>%
  inner_join(probe2symbol,by='ID')%>%
  dplyr::select(-ID)%>%     ## 去除多余信息
  dplyr::select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('TCGA',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
dim(tcga.all)
tcga.all<-tcga.all[,colnames(dat.final)]
write.table(tcga.all,file = 'dat.all(fpkm).xls',sep = '\t',row.names = T,quote = F)
                                                                                                                                                              
## 验证集ICGC--------
library(lance)
dat.va<-read.delim2('Merge_RNAseq_RPKM_ENSG.txt',header = T)%>%column_to_rownames(var = 'gene_id')%>%lc.tableToNum()
# dat.va<-edgeR::cpm(dat.va)
group.va<-read.delim2('Merge_clinical.txt')
group.va <- data.frame(sample=group.va$icgc_sample_id,group=group.va$specimen_type)
table(group.va$group)
group.va$group <- ifelse(group.va$group=='Primary tumour - solid tissue','Tumor','Normal')
group.va <- group.va[group.va$sample%in%colnames(dat.va),]
group.va <- group.va[order(group.va$group),]
table(group.va$group)
# control   tumor 
# 45      91 
dat.va <- dat.va[,group.va$sample]
hubgene <- read.delim2('../07_hubgene/hubgene.xls')
hubid <- probe2symbol[probe2symbol$symbol%in%hubgene$symbol,]
hubid$ID <- substr(hubid$ID,1,15)
dat.va <- dat.va[hubid$ID,]
rownames(dat.va) <- hubid$symbol
write.table(dat.va,file = 'dat.va.xls',sep = '\t',row.names = T,quote = F)
write.table(group.va,file = 'group.va.xls',sep = '\t',row.names = F,quote = F)
