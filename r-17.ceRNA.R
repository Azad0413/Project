rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-321")
if (! dir.exists("./17_ceRNA")){
  dir.create("./17_ceRNA")
}
setwd("./17_ceRNA")
library(GEOquery)
library(Biobase)
library(tidyverse)
library(dplyr)
library(lance)
library(readxl)
expr<-read_xlsx('GSE195599_gene_expression_anno.xlsx')
colnames(expr)

library("rtracklayer")
gtf_data = import('gencode.v40.annotation.gtf') #gtf的路径
gtf_data = as.data.frame(gtf_data)
#提取表达谱内的lncRNA,记得先把表达谱内的Ensembl_ID改成"gene_id")
table(gtf_data$gene_type)
lncRNA=gtf_data%>%
 dplyr::filter(type=="gene",gene_type=="lncRNA")%>%
 dplyr::select(gene_id,gene_type,gene_name)
colnames(lncRNA)
e2s<-lncRNA%>%
  separate("gene_id",into = 'ensemble',sep = '\\.')

dat<-expr[expr$gene_id%in%e2s$ensemble,]
colnames(dat)
dat<-dat[,c(2,18:24)]%>%as.data.frame()
dat<-dat[!duplicated(dat$gene_name),]
rownames(dat)<-dat$gene_name
dat<-dat[,-1]
colnames(dat)
group<-data.frame(sample=colnames(dat),group=c(rep('asthma',3),rep('control',2),rep('asthma',1),rep('control',1)))
colData<-group
colData$group<-factor(colData$group,levels = c('control','asthma'))
library(DESeq2)
dds<-DESeqDataSetFromMatrix(countData = dat,colData=colData,design = ~group)
dds = dds[rownames(counts(dds)) > 1,]
dds<-estimateSizeFactors(dds)
##提取标准化后的数据 
#normalized_counts <- counts(dds,normalized=T) 
dds<-DESeq(dds)
#write.table(normalized_counts,file = 'normalized.counts.xls',sep = '\t',row.names = T,quote = F)
## 提取差异结果
## deseq2默认的是后一个比前一个
res =results(dds, contrast = c("group","control","asthma"))
res =res[order(res$padj),]
head(res)
summary(res)
table(res$padj<0.05)
DEG <- subset(res, padj < 0.05 & abs(log2FoldChange) >1 )
DEG<-as.data.frame(res)
DEG<-na.omit(DEG)
dim(DEG)
head(DEG)
## 添加change列
logFC_cutoff<-1
DEG$change=as.factor(
  ifelse(DEG$padj<0.05&abs(DEG$log2FoldChange)>logFC_cutoff,
         ifelse(DEG$log2FoldChange>logFC_cutoff,'UP','DOWN'),'NOT'))
table(DEG$change)
## DOWN   NOT    UP 
## 231 3614  229
sig_diff <- subset(DEG,
                   DEG$padj < 0.05 & abs(DEG$log2FoldChange) >= logFC_cutoff)
## 460
DEG_write <- cbind(GeneSymbol=rownames(DEG), DEG)
write.table(DEG_write, file = "DEG_all.xls",
            quote = F,
            sep = "\t",
            row.names = F)
sig_diff_write <- cbind(GeneSymbol=rownames(sig_diff), sig_diff)
write.table(sig_diff_write, file = "DEG_sig.xls",
            quote = F,
            sep = "\t",
            row.names = F)


gset<-getGEO("GSE143192",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
a=gset[[1]]
pd<-pData(a)
gpl<-getGEO("GPL22120",destdir = '.')
gpl<-Table(gpl)    
colnames(gpl)
table(gpl$TYPE)
gpl<-gpl[which(gpl$TYPE=='lncRNA'),]
probe2symbol<-gpl %>%
  dplyr::select('ID','GENE_SYMBOL')%>%
  filter('GENE_SYMBOL'!='')%>%
  separate('GENE_SYMBOL',c('symbol','drop'),sep = '///')%>%
  dplyr::select(-drop)
probe2symbol=probe2symbol[probe2symbol$symbol!='',]
probe2symbol=probe2symbol[probe2symbol$symbol!='-',]
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

library(magrittr)
library(stringr)
library(lance)
library(DESeq2)
library(limma)
##02-1 差异基因鉴定
df = dat
df.group = data.frame(sample=pd$geo_accession,group=pd$`diagnosis:ch1`)
table(df.group$group)
df.group$group<-ifelse(df.group$group=='healthy','control','asthma')
df = df[df.group$sample]
df.group$group = factor(df.group$group, levels = c("control", "asthma"))
design.mat = cbind(control = ifelse(df.group$group == "control", 1, 0), 
                   asthma = ifelse(df.group$group == "control", 0, 1))
contrast.mat = makeContrasts(contrasts="asthma-control", levels=design.mat)

fit = lmFit(df, design.mat)
fit = contrasts.fit(fit, contrast.mat)
fit = eBayes(fit)
fit = topTable(fit, coef = 1, number = Inf, adjust.method = "fdr")
#fit = fit[c(1,4,5)]
DEG=na.omit(fit)
logFC_cutoff <- 0
DEG$change = as.factor(
  ifelse(DEG$adj.P.Val <0.05 & abs(DEG$logFC) > logFC_cutoff,
         ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
sig_diff <- subset(DEG,
                   DEG$adj.P.Val < 0.05 & abs(DEG$logFC) > logFC_cutoff)

dim(DEG)
dim(sig_diff)
## 3630    6
summary(sig_diff$change)
# DOWN  NOT   UP 
# 1470    0 2160 
write.table(DEG,file = "DEG_all.xls",quote = F,sep = "\t",row.names = T)
write.table(sig_diff,file = "DEG_sig.xls",quote = F,sep = "\t",row.names = T)
