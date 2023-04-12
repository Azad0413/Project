rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-300-8/")
if (! dir.exists("./02_DEGs(GSE99671)")){
  dir.create("./02_DEGs(GSE99671)")
}
setwd("./02_DEGs(GSE99671)")
library(lance)
dat<-read.delim2('../00_rawdata/dat(GSE99671).xls',row.names = 1)%>%lc.tableToNum()
library(DESeq2)
colData <- read.delim2('../00_rawdata/group(GSE99671).xls')
table(colData$group)
colData$group <- ifelse(colData$group=='control','Normal','Tumor')
colData$group<-factor(colData$group,levels = c('Normal','Tumor'))
dds<-DESeqDataSetFromMatrix(countData = dat,colData=colData,design = ~group)
dds = dds[rownames(counts(dds)) > 1,]
dds<-estimateSizeFactors(dds)
##提取标准化后的数据 
#normalized_counts <- counts(dds,normalized=T) 
#write.table(normalized_counts,file = 'normalized.counts.xls',sep = '\t',row.names = T,quote = F)
dds<-DESeq(dds)
## 提取差异结果
res =results(dds, contrast = c("group","Tumor","Normal"))
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
## 637 16862   265 
sig_diff <- subset(DEG,
                   DEG$padj < 0.05 & abs(DEG$log2FoldChange) >= logFC_cutoff)
## 902
DEG_write <- cbind(GeneSymbol=rownames(DEG), DEG)
write.table(DEG_write, file = "DEG_all(GSE99671).xls",
            quote = F,
            sep = "\t",
            row.names = F)
sig_diff_write <- cbind(GeneSymbol=rownames(sig_diff), sig_diff)
write.table(sig_diff_write, file = "DEG_sig(GSE99671).xls",
            quote = F,
            sep = "\t",
            row.names = F)

