rm(list = ls())
setwd("/data/nas1/luchunlin/project/LLZK-505")
if (! dir.exists("./10_GSEA(deg)")){
  dir.create("./10_GSEA(deg)")
}
setwd("./10_GSEA(deg)")

dat<-read.delim2('/data/nas1/luchunlin/project/LLZK-505/00_rawdata/dat.tcga.xls',row.names = 1)%>%lc.tableToNum()
colnames<-data.frame(sample=colnames(dat))
colnames$sample<-gsub('.','-',colnames$sample,fixed = T)
colnames(dat)<-colnames$sample
dat<-round(dat,digits = 0)
## 11-1 GRIM-19  NDUFS3------
cluster1<-read.delim2('/data/nas1/luchunlin/project/LLZK-505/05_Consensus/cluster(GRIM-19&NDUFS3).xls')
cluster1$cluster<-ifelse(cluster1$cluster=='1','Cluster1','Cluster2')
cluster1$sample<-rownames(cluster1)

library(DESeq2)
colData<-data.frame(sample=cluster1$sample,group=cluster1$cluster)
table(colData$group)
colData$group<-factor(colData$group,levels = c('Cluster1','Cluster2'))
dds<-DESeqDataSetFromMatrix(countData = dat,colData=colData,design = ~group)
dds = dds[rownames(counts(dds)) > 1,]
dds<-DESeq(dds)
dds
res =results(dds, contrast = c("group","Cluster1","Cluster2"))
res =res[order(res$padj),]
head(res)
summary(res)
table(res$padj<0.05)
allGeneSets<-as.data.frame(res)
allGeneSets<-na.omit(allGeneSets)
logFCcutoff <- 1
allGeneSets$change = as.factor(
  ifelse(allGeneSets$padj < 0.05 & abs(allGeneSets$log2FoldChange) > logFCcutoff,
         ifelse(allGeneSets$log2FoldChange > logFCcutoff,'UP','DOWN'),'NOT')
)

DEGeneSets <- subset(allGeneSets,
                     allGeneSets$padj < 0.05 & abs(allGeneSets$log2FoldChange) > logFCcutoff)
genelist <- allGeneSets$log2FoldChange
names(genelist) <- rownames(allGeneSets)
geneList <- sort(genelist, decreasing = T)
DEGeneSets <- DEGeneSets[order(DEGeneSets$padj),]
dim(DEGeneSets)
## GSEA KEGG
library(clusterProfiler)
library(enrichplot)
kegg_set<- read.gmt("c2.cp.kegg.v7.5.1.symbols.gmt")
set.seed(1)
kegg_gsea <- GSEA(geneList, TERM2GENE = kegg_set, pvalueCutoff = 0.05)
kegg_result <- kegg_gsea@result
dim(kegg_result)
write.table(kegg_result,file = 'GSEA(GRIM19&NDUFS3).xls',sep = '\t',quote = F,row.names = F)
gseaplot2(kegg_gsea,c(1:5),color = c('#7B68EE','#CD3333','#20B2AA','#FF8C00','#FF6666'),
          title = 'KEGG GSEA(GRIM19&NDUFS3)',
          base_size = 11,
          rel_heights = c(1.5, 0.3, 0.5))

## 11-2 NDUFA4,LRPPRC-----
cluster2<-read.delim2('/data/nas1/luchunlin/project/LLZK-505/05_Consensus/cluster(NDUFA4&LRPPRC).xls')
cluster2$cluster<-ifelse(cluster2$cluster=='1','Cluster1','Cluster2')
cluster2$sample<-rownames(cluster2)
colData<-data.frame(sample=cluster2$sample,group=cluster2$cluster)
table(colData$group)
colData$group<-factor(colData$group,levels = c('Cluster1','Cluster2'))
dds<-DESeqDataSetFromMatrix(countData = dat,colData=colData,design = ~group)
dds = dds[rownames(counts(dds)) > 1,]
dds<-DESeq(dds)
dds
res =results(dds, contrast = c("group","Cluster1","Cluster2"))
res =res[order(res$padj),]
head(res)
summary(res)
table(res$padj<0.05)
allGeneSets<-as.data.frame(res)
allGeneSets<-na.omit(allGeneSets)
logFCcutoff <- 1
allGeneSets$change = as.factor(
  ifelse(allGeneSets$padj < 0.05 & abs(allGeneSets$log2FoldChange) > logFCcutoff,
         ifelse(allGeneSets$log2FoldChange > logFCcutoff,'UP','DOWN'),'NOT')
)

DEGeneSets <- subset(allGeneSets,
                     allGeneSets$padj < 0.05 & abs(allGeneSets$log2FoldChange) > logFCcutoff)
genelist <- allGeneSets$log2FoldChange
names(genelist) <- rownames(allGeneSets)
geneList <- sort(genelist, decreasing = T)
DEGeneSets <- DEGeneSets[order(DEGeneSets$padj),]
dim(DEGeneSets)
## GSEA KEGG
library(clusterProfiler)
library(enrichplot)
kegg_set<- read.gmt("c2.cp.kegg.v7.5.1.symbols.gmt")
set.seed(1)
kegg_gsea <- GSEA(geneList, TERM2GENE = kegg_set, pvalueCutoff = 0.05,eps = 0)
kegg_result <- kegg_gsea@result
dim(kegg_result)
write.table(kegg_result,file = 'GSEA(NDUFA4&LRPPRC).xls',sep = '\t',quote = F,row.names = F)
gseaplot2(kegg_gsea,c(1:5),color = c('#7B68EE','#CD3333','#20B2AA','#FF8C00','#FF6666'),
          title = 'KEGG GSEA(NDUFA4&LRPPRC)',
          base_size = 11,
          rel_heights = c(1.5, 0.3, 0.5))
