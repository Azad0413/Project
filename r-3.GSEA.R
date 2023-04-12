rm(list = ls())
# GSEA-----------
setwd("/data/nas1/luchunlin/project/LZZK-503")
if (! dir.exists("./03_GSEA")){
  dir.create("./03_GSEA")
}
setwd("./03_GSEA")
dat<-read.delim2('/data/nas1/luchunlin/project/LLZK-505/00_rawdata/dat.tcga.xls',row.names = 1)%>%lc.tableToNum()
colnames<-data.frame(sample=colnames(dat))
colnames$sample<-gsub('.','-',colnames$sample,fixed = T)
colnames(dat)<-colnames$sample
dat<-round(dat,digits = 0)
## 05-1 GRIM-19(NDUFA13)-------
group<-read.delim2('/data/nas1/luchunlin/project/LZZK-503/01_DEGs/group.xls')
library(DESeq2)
colData<-data.frame(sample=group$sample,
                    group=group$group)
table(colData$group)
colData$group<-factor(colData$group,levels = c('High','Low'))
dds<-DESeqDataSetFromMatrix(countData = dat,colData=colData,design = ~group)
dds = dds[rownames(counts(dds)) > 1,]
dds<-DESeq(dds)
dds
res =results(dds, contrast = c("group","Low","High"))
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
write.table(kegg_result,file = 'GSEA(ZEB2-AS1).xls',sep = '\t',quote = F,row.names = F)
gseaplot2(kegg_gsea,c(1:5),color = c('#7B68EE','#CD3333','#20B2AA','#FF8C00','#FF6666'),
          title = 'KEGG GSEA(ZEB2-AS1)',
          base_size = 11,
          rel_heights = c(1.5, 0.3, 0.5))
