## 05 GSEA-----
rm(list = ls())
setwd("/data/nas1/luchunlin/project/LLZK-505")
if (! dir.exists("./04_GSEA")){
  dir.create("./04_GSEA")
}
setwd("./04_GSEA")
dat<-read.delim2('/data/nas1/luchunlin/project/LLZK-505/00_rawdata/dat.tcga.xls',row.names = 1)%>%lc.tableToNum()
colnames<-data.frame(sample=colnames(dat))
colnames$sample<-gsub('.','-',colnames$sample,fixed = T)
colnames(dat)<-colnames$sample
dat<-round(dat,digits = 0)
## 05-1 GRIM-19(NDUFA13)-------
group1<-read.delim2('/data/nas1/luchunlin/project/LLZK-505/02_clinical/group(GRIM-19).xls')
library(DESeq2)
colData<-data.frame(sample=group1$sample,
                    group=group1$group)
colData$group<-ifelse(colData$group=='High GRIM19','High_GRIM19','Low_GRIM19')
table(colData$group)
colData$group<-factor(colData$group,levels = c('High_GRIM19','Low_GRIM19'))
dds<-DESeqDataSetFromMatrix(countData = dat,colData=colData,design = ~group)
dds = dds[rownames(counts(dds)) > 1,]
dds<-DESeq(dds)
dds
res =results(dds, contrast = c("group","High_GRIM19","Low_GRIM19"))
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
write.table(kegg_result,file = 'GSEA(GRIM-19).xls',sep = '\t',quote = F,row.names = F)
gseaplot2(kegg_gsea,c(1:5),color = c('#7B68EE','#CD3333','#20B2AA','#FF8C00','#FF6666'),
          title = 'KEGG GSEA(GRIM-19)',
          base_size = 11,
          rel_heights = c(1.5, 0.3, 0.5))
## 05-2 NDUFS3------
group2<-read.delim2('/data/nas1/luchunlin/project/LLZK-505/02_clinical/group(NDUFS3).xls')
colData<-data.frame(sample=group2$sample,
                    group=group2$group)
colData$group<-ifelse(colData$group=='High NDUFS3','High_NDUFS3','Low_NDUFS3')
table(colData$group)
colData$group<-factor(colData$group,levels = c('High_NDUFS3','Low_NDUFS3'))
dds<-DESeqDataSetFromMatrix(countData = dat,colData=colData,design = ~group)
dds = dds[rownames(counts(dds)) > 1,]
dds<-DESeq(dds)
dds
res =results(dds, contrast = c("group","High_NDUFS3","Low_NDUFS3"))
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
set.seed(1)
kegg_gsea <- GSEA(geneList, TERM2GENE = kegg_set, pvalueCutoff = 0.05,eps = 0)
kegg_result <- kegg_gsea@result
dim(kegg_result)
write.table(kegg_result,file = 'GSEA(NDUFS3).xls',sep = '\t',quote = F,row.names = F)
gseaplot2(kegg_gsea,c(1:5),color = c('#7B68EE','#CD3333','#20B2AA','#FF8C00','#FF6666'),
          title = 'KEGG GSEA(NDUFS3)',
          base_size = 11,
          rel_heights = c(1.5, 0.3, 0.5))
## 05-3 NDUFA4------
group3<-read.delim2('/data/nas1/luchunlin/project/LLZK-505/02_clinical/group(NDUFA4).xls')
colData<-data.frame(sample=group3$sample,
                    group=group3$group)
colData$group<-ifelse(colData$group=='High NDUFA4','High_NDUFA4','Low_NDUFA4')
table(colData$group)
colData$group<-factor(colData$group,levels = c('High_NDUFA4','Low_NDUFA4'))
dds<-DESeqDataSetFromMatrix(countData = dat,colData=colData,design = ~group)
dds = dds[rownames(counts(dds)) > 1,]
dds<-DESeq(dds)
dds
res =results(dds, contrast = c("group","High_NDUFA4","Low_NDUFA4"))
res =res[order(res$padj),]
head(res)
summary(res)
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
set.seed(1)
kegg_gsea <- GSEA(geneList, TERM2GENE = kegg_set, pvalueCutoff = 0.05,eps = 0)
kegg_result <- kegg_gsea@result
dim(kegg_result)
write.table(kegg_result,file = 'GSEA(NDUFA4).xls',sep = '\t',quote = F,row.names = F)
gseaplot2(kegg_gsea,c(1:5),color = c('#7B68EE','#CD3333','#20B2AA','#FF8C00','#FF6666'),
          title = 'KEGG GSEA(NDUFA4)',
          base_size = 11,
          rel_heights = c(1.5, 0.3, 0.5))
## 05-4 LRPPRC------
group4<-read.delim2('/data/nas1/luchunlin/project/LLZK-505/02_clinical/group(LRPPRC).xls')
colData<-data.frame(sample=group4$sample,
                    group=group4$group)
colData$group<-ifelse(colData$group=='High LRPPRC','High_LRPPRC','Low_LRPPRC')
table(colData$group)
colData$group<-factor(colData$group,levels = c('High_LRPPRC','Low_LRPPRC'))
dds<-DESeqDataSetFromMatrix(countData = dat,colData=colData,design = ~group)
dds = dds[rownames(counts(dds)) > 1,]
dds<-DESeq(dds)
dds
res =results(dds, contrast = c("group","High_LRPPRC","Low_LRPPRC"))
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
set.seed(1)
kegg_gsea <- GSEA(geneList, TERM2GENE = kegg_set, pvalueCutoff = 0.05)
kegg_result <- kegg_gsea@result
dim(kegg_result)
write.table(kegg_result,file = 'GSEA(LRPPRC).xls',sep = '\t',quote = F,row.names = F)
gseaplot2(kegg_gsea,c(1:5),color = c('#7B68EE','#CD3333','#20B2AA','#FF8C00','#FF6666'),
          title = 'KEGG GSEA(LRPPRC)',
          base_size = 11,
          rel_heights = c(1.5, 0.3, 0.5))
