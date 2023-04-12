rm(list = ls())
setwd("/data/nas1/luchunlin/project/NN-0118-2/")
if (! dir.exists("./11_GSEA")){
  dir.create("./11_GSEA")
}
setwd("./11_GSEA")

dat<-read.delim2('../00_rawdata/dat.tpm.xls',row.names = 1)%>%lc.tableToNum()
colnames(dat)<-gsub('.','-',colnames(dat),fixed = T)
group <- read.delim2('../01_DEGs(TCGA)/group.xls')
tumor.sample <- group[which(group$group=='Tumor'),]
dat <- dat[,tumor.sample$sample]
dat<-round(dat,digits = 0)
## 1 CCNB1-------
colData<-read.delim2('../10_KM/group(CCNB1).xls')
library(DESeq2)
colData$group<-ifelse(colData$group=='High CCNB1','High_CCNB1','Low_CCNB1')
table(colData$group)
colData$group<-factor(colData$group,levels = c('High_CCNB1','Low_CCNB1'))
dds<-DESeqDataSetFromMatrix(countData = dat,colData=colData,design = ~group)
dds = dds[rownames(counts(dds)) > 1,]
dds<-DESeq(dds)
dds
res =results(dds, contrast = c("group","High_CCNB1","Low_CCNB1"))
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
kegg_set<- read.gmt("/data/nas1/luchunlin/pipeline/GSVA/c2.cp.kegg.v7.5.1.symbols.gmt")
set.seed(1)
kegg_gsea <- GSEA(geneList, TERM2GENE = kegg_set, pvalueCutoff = 0.05)
kegg_result <- kegg_gsea@result
dim(kegg_result)
write.table(kegg_result,file = 'GSEA(CCNB1).xls',sep = '\t',quote = F,row.names = F)
pdf(file = '01.GSEA.pdf',w=9,h=5)
gseaplot2(kegg_gsea,c(1:5),color = c('#7B68EE','#CD3333','#20B2AA','#FF8C00','#FF6666'),
          title = 'KEGG GSEA(CCNB1)',
          base_size = 11,
          rel_heights = c(1.5, 0.3, 0.5))
dev.off()

png(file = '01.GSEA.png',w=9,h=5)
gseaplot2(kegg_gsea,c(1:5),color = c('#7B68EE','#CD3333','#20B2AA','#FF8C00','#FF6666'),
          title = 'KEGG GSEA(CCNB1)',
          base_size = 11,
          rel_heights = c(1.5, 0.3, 0.5))
dev.off()
