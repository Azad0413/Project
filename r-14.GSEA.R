rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-406-12/")
if (! dir.exists("./14_GSEA")){
  dir.create("./14_GSEA")
}
setwd("./14_GSEA")

dat.tcga<-read.delim2("../00_rawdata/dat.tcga.xls", row.names = 1)%>% lc.tableToNum
colnames(dat.tcga)<-gsub('.','-',colnames(dat.tcga),fixed = T)
dat.tcga<-round(dat.tcga,digits = 0)
library(GSVA)
library(GSEABase)
library(DESeq2)
group<-read.delim2('../09_risk/risk.xls')%>%dplyr::select(c('id','risk'))
#BiocManager::install('DESeq2')
group<-group[order(group$risk),]
group$risk<-ifelse(group$risk=='0','High_risk','Low_risk')
gsea_exp<-dat.tcga[,group$id]
library(DESeq2)
colData<-group
colnames(colData)<-c('sample','group')
colData$group<-factor(colData$group,levels = c('High_risk','Low_risk'))
dds<-DESeqDataSetFromMatrix(countData = gsea_exp,colData=colData,design = ~group)
dds = dds[rownames(counts(dds)) > 1,]
dds<-DESeq(dds)
dds
res =results(dds, contrast = c("group","High_risk","Low_risk"))
res =res[order(res$padj),]
head(res)
summary(res)
table(res$padj<0.05)
allGeneSets<-as.data.frame(res)
allGeneSets<-na.omit(allGeneSets)
logFCcutoff <- 0
allGeneSets$change = as.factor(
  ifelse(allGeneSets$padj < 0.05 & abs(allGeneSets$log2FoldChange) > logFCcutoff,
         ifelse(allGeneSets$log2FoldChange > logFCcutoff,'UP','DOWN'),'NOT')
)

DEGeneSets <- subset(allGeneSets,
                     allGeneSets$padj < 0.05 & abs(allGeneSets$log2FoldChange) > 0)
genelist <- allGeneSets$log2FoldChange
names(genelist) <- rownames(allGeneSets)
geneList <- sort(genelist, decreasing = T)
DEGeneSets <- DEGeneSets[order(DEGeneSets$padj),]
dim(DEGeneSets)
write.table(DEGeneSets,'DEGeneSet.xls',sep = '\t',row.names = T,quote = F)
write.table(allGeneSets,'allGeneSet.xls',sep = '\t',row.names = T,quote = F)
## GSEA KEGG-----
library(clusterProfiler)
library(enrichplot)
kegg_set<- read.gmt("/data/nas1/luchunlin/pipeline/GSVA/c2.cp.kegg.v7.4.symbols.gmt")
set.seed(1)
kegg_gsea <- GSEA(geneList, TERM2GENE = kegg_set, pvalueCutoff = 0.05)
kegg_result <- kegg_gsea@result
write.table(kegg_result,file = 'GSEA.result.xls',sep = '\t',row.names = F,quote = F)
dim(kegg_result)
pdf(file = '01.GSEA(KEGG).pdf',w=9,h=6)
gseaplot2(kegg_gsea,c(1:5),color = c('#7B68EE','#CD3333','#20B2AA','#FF8C00','#FF6666'),
          title = 'KEGG GSEA',
          base_size = 11,
          rel_heights = c(1.5, 0.3, 0.5))
dev.off()
png(file = '01.GSEA(KEGG).png',w=600,h=400)
gseaplot2(kegg_gsea,c(1:5),color = c('#7B68EE','#CD3333','#20B2AA','#FF8C00','#FF6666'),
          title = 'KEGG GSEA',
          base_size = 11,
          rel_heights = c(1.5, 0.3, 0.5))
dev.off()
