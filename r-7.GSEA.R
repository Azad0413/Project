rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-441-3/")
if (! dir.exists("./07_GSEA")){
  dir.create("./07_GSEA")
}
setwd("./07_GSEA")
dat<-read.csv('../00_rawdata/dat.tcga.xls',sep = '\t',row.names = 1)
colnames(dat) <- gsub('.','-',colnames(dat),fixed = T)
hub_exp<-dat['UGCG',]
group <- read.delim2('../05_survival/group(UGCG).xls')
dat <- round(dat,digits = 0)
library(DESeq2)
colData<-data.frame(sample=group$sample,
                    group=group$group)
table(colData$group)
colData$group<-ifelse(colData$group=='High UGCG','High_UGCG','Low_UGCG')
table(colData$group)
colData$group<-factor(colData$group,levels = c('High_UGCG','Low_UGCG'))
dds<-DESeqDataSetFromMatrix(countData = dat,colData=colData,design = ~group)
dds = dds[rownames(counts(dds)) > 1,]
dds<-DESeq(dds)
dds
res =results(dds, contrast = c("group","High_UGCG","Low_UGCG"))
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
set.seed(8)
kegg_gsea <- GSEA(geneList, TERM2GENE = kegg_set, pvalueCutoff = 1)
kegg_result <- kegg_gsea@result
kegg_result <- kegg_result[which(kegg_result$pvalue<0.05),]
dim(kegg_result)
write.table(kegg_result,file = 'GSEA(KEGG).xls',sep = '\t',quote = F,row.names = F)
pdf(file = '01.GSEA(KEGG).pdf',w=9,h=6)
gseaplot2(kegg_gsea,c(4,5,7,9,12),color = c('#7B68EE','#CD3333','#20B2AA','#FF8C00','#FF6666'),
          title = 'KEGG GSEA',
          base_size = 11,
          rel_heights = c(1.5, 0.3, 0.5))
dev.off()
png(file = '01.GSEA(KEGG).png',w=600,h=400)
gseaplot2(kegg_gsea,c(4,5,7,9,12),color = c('#7B68EE','#CD3333','#20B2AA','#FF8C00','#FF6666'),
          title = 'KEGG GSEA',
          base_size = 11,
          rel_heights = c(1.5, 0.3, 0.5))
dev.off()


## GSEA GOBP
library(clusterProfiler)
library(enrichplot)
GOBP_set<- read.gmt("/data/nas1/luchunlin/pipeline/GSVA/c5.go.bp.v7.4.symbols.gmt")
set.seed(1)
GOBP_gsea <- GSEA(geneList, TERM2GENE = GOBP_set, pvalueCutoff = 0.05)
GOBP_result <- GOBP_gsea@result
dim(GOBP_result)
write.table(GOBP_result,file = 'GSEA(GOBP).xls',sep = '\t',quote = F,row.names = F)
pdf(file = '02.GSEA(GOBP).pdf',w=9,h=6)
gseaplot2(GOBP_gsea,c(1:5),color = c('#7B68EE','#CD3333','#20B2AA','#FF8C00','#FF6666'),
          title = 'GOBP GSEA',
          base_size = 11,
          rel_heights = c(1.5, 0.3, 0.5))
dev.off()
png(file = '02.GSEA(GOBP).png',w=600,h=400)
gseaplot2(GOBP_gsea,c(1:5),color = c('#7B68EE','#CD3333','#20B2AA','#FF8C00','#FF6666'),
          title = 'GOBP GSEA',
          base_size = 11,
          rel_heights = c(1.5, 0.3, 0.5))
dev.off()


##hallmark-------
hallmark_set<- read.gmt("/data/nas1/luchunlin/pipeline/GSVA/h.all.v2022.1.Hs.symbols.gmt")
set.seed(8)
hallmark_gsea <- GSEA(geneList, TERM2GENE = hallmark_set, pvalueCutoff = 1)
hallmark_result <- hallmark_gsea@result
hallmark_result <- hallmark_result[which(hallmark_result$pvalue<0.05),]
dim(hallmark_result)
write.table(hallmark_result,file = 'GSEA(hallmark).xls',sep = '\t',quote = F,row.names = F)
pdf(file = '03.GSEA(hallmark).pdf',w=9,h=6)
gseaplot2(hallmark_gsea,c(1,2,6,8,10:12),color = c('#7B68EE','#CD3333','#20B2AA','#FF8C00','#FF6666','darkblue','green'),
          title = 'hallmark GSEA',
          base_size = 11,
          rel_heights = c(1.5, 0.3, 0.5))
dev.off()
png(file = '03.GSEA(hallmark).png',w=600,h=400)
gseaplot2(hallmark_gsea,c(1,2,6,8,10:12),color = c('#7B68EE','#CD3333','#20B2AA','#FF8C00','#FF6666','darkblue','green'),
          title = 'hallmark GSEA',
          base_size = 11,
          rel_heights = c(1.5, 0.3, 0.5))
dev.off()
