rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/GY0324-12/")
if (! dir.exists("./04_DEGkegg")){
  dir.create("./04_DEGkegg")
}
setwd("./04_DEGkegg")

library(clusterProfiler)
library(org.Rn.eg.db)
library(ggplot2)
library(enrichplot)
library(GOplot)
## A1-----
sig.diff <- read.delim2('../03_DEGs/DEG_sig(A1).xls')
gene_transform <- bitr(rownames(sig.diff),
                       fromType = "SYMBOL",
                       toType = c("ENTREZID",'ENSEMBL'),
                       OrgDb = "org.Rn.eg.db")

## KEGG富集分析（气泡图）
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "rat",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 1)
kk <- setReadable(kk, OrgDb = org.Rn.eg.db, keyType="ENTREZID")
write.table(kk@result,file = "KEGG(A1vsC5).xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk, showCategory=15,label_format = 50)
kk_dot
ggsave('01.KEGG_dot(A1vsC5).png',kk_dot,width =7,height = 5)
ggsave('01.KEGG_dot(A1vsC5).pdf',kk_dot,width =7,height = 5)


## A2-----
sig.diff <- read.delim2('../03_DEGs/DEG_sig(A2).xls')
gene_transform <- bitr(rownames(sig.diff),
                       fromType = "SYMBOL",
                       toType = c("ENTREZID",'ENSEMBL'),
                       OrgDb = "org.Rn.eg.db")

## KEGG富集分析（气泡图）
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "rat",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 1)
kk <- setReadable(kk, OrgDb = org.Rn.eg.db, keyType="ENTREZID")
write.table(kk@result,file = "KEGG(A2vsC5).xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk, showCategory=15,label_format = 50)
kk_dot
ggsave('02.KEGG_dot(A2vsC5).png',kk_dot,width =7,height = 5)
ggsave('02.KEGG_dot(A2vsC5).pdf',kk_dot,width =7,height = 5)

## A3-----
sig.diff <- read.delim2('../03_DEGs/DEG_sig(A3).xls')
gene_transform <- bitr(rownames(sig.diff),
                       fromType = "SYMBOL",
                       toType = c("ENTREZID",'ENSEMBL'),
                       OrgDb = "org.Rn.eg.db")

## KEGG富集分析（气泡图）
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "rat",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 1)
kk <- setReadable(kk, OrgDb = org.Rn.eg.db, keyType="ENTREZID")
write.table(kk@result,file = "KEGG(A3vsC5).xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk, showCategory=15,label_format = 50)
kk_dot
ggsave('03.KEGG_dot(A3vsC5).png',kk_dot,width =7,height = 5)
ggsave('03.KEGG_dot(A3vsC5).pdf',kk_dot,width =7,height = 5)

## B4-----
sig.diff <- read.delim2('../03_DEGs/DEG_sig(B4).xls')
gene_transform <- bitr(rownames(sig.diff),
                       fromType = "SYMBOL",
                       toType = c("ENTREZID",'ENSEMBL'),
                       OrgDb = "org.Rn.eg.db")

## KEGG富集分析（气泡图）
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "rat",
                 pAdjustMethod = "none",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 1)
kk <- setReadable(kk, OrgDb = org.Rn.eg.db, keyType="ENTREZID")
write.table(kk@result,file = "KEGG(B4vsC5).xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk, showCategory=15,label_format = 50)
kk_dot
ggsave('04.KEGG_dot(B4vsC5).png',kk_dot,width =7,height = 5)
ggsave('04.KEGG_dot(B4vsC5).pdf',kk_dot,width =7,height = 5)


