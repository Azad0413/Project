rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/GY0324-12/")
if (! dir.exists("./14_DEGkegg3")){
  dir.create("./14_DEGkegg3")
}
setwd("./14_DEGkegg3")

library(clusterProfiler)
library(org.Rn.eg.db)
library(ggplot2)
library(enrichplot)
library(GOplot)
##
sig.diff <- read.delim2('../13_DEGs3/DEG_sig(A3vsB4).xls')
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
write.table(kk@result,file = "KEGG(A3vsB4).xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk, showCategory=15,label_format = 50)
kk_dot
ggsave('01.KEGG_dot(A3vsB4).png',kk_dot,width =7,height = 5)
ggsave('01.KEGG_dot(A3vsB4).pdf',kk_dot,width =7,height = 5)

