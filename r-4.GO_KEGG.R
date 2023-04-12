rm(list = ls())
# 04 GO/KEGG----------
setwd("/data/nas1/luchunlin/project/BJTC-320")
if (! dir.exists("./04_GO_KEGG")){
  dir.create("./04_GO_KEGG")
}
setwd("./04_GO_KEGG")

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(GOplot)
DEOSG<-read.delim2('/data/nas1/luchunlin/project/BJTC-320/03_DEIMG/DEOSG.xls')
gene_transform <- bitr(DEOSG$.,
                       fromType = "SYMBOL",
                       toType = c("ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
ego <- enrichGO(gene = gene_transform$ENTREZID,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
write.table(ego,file = "GO.xls",sep = "\t",quote = F,row.names = F)
# 展示富集最显著的 GO term
go_bar<-barplot(ego, showCategory=5, split="ONTOLOGY") +
  facet_grid(ONTOLOGY ~ ., scales = "free")
go_bar

## KEGG富集分析（气泡图）
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "hsa",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 2)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk@result,file = "KEGG.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk, showCategory=15)
kk_dot


