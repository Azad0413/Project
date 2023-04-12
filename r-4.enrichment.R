rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-370-8/")
if (! dir.exists("./04_enrichment")){
  dir.create("./04_enrichment")
}
setwd("./04_enrichment")

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(GOplot)
hubgene<-read.delim2('../03_DEGs/intersect.xls')
gene_transform <- bitr(hubgene$.,
                       fromType = "SYMBOL",
                       toType = c("ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
ego <- enrichGO(gene = gene_transform$ENTREZID,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "ALL",
                pAdjustMethod = "none",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
write.table(ego,file = "GO.xls",sep = "\t",quote = F,row.names = F)
# 展示富集最显著的 GO term
go_bar<-barplot(ego, showCategory=5, split="ONTOLOGY",label_format=50) +
  facet_grid(ONTOLOGY ~ ., scales = "free")
go_bar
ggsave(filename = '01.GO_bar.pdf',go_bar,w=7.5,h=6)
ggsave(filename = '01.GO_bar.png',go_bar,w=7.5,h=6)
## KEGG富集分析（气泡图）
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "hsa",
                 pAdjustMethod = "none",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 1)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk@result,file = "KEGG.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk, showCategory=15,label_format=50)
kk_dot
ggsave(filename = '02.KEGG_dot.pdf',kk_dot,w=6,h=4)
ggsave(filename = '02.KEGG_dot.png',kk_dot,w=6,h=4)
