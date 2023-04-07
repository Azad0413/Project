rm(list = ls())
setwd("/data/nas1/luchunlin/project/HZ0301-3/")
if (! dir.exists("./04_enrichment")){
  dir.create("./04_enrichment")
}
setwd("./04_enrichment")
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(GOplot)
intersect <- read.delim2('../03_intersect/intersect.xls')
gene_transform <- bitr(intersect$symbol,
                       fromType = "SYMBOL",
                       toType = c("ENTREZID",'ENSEMBL'),
                       OrgDb = "org.Hs.eg.db")
ego <- enrichGO(gene = gene_transform$ENTREZID,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 1,
                readable = TRUE)
write.table(ego,file = "GO.xls",sep = "\t",quote = F,row.names = F)
# 展示富集最显著的 GO term
go_bar<-dotplot(ego, showCategory=5, split="ONTOLOGY",label_format = 70) +
  facet_grid(ONTOLOGY ~ ., scales = "free",space = 'free')
go_bar
ggsave('01.GO_dot.png',go_bar,width =10,height = 8)
ggsave('01.GO_dot.pdf',go_bar,width =10,height = 8)
## KEGG富集分析（气泡图）
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "hsa",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 1)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk@result,file = "KEGG.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk, showCategory=15,label_format = 50)
kk_dot
# ggsave('02.KEGG_dot.png',kk_dot,width =5,height = 4)
# ggsave('02.KEGG_dot.pdf',kk_dot,width =5,height = 4)
