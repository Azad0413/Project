rm(list = ls())
setwd("/data/nas1/luchunlin/project/SJZZK-431-10/")
if (! dir.exists("./13_enrichment(risk)")){
  dir.create("./13_enrichment(risk)")
}
setwd("./13_enrichment(risk)")
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(GOplot)
DEG <- read.delim2('../12_DEGs(risk)/DEG_sig(risk).xls')
gene_transform <- bitr(rownames(DEG),
                       fromType = "SYMBOL",
                       toType = c("ENTREZID",'ENSEMBL'),
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
go_bar<-barplot(ego, showCategory=5, split="ONTOLOGY",label_format = 50) +
  facet_grid(ONTOLOGY ~ ., scales = "free",space = 'free')
go_bar
ggsave('GO_bar.png',go_bar,width =8,height = 6)
ggsave('GO_bar.pdf',go_bar,width =8,height = 6)
## KEGG富集分析（气泡图）
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "hsa",
                 pAdjustMethod = "none",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 1)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk@result,file = "KEGG.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk, showCategory=15,label_format = 50)
kk_dot
ggsave('KEGG_dot.png',kk_dot,width =7,height = 5)
ggsave('KEGG_dot.pdf',kk_dot,width =7,height = 5)
