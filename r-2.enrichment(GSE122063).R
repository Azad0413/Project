rm(list = ls())
setwd("/data/nas1/luchunlin/project/LZZK-519-10/")
if (! dir.exists("./02_enrichment(GSE122063)")){
  dir.create("./02_enrichment(GSE122063)")
}
setwd("./02_enrichment(GSE122063)")
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(GOplot)
sig.diff1 <- read.delim2('../01_DEGs(GSE122063)/DEG_sig(GSE122063).xls')
gene_transform <- bitr(rownames(sig.diff1),
                       fromType = "SYMBOL",
                       toType = c("ENTREZID",'ENSEMBL'),
                       OrgDb = "org.Hs.eg.db")
ego <- enrichGO(gene = gene_transform$ENTREZID,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "ALL",
                pAdjustMethod = "none",
                pvalueCutoff = 0.05,
                qvalueCutoff = 1,
                readable = TRUE)
bp <- enrichGO(gene = gene_transform$ENTREZID,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "BP",
                pAdjustMethod = "none",
                pvalueCutoff = 0.05,
                qvalueCutoff = 1,
                readable = TRUE)
cc <- enrichGO(gene = gene_transform$ENTREZID,
               OrgDb = org.Hs.eg.db,
               keyType = "ENTREZID",
               ont = "CC",
               pAdjustMethod = "none",
               pvalueCutoff = 0.05,
               qvalueCutoff = 1,
               readable = TRUE)
mf <- enrichGO(gene = gene_transform$ENTREZID,
               OrgDb = org.Hs.eg.db,
               keyType = "ENTREZID",
               ont = "MF",
               pAdjustMethod = "none",
               pvalueCutoff = 0.05,
               qvalueCutoff = 1,
               readable = TRUE)

write.table(ego,file = "GO.xls",sep = "\t",quote = F,row.names = F)
# 展示富集最显著的 GO term
# go_bar<-barplot(ego, showCategory=5, split="ONTOLOGY",label_format = 70) +
#   facet_grid(ONTOLOGY ~ ., scales = "free",space = 'free')
# go_bar

bp_dot <- dotplot(bp, showCategory=10,label_format = 50)+
  ggtitle(label = 'Biological Process')
bp_dot
ggsave('01.bp_dot.png',bp_dot,width =8,height = 6)
ggsave('01.bp_dot.pdf',bp_dot,width =8,height = 6)
cc_dot <- dotplot(cc, showCategory=10,label_format = 50)+
  ggtitle(label = 'Cellular Component')
cc_dot
ggsave('02.cc_dot.png',cc_dot,width =8,height = 6)
ggsave('02.cc_dot.pdf',cc_dot,width =8,height = 6)
mf_dot <- dotplot(mf, showCategory=10,label_format = 50)+
  ggtitle(label = 'Molecular Function')
mf_dot
ggsave('03.mf_dot.png',mf_dot,width =8,height = 6)
ggsave('03.mf_dot.pdf',mf_dot,width =8,height = 6)

## KEGG富集分析（气泡图）
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "hsa",
                 pAdjustMethod = "none",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 1)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk@result,file = "KEGG.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- barplot(kk, showCategory=10,label_format = 50)
kk_dot
ggsave('04.KEGG_bar.png',kk_dot,width =8,height = 6)
ggsave('04.KEGG_bar.pdf',kk_dot,width =8,height = 6)
