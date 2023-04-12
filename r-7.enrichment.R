rm(list = ls())
# 01 获取数据集--------------
setwd("/data/nas1/luchunlin/project/BJTC-356/")
if (! dir.exists("./07_enrichment")){
  dir.create("./07_enrichment")
}
setwd("./07_enrichment")
hubgene <- read.delim2('../06_quadrant/final.gene.xls')
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(GOplot)

gene_transform <- bitr(hubgene$symbol,
                       fromType = "SYMBOL",
                       toType = c("ENTREZID",'ENSEMBL'),
                       OrgDb = "org.Hs.eg.db")
# ego <- enrichGO(gene = gene_transform$ENTREZID,
#                 OrgDb = org.Hs.eg.db,
#                 keyType = "ENTREZID",
#                 ont = "ALL",
#                 pAdjustMethod = "BH",
#                 pvalueCutoff = 0.05,
#                 qvalueCutoff = 0.05,
#                 readable = TRUE)
# write.table(ego,file = "GO.xls",sep = "\t",quote = F,row.names = F)
# # 展示富集最显著的 GO term
# go_bar<-barplot(ego, showCategory=5, split="ONTOLOGY",label_format = 70) +
#   facet_grid(ONTOLOGY ~ ., scales = "free",space = 'free')
# go_bar
# ggsave('01.GO_bar.png',go_bar,width =8,height = 6)
# ggsave('01.GO_bar.pdf',go_bar,width =8,height = 6)
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

ggsave('02.KEGG_dot.png',kk_dot,width =8,height = 6)
ggsave('02.KEGG_dot.pdf',kk_dot,width =8,height = 6)

kegg.meta <- read_xlsx('kegg_meta.xlsx')
kk.result <- kk@result%>%as.data.frame()
#kk.result <- kk.result[which(kk.result$pvalue<0.05),]

intersect <- kegg.meta[kegg.meta$`#Pathway`%in% kk.result$Description,]
intersect <- intersect[,c(2,4)]
kk.result2 <- kk.result[kk.result$Description%in%intersect$`#Pathway`,]
kk.result2 <- kk.result2[,c(1,2,8)]
colnames(kk.result2)
colnames(intersect) <- c('Description','metabolin')

intersect <- merge(intersect,kk.result2,by='Description')
intersect <- intersect[,c(3,1,2,4)]
length(unique(intersect$ID))

write.table(intersect,file = 'interectkegg.xls',quote = F,sep = '\t',row.names = F)

