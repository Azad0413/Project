rm(list = ls())
# 03 GO/KEGG----------
setwd("/data/nas1/luchunlin/project/JNZK-207")
if (! dir.exists("./02_GO_KEGG")){
  dir.create("./02_GO_KEGG")
}
setwd("./02_GO_KEGG")
age<-read.csv(file = 'genage_human.csv')
sig_diff<-read.delim2("/data/nas1/luchunlin/project/JNZK-207/01_DEGs/DEG_sig.xls", row.names = 1)
AGDEs<-sig_diff[rownames(sig_diff)%in%age$symbol,]
write.table(AGDEs,file = 'AGDEs.xls',sep = '\t',quote = F,row.names = T)
library(ggvenn)
mydata<-list('AG'=age$symbol,'DEGs'=rownames(sig_diff))
ggvenn(mydata,c('AG','DEGs'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(GOplot)
gene_transform <- bitr(rownames(AGDEs),
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
kk_dot <- dotplot(kk, showCategory=20)
kk_dot

