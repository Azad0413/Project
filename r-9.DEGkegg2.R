rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/GY0324-12/")
if (! dir.exists("./09_DEGkegg2")){
  dir.create("./09_DEGkegg2")
}
setwd("./09_DEGkegg2")

library(clusterProfiler)
library(org.Rn.eg.db)
library(ggplot2)
library(enrichplot)
library(GOplot)
## A1-----
sig.diff1 <- read.delim2('../08_DEGs2/DEG_sig(A1vsA2).xls')
gene_transform <- bitr(rownames(sig.diff1),
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
write.table(kk@result,file = "KEGG(A1vsA2).xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk, showCategory=15,label_format = 50)
kk_dot
ggsave('01.KEGG_dot(A1vsA2).png',kk_dot,width =8,height = 5)
ggsave('01.KEGG_dot(A1vsA2).pdf',kk_dot,width =8,height = 5)


## A2-----
sig.diff2 <- read.delim2('../08_DEGs2/DEG_sig(A2vs.A3).xls')
gene_transform <- bitr(rownames(sig.diff2),
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
write.table(kk@result,file = "KEGG(A2vsA3).xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk, showCategory=15,label_format = 50)
kk_dot
ggsave('02.KEGG_dot(A2vsA3).png',kk_dot,width =7,height = 5)
ggsave('02.KEGG_dot(A2vsA3).pdf',kk_dot,width =7,height = 5)



###交集
interDEG <- data.frame(symbol=intersect(rownames(sig.diff1),rownames(sig.diff2)))
##1119
write.table(interDEG,file = 'interDEG.xls',sep = '\t',row.names = F,quote = F)

library(ggvenn)
mydata<-list('DEGs(A1vs.A2)'=rownames(sig.diff1),'DEGs(A2vs.A3)'=rownames(sig.diff2))
pdf('03.venn.pdf',w=5,h=5)
ggvenn(mydata,c('DEGs(A1vs.A2)','DEGs(A2vs.A3)'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       text_size = 4,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()
png('03.venn.png',w=400,h=400)
ggvenn(mydata,c('DEGs(A1vs.A2)','DEGs(A2vs.A3)'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       text_size = 4,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()






gene_transform <- bitr(interDEG$symbol,
                       fromType = "SYMBOL",
                       toType = c("ENTREZID",'ENSEMBL'),
                       OrgDb = "org.Rn.eg.db")

##富集
ego <- enrichGO(gene = gene_transform$ENTREZID,
                OrgDb = org.Rn.eg.db,
                keyType = "ENTREZID",
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
write.table(ego,file = "GO.xls",sep = "\t",quote = F,row.names = F)
# 展示富集最显著的 GO term
go_bar<-barplot(ego, showCategory=5, split="ONTOLOGY",label_format = 70) +
  facet_grid(ONTOLOGY ~ ., scales = "free",space = 'free')
go_bar
ggsave('04.GO_bar.png',go_bar,width =8,height = 6)
ggsave('04.GO_bar.pdf',go_bar,width =8,height = 6)

## KEGG富集分析（气泡图）
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "rat",
                 pAdjustMethod = "none",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 1)
kk <- setReadable(kk, OrgDb = org.Rn.eg.db, keyType="ENTREZID")
write.table(kk@result,file = "KEGG(interDEG).xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk, showCategory=15,label_format = 50)
kk_dot
ggsave('05.KEGG_dot(interDEG).png',kk_dot,width =7,height = 5)
ggsave('05.KEGG_dot(interDEG).pdf',kk_dot,width =7,height = 5)
