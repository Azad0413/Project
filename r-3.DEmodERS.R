rm(list = ls())
setwd("/data/nas1/luchunlin/project/GY0315/")
if (! dir.exists("./03_DEmodERS")){
  dir.create("./03_DEmodERS")
}
setwd("./03_DEmodERS")
modgene<-read.delim2('../02_WGCNA/modGene.xls')
diff<-read.delim2('../01_DEGs/DEG_sig.xls',row.names = 1)
DEmod<-modgene[modgene$modgene%in%rownames(diff),]%>%as.data.frame()
ERS<-read.csv('GeneCards-SearchResults.csv')
ERS<-ERS[which(ERS$Relevance.score>5),]
intersect<-DEmod[DEmod$.%in%ERS$Gene.Symbol,]%>%as.data.frame()
write.table(intersect,file = 'intersect.xls',sep = '\t',row.names = F,quote = F)
library(ggvenn)
mydata<-list(modGene=modgene$modgene,ERS=ERS$Gene.Symbol,DEGs=rownames(diff))
pdf(file = '01.hubgene.pdf',w=6,h=6)
ggvenn(mydata,c('modGene','ERS','DEGs'),
       fill_color = c("#ffb2b2","#b2e7cb","#F0E442"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb","#F0E442"),
       text_color = 'black')
dev.off()
png(file = '01.hubgene.png',w=500,h=500)
ggvenn(mydata,c('modGene','ERS','DEGs'),
       fill_color = c("#ffb2b2","#b2e7cb","#F0E442"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb","#F0E442"),
       text_color = 'black')
dev.off()
# 04 GO/KEGG----------

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(GOplot)
hubgene<-read.delim2('../03_DEmodERS/intersect.xls')

gene_transform <- bitr(hubgene$.,
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
go_bar<-barplot(ego, showCategory=5, split="ONTOLOGY",label_format=50) +
  facet_grid(ONTOLOGY ~ ., scales = "free",space = 'free')
go_bar
ggsave(filename = '02.GO_bar.pdf',go_bar,w=7.5,h=6)
ggsave(filename = '02.GO_bar.png',go_bar,w=7.5,h=6)
## KEGG富集分析（气泡图）
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "hsa",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk@result,file = "KEGG.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk, showCategory=15,label_format=50)
kk_dot
ggsave(filename = '03.KEGG_dot.pdf',kk_dot,w=6.5,h=4)
ggsave(filename = '03.KEGG_dot.png',kk_dot,w=6.5,h=4)
