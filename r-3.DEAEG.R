rm(list = ls())
# 02 差异分析----------
setwd("/data/nas1/luchunlin/project/BJTC-321")
if (! dir.exists("./03_DEAEG")){
  dir.create("./03_DEAEG")
}
setwd("./03_DEAEG")
library(lance)
library(tidyverse)
library(readxl)
diff<-read.delim2('/data/nas1/luchunlin/project/BJTC-321/01_DEGs/DEG_sig.xls',row.names = 1)
ecm<-read_xlsx('ECMgeneset.xlsx')
modgene<-read.delim2('/data/nas1/luchunlin/project/BJTC-321/02_WGCNA/modGene.xls')
DEmod<-modgene[modgene$modgene%in%rownames(diff),]%>%as.data.frame()
##1347

DEAEG<-DEmod[DEmod$.%in%ecm$geneset,]%>%as.data.frame()
## 62
write.table(DEAEG,file = 'DEAEG.xls',sep = '\t',row.names = F,quote = F)
library(ggvenn)
mydata<-list('DEGs'=rownames(diff),'ECM'=ecm$geneset,'MOD'=modgene$modgene)
pdf('DEAEG.venn.pdf',w=5,h=5)
ggvenn(mydata,c('DEGs','ECM','MOD'),
       fill_color = c("#ffb2b2","#b2e7cb","#F0E442"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb","#F0E442"),
       text_color = 'black')
dev.off()
png('DEAEG.venn.png',w=400,h=400)
ggvenn(mydata,c('DEGs','ECM','MOD'),
       fill_color = c("#ffb2b2","#b2e7cb","#F0E442"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb","#F0E442"),
       text_color = 'black')
dev.off()
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(GOplot)
gene_transform <- bitr(DEAEG$.,
                       fromType = "SYMBOL",
                       toType = c("ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
ego <- enrichGO(gene = gene_transform$ENTREZID,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                readable = TRUE)
write.table(ego,file = "GO.xls",sep = "\t",quote = F,row.names = F)
# 展示富集最显著的 GO term
go_bar<-barplot(ego, showCategory=5, split="ONTOLOGY",label_format = 100) +
  facet_grid(ONTOLOGY ~ ., scales = "free",space = 'free')
go_bar
ggsave('GO_bar.pdf',go_bar,height = 7,width = 7)
ggsave('GO_bar.png',go_bar,height = 7,width = 7)
## KEGG富集分析（气泡图）
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "hsa",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 2)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk@result,file = "KEGG.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk, showCategory=20,label_format = 50)
kk_dot
ggsave('KEGG_dot.pdf',kk_dot,height = 4,width = 7)
ggsave('KEGG_dot.png',kk_dot,height = 4,width = 7)
kk2<-data.frame(Category = "ALL",ID = kk$ID,Term = kk$Description, Genes = gsub("/", ", ", kk$geneID), adj_pval = kk$p.adjust)
sigind<-diff[rownames(diff)%in%DEAEG$.,]
genelist<-data.frame(ID=rownames(sigind),
                     logFC=sigind$logFC)
genelist$logFC<-as.numeric(genelist$logFC)
rownames(genelist)<-genelist$ID
genelist<-genelist[order(genelist$logFC,decreasing = T),]

## 读取logFC文件
circ<-circle_dat(kk2,genelist)
# GOBubble(circ, labels = 3,table.legend =F)
# GOCircle(circ,rad1=2.5,rad2=3.5,label.size=4,nsub=10) 
#termNum<-20
#geneNum<-nrow(genelist)
chord<-chord_dat(circ,genelist)
kegg_chord<-GOChord(chord,
                    gene.order = 'logFC',
                    gene.space = 0.25,
                    gene.size = 8,
                    space = 0.01,
                    lfc.col = c('red','white','blue'),
                    process.label = 13)
kegg_chord
