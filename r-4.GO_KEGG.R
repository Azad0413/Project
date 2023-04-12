rm(list = ls())
# 04 GO/KEGG----------
setwd("/data/nas1/luchunlin/project/BJTC-317")
if (! dir.exists("./04_GO_KEGG")){
  dir.create("./04_GO_KEGG")
}
setwd("./04_GO_KEGG")

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(GOplot)
diff<-read.delim2('/data/nas1/luchunlin/project/BJTC-317/01_DEGs/DEG_sig.xls')
DEMID<-read.delim2('/data/nas1/luchunlin/project/BJTC-317/03_DEMDG/DEMDG.xls')
gene_transform <- bitr(DEMID$.,
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
go_dot<-dotplot(ego, showCategory=5, split="ONTOLOGY",color='p.adjust') +
  facet_grid(ONTOLOGY ~ ., scales = "free")
go_dot

ggsave('GO.dot.pdf',height = 10,width = 7)
ggsave('GO.dot.png',height = 10,width = 7)
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
ggsave('KEGG_dot.pdf',height = 6,width = 6)
ggsave('KEGG_dot.png',height = 6,width = 6)

kk2<-data.frame(Category = "ALL",ID = kk$ID,Term = kk$Description, Genes = gsub("/", ", ", kk$geneID), adj_pval = kk$p.adjust)
## 读取logFC文件
sigind<-diff[rownames(diff)%in%DEMID$.,]
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
