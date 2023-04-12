rm(list = ls())
# 01 获取数据集--------------
setwd("/data/nas1/luchunlin/project/BJTC-356/")
if (! dir.exists("./03_enrichment(DEGs)")){
  dir.create("./03_enrichment(DEGs)")
}
setwd("./03_enrichment(DEGs)")

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(GOplot)
deg <- read.delim2('../02_DEGs/DEG_sig.xls')

gene_transform <- bitr(rownames(deg),
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
go_bar<-dotplot(ego, showCategory=5, split="ONTOLOGY",label_format = 50) +
  facet_grid(ONTOLOGY ~ ., scales = "free",space = 'free')
go_bar
ggsave('01.GO_dot.png',go_bar,width =8,height = 6)
ggsave('01.GO_dot.pdf',go_bar,width =8,height = 6)
## GO 圈图
go_result<-read.table('GO.xls',header = T,sep = '\t',check.names = F)
go2=data.frame(Category='ALL',ID=go_result$ID,Term=go_result$Description,Genes=gsub("/", ", ", go_result$geneID), adj_pval = go_result$p.adjust)
sigind<-deg
genelist<-data.frame(ID=rownames(sigind),
                     logFC=sigind$logFC)
genelist$logFC <- as.numeric(genelist$logFC)
rownames(genelist)<-genelist$ID
genelist<-genelist[order(genelist$logFC,decreasing = T),]

## 读取logFC文件
circ<-circle_dat(go2,genelist)
circ <- circ[c(1:196),]
# GOBubble(circ, labels = 3,table.legend =F)
# GOCircle(circ,rad1=2.5,rad2=3.5,label.size=4,nsub=10) 
#termNum<-20
#geneNum<-nrow(genelist)
chord<-chord_dat(circ,genelist)
trace(GOChord,edit = T)

kegg_chord<-GOChord(chord,
                    gene.order = 'logFC',
                    gene.space = 0.3,
                    gene.size = 4,
                    space = 0.01,
                    lfc.col = c('red','white','blue'),
                    process.label = 10)

kegg_chord
ggsave(kegg_chord,filename = '02.GO.chord.pdf',w=12,h=13)
ggsave(kegg_chord,filename = '02.GO.chord.png',w=12,h=13)
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
ggsave('03.KEGG_dot.png',kk_dot,width =7,height = 5)
ggsave('03.KEGG_dot.pdf',kk_dot,width =7,height = 5)
