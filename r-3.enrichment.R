rm(list = ls())
setwd("/data/nas1/luchunlin/project/JNZK-255-2/")
if (! dir.exists("./03_enrichment")){
  dir.create("./03_enrichment")
}
setwd("./03_enrichment")
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(GOplot)
intersect <- read.delim2('../02_intersect/intersect.xls')
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
                qvalueCutoff = 0.2,
                readable = TRUE)
write.table(ego,file = "GO.xls",sep = "\t",quote = F,row.names = F)
# 展示富集最显著的 GO term
go_bar<-dotplot(ego, showCategory=10, split="ONTOLOGY",label_format = 70) +
  facet_grid(ONTOLOGY ~ ., scales = "free",space = 'free')
go_bar
ggsave('01.GO_dot.png',go_bar,width =10,height = 8)
ggsave('01.GO_dot.pdf',go_bar,width =10,height = 8)

## GO 圈图
go_result<-read.table('GO.xls',header = T,sep = '\t',check.names = F)
go2=data.frame(Category='ALL',ID=go_result$ID,Term=go_result$Description,Genes=gsub("/", ", ", go_result$geneID), adj_pval = go_result$p.adjust)
df <- read.delim2('../01_DEGs/DEG_sig.xls',row.names = 1)%>%lc.tableToNum()
idfc2<-df
genelist<-data.frame(ID=rownames(idfc2),logFC=idfc2$log2FoldChange)
rownames(genelist)=genelist[,1]
circ<-circle_dat(go2,genelist)
go_cir<-GOCircle(circ,rad1 = 2.5,rad2 = 3.5,label.size = 4,nsub = 10,table.legend = F)
go_cir
ggsave('02.GO_cir.png',go_cir,width =7,height = 7)
ggsave('02.GO_cir.pdf',go_cir,width =6,height = 6)
## KEGG富集分析（气泡图）
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "hsa",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.2)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk@result,file = "KEGG.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- barplot(kk, showCategory=15,label_format = 50)
kk_dot
ggsave('04.KEGG_dot.png',kk_dot,width =7,height = 4)
ggsave('04.KEGG_dot.pdf',kk_dot,width =7,height = 4)
kk2<-data.frame(Category = "ALL",ID = kk$ID,Term = kk$Description, Genes = gsub("/", ", ", kk$geneID), adj_pval = kk$p.adjust)
## 读取logFC文件
circ<-circle_dat(kk2,genelist)
# GOBubble(circ, labels = 3,table.legend =F)
# GOCircle(circ,rad1=2.5,rad2=3.5,label.size=4,nsub=10) 
#termNum<-20
#geneNum<-nrow(genelist)
chord<-chord_dat(circ,genelist)
trace(GOChord,edit = T)
kegg_chord<-GOChord(chord,
                    gene.order = 'logFC',
                    gene.space = 0.25,
                    gene.size = 5,
                    space = 0.01,
                    lfc.col = c('red','white','blue'),
                    process.label = 10)
kegg_chord
