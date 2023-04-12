rm(list = ls())
setwd("/data/nas1/luchunlin/project/LZZK-512/")
if (! dir.exists("./02_DEGRG")){
  dir.create("./02_DEGRG")
}
setwd("./02_DEGRG")
sig.diff<-read.delim2('/data/nas1/luchunlin/project/LZZK-512/01_DEGs/DEG_sig.xls',row.names = 1)
gtgeneset<-read_xlsx('GT_Geneset.xlsx')

DEGRG<-gtgeneset[gtgeneset$Geneset%in%rownames(sig.diff),]%>%as.data.frame()
##33
write.table(DEGRG,file = 'DEGRG.xls',sep = '\t',quote = F,row.names = F)
library(ggvenn)
mydata<-list('DEGs'=rownames(sig.diff),'Glycosyl transferase'=gtgeneset$Geneset)
pdf("DEGRG.venn.pdf",w=5,h=5)
ggvenn(mydata,c('DEGs','Glycosyl transferase'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()
png("DEGRG.venn.png",w=500,h=500)
ggvenn(mydata,c('DEGs','Glycosyl transferase'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()
## GO/KEGG
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(GOplot)
gene_transform <- bitr(DEGRG$Geneset,
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
                 qvalueCutoff = 0.05)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk@result,file = "KEGG.xls",sep = "\t",quote = F,row.names = F)
kk2<-data.frame(Category = "ALL",ID = kk$ID,Term = kk$Description, Genes = gsub("/", ", ", kk$geneID), adj_pval = kk$p.adjust)
## 读取logFC文件
df<-read.delim2('/data/nas1/luchunlin/project/LZZK-512/01_DEGs/DEG_sig.xls',row.names = 1)%>%lc.tableToNum()
df<-df[DEGRG$Geneset,]
idfc2<-df
genelist<-data.frame(ID=rownames(idfc2),logFC=idfc2$log2FoldChange)
circ<-circle_dat(kk2,genelist)
#termNum<-20
#geneNum<-nrow(genelist)
chord<-chord_dat(circ,genelist)
trace(GOChord,edit = T)
kegg_chord<-GOChord(chord,
                    gene.order = 'logFC',
                    gene.space = 0.25,
                    gene.size = 4,
                    space = 0.01,
                    lfc.col = c('red','white','blue'),
                    process.label = 12)
kegg_chord
