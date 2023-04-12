## 02 差异基因鉴定-----
rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-308/")
if (! dir.exists("./02_DENRG")){
  dir.create("./02_DENRG")
}
setwd("./02_DENRG")
sig.diff<-read.delim2('/data/nas1/luchunlin/project/BJTC-308/01_DEGs/DEG_sig.xls',row.names = 1)
#necroptosis<-read_xlsx(file = 'Necroptosis.xlsx')
#inflam<-read_xlsx('inflammatory.xlsx')
lipid<-read_xlsx('LIPIDS.xlsx')
lipid<-lipid[!duplicated(lipid$Geneset),]
#invasion<-read_xlsx('invasion.xlsx')
#DENRG<-necroptosis[necroptosis$Symbol%in%rownames(sig.diff),]
#DENRG<-inflam[inflam$MAPPED_SYMBOLS%in%rownames(sig.diff),]
DENRG<-lipid[lipid$Geneset%in%rownames(sig.diff),]

library(ggvenn)
mydata<-list('DEGs'=rownames(sig.diff),'necroptosis'=necroptosis$Symbol)
ggvenn(mydata,c('DEGs','necroptosis'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
write.table(DENRG,file = 'DENRG.xls',sep = '\t',quote = F,row.names = F)

## GO/KEGG
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(GOplot)
gene_transform <- bitr(DENRG$Symbol,
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
kk_dot <- dotplot(kk, showCategory=15)
kk_dot