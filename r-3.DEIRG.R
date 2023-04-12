## 02 差异基因鉴定-----
rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-302/")
if (! dir.exists("./03_DEIRG")){
  dir.create("./03_DEIRG")
}
setwd("./03_DEIRG")
library(readxl)
## 与免疫基因集取交集
immport<-read_xlsx('Immport_IRGs.xlsx')
innatedb<-read_xls('innatedb.xls')
immport<-data.frame(symbol=immport$Symbol)
innatedb<-data.frame(symbol=innatedb$`Gene Symbol`)
immune<-rbind(immport,innatedb)
immune<-immune[!duplicated(immune$symbol),]%>%as.data.frame()
##2533
modgene<-read.delim2('/data/nas1/luchunlin/project/BJTC-302/02_WGCNA/DEmod.xls')

DEIRG<-immune[immune$.%in%modgene$.,]%>%as.data.frame()
##113
write.table(DEIRG,file = 'DEIRG.xls',sep = '\t',row.names = F,quote = F)
library(ggvenn)
mydata<-list('DE-mod'=modgene$.,'Immune'=immune$.)
ggvenn(mydata,c('DE-mod','Immune'),
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

gene_transform <- bitr(DEIRG$.,
                       fromType = "SYMBOL",
                       toType = c("ENTREZID",'ENSEMBL'),
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


