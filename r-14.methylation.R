rm(list = ls())
#10 诊断基因功能相似性分析（GOSemSim包）----------
setwd("/data/nas1/luchunlin/project/BJTC-321")
if (! dir.exists("./14_methylation")){
  dir.create("./14_methylation")
}
setwd("./14_methylation")
library(GEOquery)
library(Biobase)
library(tidyr)
library(tidyverse)
library(readxl)
## 01-1 GSE67472------
gset<-getGEO("GSE104472",destdir = '.',GSEMatrix = T,getGPL = F)
a=gset[[1]]
pd<-pData(a)
expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL13534",destdir = '.')
gpl<-Table(gpl) 
colnames(gpl)
probe2symbol<-gpl %>%
  dplyr::select('ID','UCSC_RefGene_Name')%>%
  filter('UCSC_RefGene_Name'!='')%>%
  separate('UCSC_RefGene_Name',c('symbol','drop'),sep = ';')%>%
  dplyr::select(-drop)
probe<-rownames(expr)%>%as.data.frame()

### hub基因的甲基化表达矩阵提取出来
hubgene <- read.delim2('/data/nas1/luchunlin/project/BJTC-321/08_va/hub_final.xls')
hub.probe<-probe2symbol[probe2symbol$symbol%in%hubgene$hubgene,]
table(hub.probe$symbol)
hub.me<-expr[rownames(expr)%in%hub.probe$ID,]
hub.probe<-hub.probe[hub.probe$ID%in%rownames(hub.me),]
hub.me<-hub.me[hub.probe$ID,]
group<-data.frame(sample=pd$geo_accession,group=pd$`disease state:ch1`)
table(group$group)
group<-group[order(group$group),]
hub.me<-hub.me[,group$sample]
## 热图-------
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(gplots)

group_rt<-group$group%>%as.data.frame()
colnames(group_rt)<-'group'
rownames(group_rt)<-group$sample
heat<-hub.me

ann_colors<-list(
  group = c(Normal="lightblue",Asthma="darkorange"),
  Change=c(CCL26="#FF0000",CST1="#436EEE",CST2='#FFCC66',CST4='#33FFCC',POSTN='#F0BBFF',SERPINB2='#7B68EE')
)
table(hub.probe$symbol)
annotation_raw<-data.frame(row.names = rownames(heat),
                           Change=factor(rep(c('CCL26','CST1','CST2','CST4','POSTN','SERPINB2'),c(6,5,5,6,6,5))))
pheatmap(mat=heat,
         annotation_col = group_rt,
         color=bluered(100),
         annotation_row = annotation_raw,
         scale = "row",
         annotation_colors = ann_colors,
         fontsize = 10,
         show_colnames = F,
         cluster_cols = F,
         cluster_rows = F,
         annotation_names_row = F)
