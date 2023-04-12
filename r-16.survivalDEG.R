rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-302")
if (! dir.exists("./16_survival")){
  dir.create("./16_survival")
}
setwd("./16_survival")

### GSE19743 存活和死亡差异
library(GEOquery)
library(Biobase)
library(tidyverse)
library(dplyr)
gset<-getGEO("GSE19743",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
a=gset[[1]]
pd<-pData(a)
gpl<-getGEO("GPL570",destdir = '.')
gpl<-Table(gpl)    
colnames(gpl)
probe2symbol<-gpl %>%
  dplyr::select('ID','Gene Symbol')%>%
  filter('Gene Symbol'!='')%>%
  separate('Gene Symbol',c('symbol','drop'),sep = '///')%>%
  dplyr::select(-drop)
#probe2symbol=probe2symbol[probe2symbol$symbol!='',]
dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symbol$ID<-as.character(probe2symbol$ID)
dat<-dat %>%
  inner_join(probe2symbol,by='ID')%>% 
  dplyr::select(-ID)%>%     ## 去除多余信息
  dplyr::select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
group<-data.frame(sample=pd$geo_accession,group=pd$characteristics_ch1.8)
group<-group[c(1:114),]
table(group$group)
group$group<-ifelse(group$group=='survival: No','Dead','Alive')
group<-group[order(group$group),]
dat.final<-dat[,group$sample]
write.table(dat.final,file = 'dat_survival.xls',sep = '\t',row.names = T,quote = F)
write.table(group,file = 'gorup_survival.xls',sep = '\t',row.names = F,quote = F)

hubgene<-read.delim2('/data/nas1/luchunlin/project/BJTC-302/10_validation/hubgene.final.xls')
hub.exp<-dat.final[hubgene$hubgene,]
## 热图-------
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(gplots)
group_rt<-group
group_rt<-group_rt$group%>%as.data.frame()
colnames(group_rt)<-'group'
rownames(group_rt)<-colnames(hub.exp)
heat<-hub.exp
x<-log2(heat+1)
#x<-t(scale(t(heat)))
ann_colors<-list(
  group = c(Alive="#00CED1",Dead="#F08080"))

pheatmap(mat=x,
         annotation_col = group_rt,
         color=bluered(100),
         scale = "row",
         annotation_colors = ann_colors,
         fontsize = 10,
         show_colnames = F,
         cluster_cols = F,
         cluster_rows = T,
         show_rownames = T,
         annotation_names_row = F,
         cellwidth = 2.2)

