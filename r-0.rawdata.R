rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/LZZK-503")
if (! dir.exists("./00_rawdata")){
  dir.create("./00_rawdata")
}
setwd("./00_rawdata")
library(GEOquery)
library(Biobase)
gset<-getGEO("GSE119054",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))

##GSE119054
library(AnnoProbe)
gpl<-'GPL19615'
probe2symbol<-idmap(gpl,type = 'pipe')

a=gset[[1]]
gpl<-Table(gpl)    
colnames(probe2symbol)<-c('ID','symbol')
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
pd<-pData(a)

write.table(dat,file = 'dat.xls',sep = '\t',quote = F,row.names = T)

group<-data.frame(sample=pd$geo_accession,group=c(rep('OV',6),rep('control',3)))
write.table(group,file = 'group.xls',sep = '\t',row.names = F,quote = F)
