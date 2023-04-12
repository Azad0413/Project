rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/BJTC-317")
if (! dir.exists("./00_rawdata")){
  dir.create("./00_rawdata")
}
setwd("./00_rawdata")
### GEO数据库 GSE7621
library(GEOquery)
library(Biobase)
library(tidyr)
library(AnnoProbe)
library(tidyverse)
## 01-1 GSE26927集------
gset<-getGEO("GSE26927",destdir = '.',GSEMatrix = T,getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
a=gset[[1]]
gpl<-getGEO("GPL6255",destdir = '.')
gpl<-Table(gpl)    
colnames(gpl)
probe2symbol<-gpl %>%
  select('ID','SYMBOL')%>%
  filter('SYMBOL'!='')%>%
  separate('SYMBOL',c('symbol','drop'),sep = '///')%>%
  select(-drop)
dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symbol$ID<-as.character(probe2symbol$ID)
dat<-dat %>%
  inner_join(probe2symbol,by='ID')%>% 
  select(-ID)%>%     ## 去除多余信息
  select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
## 16个PD实验 9个control对照
pd<-pData(a)
group<-data.frame(sample=pd$geo_accession,
                  group=pd$title)
table(group$group)
group<-group[c(79:98),]

group$group<-c(rep('control',8),rep('PD',12))
dat.final<-dat[,colnames(dat)%in%group$sample]
#dat.final<-log2(dat.final+1)
write.table(dat.final,file = 'dat.final.xls',quote = F,sep = '\t',row.names = T)
write.table(group,file = 'group.xls',quote = F,sep = '\t',row.names = F)

