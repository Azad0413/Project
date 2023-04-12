rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/BJTC-320")
if (! dir.exists("./00_rawdata")){
  dir.create("./00_rawdata")
}
setwd("./00_rawdata")
library(GEOquery)
library(Biobase)
library(tidyr)
library(AnnoProbe)
library(tidyverse)
## 01-1 GSE138198用做训练集------
gset<-getGEO("GSE138198",destdir = '.',GSEMatrix = T,getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
a=gset[[1]]
gpl='GPL6244'
probe2symbol=idmap(gpl)
dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
colnames(probe2symbol)<-c('ID','symbol')
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
## 13个HT实验 3个normal对照
pd<-pData(a)
group<-data.frame(sample=pd$geo_accession,
               group=pd$title)
group<-group[-(14:33),]
group$group<-c(rep('HT',13),rep('control',3))
dat.final<-dat[,colnames(dat)%in%group$sample]
write.table(dat.final,file = 'dat.final.xls',quote = F,sep = '\t',row.names = T)
write.table(group,file = 'group.xls',quote = F,sep = '\t',row.names = F)
## 01-2 GSE29315验证集 -----
gset_va<-getGEO("GSE29315",destdir = '.',GSEMatrix = T,getGPL = F)
expr_va<-as.data.frame(exprs(gset_va[[1]]))
a_va=gset_va[[1]]
pd_va<-pData(a_va)
gpl2<-getGEO("GPL8300",destdir = '.')

gpl2<-Table(gpl2)    
colnames(gpl2)
probe2symobl2<-gpl2 %>%
  select('ID','Gene Symbol')%>%
  filter('Gene Symbol'!='')%>%
  separate('Gene Symbol',c('symbol','drop'),sep = '///')%>%
  select(-drop)
probe2symobl2=probe2symobl2[probe2symobl2$symbol!='',]

dat_va<-expr_va
dat_va$ID<-rownames(dat_va)
dat_va$ID<-as.character(dat_va$ID)
probe2symobl2$ID<-as.character(probe2symobl2$ID)
dat_va<-dat_va %>%
  inner_join(probe2symobl2,by='ID')%>% 
  select(-ID)%>%     ## 去除多余信息
  select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
## 6个实验8个对照
table(pd_va$`clinical description:ch1`)
group_va<-data.frame(sample=pd_va$geo_accession,group=pd_va$`clinical description:ch1`)
group_va<-group_va[c(40:53),]
dat_va<-dat_va[,colnames(dat_va)%in%group_va$sample]
table(group_va$group)
group_va$group<-ifelse(group_va$group=='hyperplasia','TPH','HT')
write.table(dat_va,'dat_va.xls',sep = '\t',quote = F,row.names = T)
write.table(group_va,'group_va.xls',sep = '\t',row.names = F,quote = F)
