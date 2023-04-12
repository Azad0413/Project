rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/JNZK-207")
if (! dir.exists("./00_rawdata")){
  dir.create("./00_rawdata")
}
setwd("./00_rawdata")
library(GEOquery)
library(Biobase)
library(tidyr)
library(AnnoProbe)
library(tidyverse)
## 01-1 GSE47552------
gset<-getGEO("GSE47552",destdir = '.',GSEMatrix = T,getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
a=gset[[1]]
gpl<-getGEO("GPL6244",destdir = '.')
gpl<-Table(gpl)    
colnames(gpl)
gpl<-"GPL6244"
probe2symbol=idmap(gpl)
colnames(probe2symbol)<-c('ID','symbol')
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
## MM样本和正常样本
pd<-pData(a)

group<-data.frame(sample=pd$geo_accession,
                  group=pd$title)

group<-group[c(21:66),]
group<-group[order(group$group,decreasing = T),]
group$group<-c(rep('control',5),rep('MM',41))
dat.final<-dat[,group$sample]
pd <- pd[pd$geo_accession%in%group$sample,]
write.table(pd,file = 'clinical(GSE6244).xls',sep = '\t',row.names = F,quote = F)

write.table(dat.final,file = 'dat.final.xls',quote = F,sep = '\t',row.names = T)
write.table(group,file = 'group.xls',quote = F,sep = '\t',row.names = F)

## 01-2 MMRF-------
expr2<-read_tsv(file = 'MMRF-COMMPASS.htseq_fpkm.tsv')
expr2<-as.data.frame(expr2)
rownames(expr2)<-expr2[,1]
expr2<-expr2[,-1]
genecode<-read.table(file = 'gencode.v22.annotation.gene.probeMap')
probe2symbol2<-genecode[,(1:2)]
probe2symbol2<-probe2symbol2[-1,]
colnames(probe2symbol2)<-c('ID','symbol')
dat2<-expr2
#dat2<-as.matrix(dat2)
dat2$ID<-rownames(dat2)
dat2$ID<-as.character(dat2$ID)
probe2symbol2$ID<-as.character(probe2symbol2$ID)
dat2<-dat2 %>%
  inner_join(probe2symbol2,by='ID')%>% 
  select(-ID)%>%     ## 去除多余信息
  select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('MMRF_',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
write.table(dat2,file = 'dat_tcga.xls',sep = '\t',row.names = T,quote = F)
## 01-3 验证集GSE4581----------
gset_va<-getGEO("GSE4581",
                destdir = '.',
                GSEMatrix = T,
                getGPL = F)
expr_va<-as.data.frame(exprs(gset_va[[1]]))
a_va=gset_va[[1]]
pd_va<-pData(a_va)
gpl2<-getGEO("GPL570",destdir = '.')
gpl2<-Table(gpl2)    
colnames(gpl2)
probe2symobl2<-gpl2 %>%
  select('ID','Gene Symbol')%>%
  filter('Gene Symbol'!='')%>%
  separate('Gene Symbol',c('symbol','drop'),sep = '//')%>%
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
survival_va<-data.frame(sample=pd_va$geo_accession,OS=pd_va$characteristics_ch1,OS.time=pd_va$characteristics_ch1.2)
survival_va$OS<-gsub('[SURIND=','',survival_va$OS,fixed = T)
survival_va$OS<-gsub(' (Indicator of disease-related death; integer, 0=alive or death by other cause, 1=disease related death, na=death cause undetermined)]','',survival_va$OS,fixed = T)
survival_va$OS.time<-gsub('[SURTIM=','',survival_va$OS.time,fixed = T)
survival_va$OS.time<-gsub(' (Follow-up time in months from Pre-Treatment baseline; integer)]','',survival_va$OS.time,fixed = T)
survival_va$OS<-as.numeric(survival_va$OS)
survival_va$OS.time<-as.numeric(survival_va$OS.time)*30
write.table(dat_va,file = 'exp_va.xls',sep = '\t',quote = F,row.names = T)
write.table(survival_va,file = 'survival_va.xls',sep = '\t',row.names = F,quote = F)
