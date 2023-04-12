rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/CD-0601-2/")
if (! dir.exists("./00_rawdata")){
  dir.create("./00_rawdata")
}
setwd("./00_rawdata")
library(GEOquery)
library(tidyverse)
library(lance)
### GSE65682-------
gset<-getGEO("GSE65682",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL13667",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)    
colnames(gpl)
probe2symobl<-gpl %>%
  dplyr::select('ID','Gene Symbol')%>%
  filter('Gene Symbol'!='')%>%
  separate('Gene Symbol',c('symbol','drop'),sep = ' /// ')%>%
  dplyr::select(-drop)
probe2symobl=probe2symobl[probe2symobl$symbol!='',]
probe2symobl=probe2symobl[probe2symobl$symbol!='---',]

dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symobl$ID<-as.character(probe2symobl$ID)
dat<-dat %>%
  inner_join(probe2symobl,by='ID')%>% 
  dplyr::select(-ID)%>%     ## 去除多余信息
  dplyr::select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
# pcg <- read.delim2('/data/nas1/luchunlin/pipeline/PCG/PCG.xls(v22)')
# dat <- dat[pcg$gene_name,]
# dat <- na.omit(dat)
pd <- pData(a)
table(pd$`icu_acquired_infection:ch1`)
group <- data.frame(sample=pd$geo_accession,group=pd$`icu_acquired_infection:ch1`)
table(group$group)
group$group <- ifelse(group$group=='healthy','control','Sepsis')

dat<-dat[,group$sample]

write.table(dat,file = 'dat(GSE65682).xls',sep = '\t',row.names = T,quote = F)
write.table(group,file = 'group(GSE65682).xls',sep = '\t',row.names = F,quote = F)


###GSE28750------
gset<-getGEO("GSE28750",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL570",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)    
colnames(gpl)
probe2symobl<-gpl %>%
  dplyr::select('ID','Gene Symbol')%>%
  filter('Gene Symbol'!='')%>%
  separate('Gene Symbol',c('symbol','drop'),sep = '///')%>%
  dplyr::select(-drop)
probe2symobl=probe2symobl[probe2symobl$symbol!='',]

dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symobl$ID<-as.character(probe2symobl$ID)
dat<-dat %>%
  inner_join(probe2symobl,by='ID')%>% 
  dplyr::select(-ID)%>%     ## 去除多余信息
  dplyr::select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
# pcg <- read.delim2('/data/nas1/luchunlin/pipeline/PCG/PCG.xls(v22)')
# dat <- dat[pcg$gene_name,]
# dat <- na.omit(dat)
pd <- pData(a)
table(pd$`health status:ch1`)
group <- data.frame(sample=pd$geo_accession,group=pd$`health status:ch1`)
table(group$group)
group <- subset(group,group=='HEALTHY'|group=='SEPSIS')
group$group <- ifelse(group$group=='HEALTHY','control','Sepsis')

dat<-dat[,group$sample]

write.table(dat,file = 'dat(GSE28750).xls',sep = '\t',row.names = T,quote = F)
write.table(group,file = 'group(GSE28750).xls',sep = '\t',row.names = F,quote = F)
