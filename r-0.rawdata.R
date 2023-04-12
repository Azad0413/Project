rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/SJZZK-431-10/")
if (! dir.exists("./00_rawdata")){
  dir.create("./00_rawdata")
}
setwd("./00_rawdata")
library(tidyverse)
library(lance)
library(data.table)
exp <- fread('data_exp_icc.txt',header = F)%>%column_to_rownames(var = 'V1')
expr <- exp
library(readxl)
library(GEOquery)
group <- read_xlsx('group.xlsx')
colnames(expr) <- group$sample

gpl<-getGEO("GPL17585",destdir = '.')
gpl<-Table(gpl)    
colnames(gpl)
gpl$gene_symbols <- gsub('|','//',gpl$gene_symbols,fixed = T)
probe2symbol<-gpl %>%
  dplyr::select('ID','gene_symbols')%>%
  filter('gene_symbols'!='')%>%
  separate('gene_symbols',c('symbol','drop'),sep = '//')%>%
  dplyr::select(-drop)
probe2symbol=probe2symbol[probe2symbol$symbol!='---',]
probe2symbol <- na.omit(probe2symbol)
dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symbol$ID<-as.character(probe2symbol$ID)
dat<-dat %>%
  inner_join(probe2symbol,by='ID')%>% 
  dplyr::select(-ID)%>%     ## 去除多余信息
  dplyr::select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('ICC',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
write.table(dat,file = 'dat.xls',sep = '\t',row.names = T,quote = F)

