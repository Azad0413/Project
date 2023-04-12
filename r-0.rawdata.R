rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/HF-0106-2/")
if (! dir.exists("./00_rawdata")){
  dir.create("./00_rawdata")
}
setwd("./00_rawdata")
library(GEOquery)
library(tidyverse)
### -------
gset<-getGEO("GSE169568",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL10558",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)    
colnames(gpl)
gpl <- na.omit(gpl)

probe2symobl<-gpl %>%
  dplyr::select('ID','Symbol')%>%
  filter('Symbol'!='')%>%
  separate('Symbol',c('symbol','drop'),sep = '///')%>%
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

pd <- pData(a)
table(pd$`diagnosis:ch1`)
group<-data.frame(sample=pd$geo_accession,group=pd$`diagnosis:ch1`)
table(group$group)
group <- subset(group,group=='Healthy control'|group=='Ulcerative colitis')

group$group<-ifelse(group$group=='Healthy control','control','UC')
group<-group[order(group$group),]
dat<-dat[,group$sample]

write.table(dat,file = 'dat(GSE169568).xls',sep = '\t',row.names = T,quote = F)
write.table(group,file = 'group(GSE169568).xls',sep = '\t',row.names = F,quote = F)

### -------
gset<-getGEO("GSE94648",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL19109",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)    
colnames(gpl)
gpl <- na.omit(gpl)

probe2symobl<-gpl %>%
  dplyr::select('ID','ENTREZ_GENE_ID')%>%
  filter('ENTREZ_GENE_ID'!='')
write.table(probe2symobl,file = 'id.xls',sep = '\t',row.names = F,quote = F)
probe2symobl <- read.delim2('probe2symbol.xls')
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

pd <- pData(a)
table(pd$`case_phenotype:ch1`)
group<-data.frame(sample=pd$geo_accession,group=pd$`case_phenotype:ch1`)
table(group$group)
group <- subset(group,group=='Control'|group=='Colitis')

group$group<-ifelse(group$group=='Control','control','UC')
group<-group[order(group$group),]
dat<-dat[,group$sample]

write.table(dat,file = 'dat(GSE94648).xls',sep = '\t',row.names = T,quote = F)
write.table(group,file = 'group(GSE94648).xls',sep = '\t',row.names = F,quote = F)



### -------
gset<-getGEO("GSE186507",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-read.table('GSE186507_MSCCR_Blood_adjcounts.txt.gz')

id <- data.frame(rownames(expr))
write.table(id,file = 'id.xls',sep = '\t',row.names = F,quote = F)

probe2symobl <- read.delim2('probe2symbol.xls')
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

pd <- pData(gset[[1]])
table(pd$`diagnosis:ch1`)
group<-data.frame(sample=pd$geo_accession,group=pd$`ibd_disease:ch1`,title=pd$title)
table(group$group)
group <- subset(group,group=='Control'|group=='UC')

group$group<-ifelse(group$group=='Control','control','UC')
group<-group[order(group$group),]
group$title <- gsub(', UC participants','',group$title)
group$title <- gsub(', Control participants','',group$title)

dat<-dat[,group$title]
colnames(dat) <- group$sample

write.table(dat,file = 'dat(GSE186507).xls',sep = '\t',row.names = T,quote = F)
write.table(group,file = 'group(GSE186507).xls',sep = '\t',row.names = F,quote = F)

### -------
gset<-getGEO("GSE119600",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))

gpl<-getGEO("GPL10558",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)    
colnames(gpl)
gpl <- na.omit(gpl)

probe2symobl<-gpl %>%
  dplyr::select('ID','Symbol')%>%
  filter('Symbol'!='')%>%
  separate('Symbol',c('symbol','drop'),sep = '///')%>%
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

pd <- pData(a)
table(pd$`condition:ch1`)
group<-data.frame(sample=pd$geo_accession,group=pd$`condition:ch1`)
table(group$group)
group <- subset(group,group=='control'|group=='ulcerative colitis')

group$group<-ifelse(group$group=='control','control','UC')
group<-group[order(group$group),]
dat<-dat[,group$sample]

write.table(dat,file = 'dat(GSE119600).xls',sep = '\t',row.names = T,quote = F)
write.table(group,file = 'group(GSE119600).xls',sep = '\t',row.names = F,quote = F)


### -------
gset<-getGEO("GSE126124",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL6244",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)    
colnames(gpl)
library(AnnoProbe)
gpl <- idmap(gpl = 'GPL6244')

probe2symobl<-gpl 
colnames(probe2symobl) <- c('ID','symbol')
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

pd <- pData(a)

table(pd$`tissue:ch1`)
pd <- pd[which(pd$`tissue:ch1`=='peripheral whole blood'),]
table(pd$`disease type:ch1`)

group<-data.frame(sample=pd$geo_accession,group=pd$`disease type:ch1`)
table(group$group)
group <- subset(group,group=='Control'|group=='Control - celiac'|group=='Control - infect'|group=='Control (E.coli) --> CD'|group=='Control--> CD'|
                  group=='UC - L colitis'|group=='UC - pancolitis'|group=='UC - proctitis')
# group <- subset(group,group=='Control'|
#                   group=='UC - L colitis'|group=='UC - pancolitis'|group=='UC - proctitis')
table(group$group)
group <- group[order(group$group),]
group$group <- c(rep('control',39),rep('UC',18))
# group$group<-ifelse(group$group=='Healthy control','control','UC')
# group<-group[order(group$group),]
dat<-dat[,group$sample]

write.table(dat,file = 'dat(GSE126124).xls',sep = '\t',row.names = T,quote = F)
write.table(group,file = 'group(GSE126124).xls',sep = '\t',row.names = F,quote = F)

