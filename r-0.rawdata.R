rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/BJTC-302")
if (! dir.exists("./00_rawdata")){
  dir.create("./00_rawdata")
}
setwd("./00_rawdata")

library(GEOquery)
library(Biobase)
library(tidyverse)
library(dplyr)
gset<-getGEO("GSE77791",
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
## 30 个brun治疗前。13个健康
group.all<-data.frame(sample=pd$geo_accession,group=pd$description)
mete=data.frame(group.all$group)
for (i in 1:length(mete[,1])) {
  num=as.numeric(as.character(substring(mete[i,1],8,8)))
  if(num %in% seq(1)){mete[i,2]="T"}
  if(num %in% seq(2,4)){mete[i,2]="N"}
}
names(mete)=c("sample","group")
mete[101,2]<-'T'
mete$group=as.factor(mete$group)
mete=subset(mete,mete$group=="T")
brun.sample<-group.all[group.all$group%in%mete$sample,]
brun.sample$group<-'T'
##30
control.sample<-group.all[c(105:117),]
group<-rbind(brun.sample,control.sample)
table(group$group)
group$group<-ifelse(group$group=='T','Burn','control')
write.table(group,file = 'group.xls',sep = '\t',row.names = F,quote = F)
dat.final<-dat[,group$sample]
dat.final<-dat.final[-4,]
write.table(dat.final,'dat.xls',sep = '\t',row.names = T,quote = F)
pd<-pd[pd$geo_accession%in%colnames(dat.final),]
# phenotype<-data.frame(sample=pd$geo_accession,tbsa=pd$`tbsa:ch1`)
# write.table(phenotype,'tbsa.xls',sep = '\t',row.names = F,quote = F)
survival<-data.frame(sample=pd$geo_accession,OS=pd$`survival (d28):ch1`)
write.table(survival,file = 'survival.xls',sep = '\t',row.names = F,quote = F)
