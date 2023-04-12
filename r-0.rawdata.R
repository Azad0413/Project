rm(list = ls())
# 设置当前工作环境 ---------------------------------------------------------------
setwd("/data/nas1/luchunlin/project/BJTC-334")
if (! dir.exists("./00_rawdata")){
  dir.create("./00_rawdata")
}
setwd("./00_rawdata")
library(GEOquery)    ##加载包
gset <- getGEO('GSE23561',   ##GSE73094
               AnnotGPL = F,
               getGPL = F)
#save(gset,file = '00_Rawdata/GSE23561.gset.Rdata')   
eset<-exprs(gset[[1]])   ##表达矩阵  
metadata<-pData(gset[[1]]) ##临床信息

# ID转换 --------------------------------------------------------------------
gpl <- getGEO('GPL10775')    ##数据集对应的平台文件，若无法下载可直接去GEO官网下载
anno1<-Table(gpl)[,c(1,5)]    #提取id和gene symbol 
colnames(anno1)<-c('ID','Gene Symbol')  
head(anno1)
eset<-as.data.frame(eset)  
eset$ID <- rownames(eset)   
head(eset)
merg<-merge(eset,anno1,by="ID")  ##平台ID和表达量矩阵ID合并
y<-merg$`Gene Symbol`
gene<-unlist(lapply(y,function(y) strsplit(as.character(y)," /// ")[[1]][1]))
merg$gene <- gene  
aggr<-aggregate(merg[,2:36],by=list(merg$gene),mean)   ##371为样本数+1
aggr<-aggr[!duplicated(aggr$Group.1),]
rownames(aggr)<-aggr[,1]    ###以基因名为列名 
aggr<-aggr[,-1]       
expr<-aggr
condition<-data.frame(sample=rownames(metadata),group=metadata$characteristics_ch1)
condition<-subset(condition,group=='disease state: Control'|group=='disease state: Coronary Artery Disease')
condition$group<-ifelse(condition$group=='disease state: Control','Control','CAD')
group<-data.frame(sample=condition$sample,group=condition$group)
expr<-subset(expr,select=group$sample)
mRNA<-read.table('mRNA.txt')
c<-intersect(rownames(expr),mRNA$gene_name)
dat.final = expr[c,]   ###提取mrna  ##13717个基因  15个样本（正常9，患病6）
write.table(dat.final,file = 'dat.final.xls',quote = F,sep = '\t',row.names = T)

write.table(group,file = 'group.xls',quote = F,sep = '\t',row.names = F)
# 
# # 01 获取数据集----------
# setwd("/data/nas1/luchunlin/project/BJTC-334")
# if (! dir.exists("./00_rawdata")){
#   dir.create("./00_rawdata")
# }
# setwd("./00_rawdata")
# ### GEO数据库 GSE23561
# library(GEOquery)
# library(Biobase)
# library(tidyr)
# library(AnnoProbe)
# library(tidyverse)
# gset<-getGEO("GSE23561",destdir = '.',GSEMatrix = T,getGPL = F)
# expr<-as.data.frame(exprs(gset[[1]]))
# a=gset[[1]]
# gpl<-getGEO("GPL10775",destdir = '.')
# gpl<-Table(gpl)  
# colnames(gpl)
# gpl$`Symbol v12`<-gsub('.','',gpl$`Symbol v12`,fixed = T)
# probe2symbol<-gpl %>%
#   select('ID','Symbol v12')%>%
#   filter('Symbol v12'!='')%>%
#   separate('Symbol v12',c('symbol','drop'),sep = '///')%>%
#   select(-drop)
# length(unique(probe2symbol$symbol))
# #probe2symbol$symbol<-gsub('.','',probe2symbol$symbol,fixed = T)
# probe2symbol<-probe2symbol[!probe2symbol$symbol=='',]
# dat<-expr
# dat$ID<-rownames(dat)
# dat$ID<-as.character(dat$ID)
# probe2symbol$ID<-as.character(probe2symbol$ID)
# dat<-dat %>%
#   inner_join(probe2symbol,by='ID')%>% 
#   select(-ID)%>%     ## 去除多余信息
#   select(symbol,everything())%>%     ## 重新排列
#   # mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
#   # arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
#   distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
#   # select(-rowMean)%>%     ## 反向选择去除rowMean这一列
#   tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
# ## 6个冠心病和9个对照样本
# pd<-pData(a)
# group<-data.frame(sample=pd$geo_accession,
#                   group=pd$title)
# group<-group[c(1:9,22:27),]
# group$group<-c(rep('control',9),rep('CAD',6))
# dat.final<-dat[,colnames(dat)%in%group$sample]
# write.table(dat.final,file = 'dat.final.xls',quote = F,sep = '\t',row.names = T)
# write.table(group,file = 'group.xls',quote = F,sep = '\t',row.names = F)
