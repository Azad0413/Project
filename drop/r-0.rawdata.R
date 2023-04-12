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
gset<-getGEO("GSE182616",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
a=gset[[1]]
pd<-pData(a)
gpl<-getGEO("GPL17077",destdir = '.')
gpl<-Table(gpl)    
colnames(gpl)
probe2symbol<-gpl %>%
  select('ID','GENE_SYMBOL')%>%
  filter('GENE_SYMBOLl'!='')%>%
  separate('GENE_SYMBOL',c('symbol','drop'),sep = '///')%>%
  select(-drop)
probe2symbol=probe2symbol[probe2symbol$symbol!='',]
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
### 提取0、2、4、8、12、24h的表达数据；
colnames(pd)
clinical<-data.frame(Sample=pd$geo_accession,
                     time=pd$characteristics_ch1,
                   #  statu=pd$`mortality:ch1`,
                     TBSA=pd$`tbsa:ch1`)
table(clinical$time)
clinical$time<-ifelse(clinical$time=='time point: hr0','0',
                      ifelse(clinical$time=='time point: hr2','2',
                             ifelse(clinical$time=='time point: hr4','4',
                                    ifelse(clinical$time=='time point: hr8','8',
                                           ifelse(clinical$time=='time point: hr12','12',
                                                  ifelse(clinical$time=='time point: hr24','24',NA))))))
table(clinical$time)
clinical$time<-as.numeric(clinical$time)
clinical$TBSA<-ifelse(clinical$TBSA<20,'Low','High')
clinical<-clinical[order(clinical$time),]
clinical<-clinical[c(1:305),]
## 把表达矩阵提取出来
data.final<-dat[,clinical$Sample]
write.table(clinical,file = 'clinical.xls',sep = '\t',quote = F,row.names = F)
write.table(data.final,file = 'data.final.xls',sep = '\t',quote = F,row.names = T)



## STEM需要准备有时间梯度的数据。第一列是基因名，之后是表达信息，需要按照时间序列排序。
## STEM对于重复可以预先处理为平均值或者按照重复数据多次导入（中位数处理）
stem.dead<-data.frame('0h'=rowMeans(dead.dat[,dead.sample$time=='0'],na.rm = T),
                      '2h'=rowMeans(dead.dat[,dead.sample$time=='2'],na.rm = T),
                      '4h'=rowMeans(dead.dat[,dead.sample$time=='4'],na.rm = T),
                      '8h'=rowMeans(dead.dat[,dead.sample$time=='8'],na.rm = T),
                      '12h'=rowMeans(dead.dat[,dead.sample$time=='12'],na.rm = T),
                      '24h'=rowMeans(dead.dat[,dead.sample$time=='24'],na.rm = T))

colnames(stem.dead)<-c('0h','2h','4h','8h','12h','24h')
write.table(stem.dead,file = 'dead_dat.xls',
            sep = '\t',
            row.names = T)
stem.alive<-data.frame('0h'=rowMeans(alive.dat[,alive.sample$time=='0'],na.rm = T),
                      '2h'=rowMeans(alive.dat[,alive.sample$time=='2'],na.rm = T),
                      '4h'=rowMeans(alive.dat[,alive.sample$time=='4'],na.rm = T),
                      '8h'=rowMeans(alive.dat[,alive.sample$time=='8'],na.rm = T),
                      '12h'=rowMeans(alive.dat[,alive.sample$time=='12'],na.rm = T),
                      '24h'=rowMeans(alive.dat[,alive.sample$time=='24'],na.rm = T))

colnames(stem.alive)<-c('0h','2h','4h','8h','12h','24h')
write.table(stem.alive,file = 'alive_dat.xls',
            sep = '\t',
            row.names = T)
# 02 STEM时序分析----------
setwd("/data/nas1/luchunlin/project/BJTC-302")
if (! dir.exists("./01_STEM")){
  dir.create("./01_STEM")
}
setwd("./01_STEM")
## 02-1 死亡时序分析---------
## 39 18 44 24 13 16 21 8 48 47 25 27 30 40 49 2 
dead.profile<-read_xlsx('dead_GENE.xlsx')
dead.sig<-c('39','18','44','24','13','16','21','8','48','47','25','27','30','40','49','2')
dead.profile<-dead.profile[dead.profile$Profile%in%dead.sig,]
## 4495
write.table(dead.profile,file = 'dead_profile.xls',
            sep = '\t',
            row.names = T)
## 02-2 存活时序分析---------
## 8 39 40 21 0 36 41 3 13 16 
alive.profile<-read_xlsx('alive_GENE.xlsx')
alive.sig<-c('8','39','40','21','0','36','41','3','13','16')
alive.profile<-alive.profile[alive.profile$Profile%in%alive.sig,]
## 2228
write.table(alive.profile,file = 'alive_profile.xls',
            sep = '\t',
            row.names = T)

# 03 GO/KEGG----------
setwd("/data/nas1/luchunlin/project/BJTC-302")
if (! dir.exists("./02_GO_KEGG")){
  dir.create("./02_GO_KEGG")
}
setwd("./02_GO_KEGG")
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(GOplot)
## 03-1 DEAD kegg----------
### 39 
dead.39 <- bitr(dead.profile$Symbol[dead.profile$Profile=='39'],
                        fromType = "SYMBOL",
                        toType = c("ENTREZID"),
                        OrgDb = "org.Hs.eg.db")

## KEGG富集分析（气泡图）
kk.39 <- enrichKEGG(gene = dead.39$ENTREZID,
                 keyType = "kegg",
                 organism = "hsa",
                 pAdjustMethod = "none",
                 pvalueCutoff = 0.05)
kk.39 <- setReadable(kk.39, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk.39@result,file = "dead.39.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk.39,showCategory=20)
kk_dot
### 18 
dead.18 <- bitr(dead.profile$Symbol[dead.profile$Profile=='18'],
                fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = "org.Hs.eg.db")
kk.18 <- enrichKEGG(gene = dead.18$ENTREZID,
                    keyType = "kegg",
                    organism = "hsa",
                    pAdjustMethod = "none",
                    pvalueCutoff = 0.05)
kk.18 <- setReadable(kk.18, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk.18@result,file = "dead.18.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk.18,showCategory=20)
kk_dot
### 44 
dead.44 <- bitr(dead.profile$Symbol[dead.profile$Profile=='44'],
                fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = "org.Hs.eg.db")
kk.44 <- enrichKEGG(gene = dead.44$ENTREZID,
                    keyType = "kegg",
                    organism = "hsa",
                    pAdjustMethod = "none",
                    pvalueCutoff = 0.05)
kk.44 <- setReadable(kk.44, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk.44@result,file = "dead.44.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk.44,showCategory=20)
kk_dot
### 24 
dead.24 <- bitr(dead.profile$Symbol[dead.profile$Profile=='24'],
                fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = "org.Hs.eg.db")
kk.24 <- enrichKEGG(gene = dead.24$ENTREZID,
                    keyType = "kegg",
                    organism = "hsa",
                    pAdjustMethod = "none",
                    pvalueCutoff = 0.05)
kk.24 <- setReadable(kk.24, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk.24@result,file = "dead.24.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk.24,showCategory=20)
kk_dot
### 13 
dead.13 <- bitr(dead.profile$Symbol[dead.profile$Profile=='13'],
                fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = "org.Hs.eg.db")
kk.13 <- enrichKEGG(gene = dead.13$ENTREZID,
                    keyType = "kegg",
                    organism = "hsa",
                    pAdjustMethod = "none",
                    pvalueCutoff = 0.05)
kk.13 <- setReadable(kk.13, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk.13@result,file = "dead.13.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk.13,showCategory=20)
kk_dot
### 16 
dead.16 <- bitr(dead.profile$Symbol[dead.profile$Profile=='16'],
                fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = "org.Hs.eg.db")
kk.16 <- enrichKEGG(gene = dead.16$ENTREZID,
                    keyType = "kegg",
                    organism = "hsa",
                    pAdjustMethod = "none",
                    pvalueCutoff = 0.05)
kk.16 <- setReadable(kk.16, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk.16@result,file = "dead.16.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk.16,showCategory=20)
kk_dot
### 21 
dead.21 <- bitr(dead.profile$Symbol[dead.profile$Profile=='21'],
                fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = "org.Hs.eg.db")
kk.21 <- enrichKEGG(gene = dead.21$ENTREZID,
                    keyType = "kegg",
                    organism = "hsa",
                    pAdjustMethod = "none",
                    pvalueCutoff = 0.05)
kk.21 <- setReadable(kk.21, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk.21@result,file = "dead.21.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk.21,showCategory=20)
kk_dot
### 8 
dead.8 <- bitr(dead.profile$Symbol[dead.profile$Profile=='8'],
                fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = "org.Hs.eg.db")
kk.8 <- enrichKEGG(gene = dead.8$ENTREZID,
                    keyType = "kegg",
                    organism = "hsa",
                    pAdjustMethod = "none",
                    pvalueCutoff = 0.05)
kk.8 <- setReadable(kk.8, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk.8@result,file = "dead.8.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk.8,showCategory=20)
kk_dot
### 48 
dead.48 <- bitr(dead.profile$Symbol[dead.profile$Profile=='48'],
                fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = "org.Hs.eg.db")
kk.48 <- enrichKEGG(gene = dead.48$ENTREZID,
                    keyType = "kegg",
                    organism = "hsa",
                    pAdjustMethod = "none",
                    pvalueCutoff = 0.05)
kk.48 <- setReadable(kk.48, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk.48@result,file = "dead.48.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk.48,showCategory=20)
kk_dot
### 47 
dead.47 <- bitr(dead.profile$Symbol[dead.profile$Profile=='47'],
                fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = "org.Hs.eg.db")
kk.47 <- enrichKEGG(gene = dead.47$ENTREZID,
                    keyType = "kegg",
                    organism = "hsa",
                    pAdjustMethod = "none",
                    pvalueCutoff = 0.05)
kk.47 <- setReadable(kk.47, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk.47@result,file = "dead.47.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk.47,showCategory=20)
kk_dot
### 25 
dead.25 <- bitr(dead.profile$Symbol[dead.profile$Profile=='25'],
                fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = "org.Hs.eg.db")
kk.25 <- enrichKEGG(gene = dead.25$ENTREZID,
                    keyType = "kegg",
                    organism = "hsa",
                    pAdjustMethod = "none",
                    pvalueCutoff = 0.05)
kk.25 <- setReadable(kk.25, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk.25@result,file = "dead.25.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk.25,showCategory=20)
kk_dot
### 27 
dead.27 <- bitr(dead.profile$Symbol[dead.profile$Profile=='27'],
                fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = "org.Hs.eg.db")
kk.27 <- enrichKEGG(gene = dead.27$ENTREZID,
                    keyType = "kegg",
                    organism = "hsa",
                    pAdjustMethod = "none",
                    pvalueCutoff = 0.05)
kk.27 <- setReadable(kk.27, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk.27@result,file = "dead.27.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk.27,showCategory=20)
kk_dot
### 30 
dead.30 <- bitr(dead.profile$Symbol[dead.profile$Profile=='30'],
                fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = "org.Hs.eg.db")
kk.30 <- enrichKEGG(gene = dead.30$ENTREZID,
                    keyType = "kegg",
                    organism = "hsa",
                    pAdjustMethod = "none",
                    pvalueCutoff = 0.05)
kk.30 <- setReadable(kk.30, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk.30@result,file = "dead.30.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk.30,showCategory=20)
kk_dot
### 40 
dead.40 <- bitr(dead.profile$Symbol[dead.profile$Profile=='40'],
                fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = "org.Hs.eg.db")
kk.40 <- enrichKEGG(gene = dead.40$ENTREZID,
                    keyType = "kegg",
                    organism = "hsa",
                    pAdjustMethod = "none",
                    pvalueCutoff = 0.05)
kk.40 <- setReadable(kk.40, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk.40@result,file = "dead.40.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk.40,showCategory=20)
kk_dot
### 49 
dead.49 <- bitr(dead.profile$Symbol[dead.profile$Profile=='49'],
                fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = "org.Hs.eg.db")
kk.49 <- enrichKEGG(gene = dead.49$ENTREZID,
                    keyType = "kegg",
                    organism = "hsa",
                    pAdjustMethod = "none",
                    pvalueCutoff = 0.05)
kk.49 <- setReadable(kk.49, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk.49@result,file = "dead.49.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk.49,showCategory=20)
kk_dot
### 2
dead.2 <- bitr(dead.profile$Symbol[dead.profile$Profile=='2'],
                fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = "org.Hs.eg.db")
kk.2 <- enrichKEGG(gene = dead.2$ENTREZID,
                    keyType = "kegg",
                    organism = "hsa",
                    pAdjustMethod = "none",
                    pvalueCutoff = 0.05)
kk.2 <- setReadable(kk.2, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk.2@result,file = "dead.2.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk.2,showCategory=20)
kk_dot
## 03-2 ALIVE kegg----------
### 8 
alive.8 <- bitr(alive.profile$Symbol[alive.profile$Profile=='8'],
                fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = "org.Hs.eg.db")

## KEGG富集分析（气泡图）
kk.8 <- enrichKEGG(gene = alive.8$ENTREZID,
                    keyType = "kegg",
                    organism = "hsa",
                    pAdjustMethod = "none",
                    pvalueCutoff = 0.05)
kk.8 <- setReadable(kk.8, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk.8@result,file = "alive.8.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk.8,showCategory=20)
kk_dot

### 39
alive.39 <- bitr(alive.profile$Symbol[alive.profile$Profile=='39'],
                fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = "org.Hs.eg.db")

## KEGG富集分析（气泡图）
kk.39 <- enrichKEGG(gene = alive.39$ENTREZID,
                   keyType = "kegg",
                   organism = "hsa",
                   pAdjustMethod = "none",
                   pvalueCutoff = 0.05)
kk.39 <- setReadable(kk.39, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk.39@result,file = "alive.39.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk.39,showCategory=20)
kk_dot
### 40 
alive.40 <- bitr(alive.profile$Symbol[alive.profile$Profile=='40'],
                fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = "org.Hs.eg.db")

## KEGG富集分析（气泡图）
kk.40 <- enrichKEGG(gene = alive.40$ENTREZID,
                   keyType = "kegg",
                   organism = "hsa",
                   pAdjustMethod = "none",
                   pvalueCutoff = 0.05)
kk.40 <- setReadable(kk.40, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk.40@result,file = "alive.40.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk.40,showCategory=20)
kk_dot
### 21 
alive.21 <- bitr(alive.profile$Symbol[alive.profile$Profile=='21'],
                fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = "org.Hs.eg.db")

## KEGG富集分析（气泡图）
kk.21 <- enrichKEGG(gene = alive.21$ENTREZID,
                   keyType = "kegg",
                   organism = "hsa",
                   pAdjustMethod = "none",
                   pvalueCutoff = 0.05)
kk.21 <- setReadable(kk.21, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk.21@result,file = "alive.21.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk.21,showCategory=20)
kk_dot
### 0 
alive.0 <- bitr(alive.profile$Symbol[alive.profile$Profile=='0'],
                fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = "org.Hs.eg.db")

## KEGG富集分析（气泡图）
kk.0 <- enrichKEGG(gene = alive.0$ENTREZID,
                   keyType = "kegg",
                   organism = "hsa",
                   pAdjustMethod = "none",
                   pvalueCutoff = 0.05)
kk.0 <- setReadable(kk.0, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk.0@result,file = "alive.0.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk.0,showCategory=20)
kk_dot
### 36 
alive.36 <- bitr(alive.profile$Symbol[alive.profile$Profile=='36'],
                fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = "org.Hs.eg.db")

## KEGG富集分析（气泡图）
kk.36 <- enrichKEGG(gene = alive.36$ENTREZID,
                   keyType = "kegg",
                   organism = "hsa",
                   pAdjustMethod = "none",
                   pvalueCutoff = 0.05)
kk.36 <- setReadable(kk.36, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk.36@result,file = "alive.36.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk.36,showCategory=20)
kk_dot
### 41 
alive.41 <- bitr(alive.profile$Symbol[alive.profile$Profile=='41'],
                fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = "org.Hs.eg.db")

## KEGG富集分析（气泡图）
kk.41 <- enrichKEGG(gene = alive.41$ENTREZID,
                   keyType = "kegg",
                   organism = "hsa",
                   pAdjustMethod = "none",
                   pvalueCutoff = 0.05)
kk.41 <- setReadable(kk.41, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk.41@result,file = "alive.41.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk.41,showCategory=20)
kk_dot
### 3 
alive.3 <- bitr(alive.profile$Symbol[alive.profile$Profile=='3'],
                fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = "org.Hs.eg.db")

## KEGG富集分析（气泡图）
kk.3 <- enrichKEGG(gene = alive.3$ENTREZID,
                   keyType = "kegg",
                   organism = "hsa",
                   pAdjustMethod = "none",
                   pvalueCutoff = 0.05)
kk.3 <- setReadable(kk.3, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk.3@result,file = "alive.3.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk.3,showCategory=20)
kk_dot
### 13 
alive.13 <- bitr(alive.profile$Symbol[alive.profile$Profile=='13'],
                fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = "org.Hs.eg.db")

## KEGG富集分析（气泡图）
kk.13 <- enrichKEGG(gene = alive.13$ENTREZID,
                   keyType = "kegg",
                   organism = "hsa",
                   pAdjustMethod = "none",
                   pvalueCutoff = 0.05)
kk.13 <- setReadable(kk.13, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk.13@result,file = "alive.13.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk.13,showCategory=20)
kk_dot
### 16 
alive.16 <- bitr(alive.profile$Symbol[alive.profile$Profile=='16'],
                fromType = "SYMBOL",
                toType = c("ENTREZID"),
                OrgDb = "org.Hs.eg.db")

## KEGG富集分析（气泡图）
kk.16 <- enrichKEGG(gene = alive.16$ENTREZID,
                   keyType = "kegg",
                   organism = "hsa",
                   pAdjustMethod = "none",
                   pvalueCutoff = 0.05)
kk.16 <- setReadable(kk.16, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk.16@result,file = "alive.16.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk.16,showCategory=20)
kk_dot

## 交集基因
intersec<-dead.profile[dead.profile$Symbol%in%alive.profile$Symbol,]
### 461

## 并集基因 去重
sum<-rbind(dead.profile,alive.profile)
sum<-sum[!duplicated(sum$Symbol),]
## 6262