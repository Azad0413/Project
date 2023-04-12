rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-386-10/")
if (! dir.exists("./03_m7Glnc")){
  dir.create(".//03_m7Glnc")
}
setwd(".//03_m7Glnc")
# pearson
library(lance)
library(tidyverse)
library(Ipaper)
m7G <- read_xlsx('01.m7G_geneset.xlsx')
dat = read.delim2("../00_rawdata/dat.fpkm.xls", row.names = 1) %>% lc.tableToNum
colnames(dat)<-gsub('.','-',colnames(dat),fixed = T)
dat.all<-read.delim2("../00_rawdata/dat.all(fpkm).xls", row.names = 1) %>% lc.tableToNum
colnames(dat.all)<-gsub('.','-',colnames(dat.all),fixed = T)
#dat<-dat[,colnames(dat1)]
dat.m7G<-dat.all[m7G$symbol,]
dat.m7G<-dat.all[rownames(dat.all)%in%m7G$symbol,]
## NCBP3 	NUDT4B
## C17orf85  NUDT4P1
unmap <- dat.all[c('C17orf85','NUDT4P1'),]
rownames(unmap) <- c('NCBP3','NUDT4B')
dat.m7G <- rbind(dat.m7G,unmap)
write.table(dat.m7G,file = '02.dat.m7G.xls',sep = '\t',row.names = T,quote = F)
library(Hmisc)
## 氧化应激相关基因进行相关性分析（cor>0.4），
corr.dat<-t(rbind(dat,dat.m7G))
#基因表达值的相关性分析，以Pearson相关系数为例
gene_cor <- rcorr(corr.dat)$r[1:nrow(dat),(ncol(corr.dat)-length(rownames(dat.m7G))+1):ncol(corr.dat)]

#将获得的相关性矩阵转换为两两对应的数据框结构
gene_cor <- reshape2::melt(gene_cor)
#gene_cor <- subset(gene_cor, value != 0)  #去除0值的相关性
head(gene_cor)  #前两列是两个基因名称，第三列为两个基因的相关性
colnames(gene_cor) <- c('lncRNA','m7G','correlation')
## 计算p值
gene_p <- rcorr(corr.dat)$P[1:nrow(dat),(ncol(corr.dat)-length(rownames(dat.m7G))+1):ncol(corr.dat)]
gene_p <- reshape2::melt(gene_p)
head(gene_p)  #前两列是两个基因名称，第三列为
gene_cor$pvalu <- gene_p$value
write.table(gene_cor,file = '03.correlation.all.xls',sep = '\t',row.names = F,quote = F)
##|Cor| > 0.4, P< 0.001
cor.final <- subset(gene_cor,gene_cor$pvalu<0.001 & abs(gene_cor$correlation)>0.6)
## 62
write.table(cor.final,file = '04.correlation.final.xls',sep = '\t',row.names = F,quote = F)
