rm(list = ls())
setwd("/data/nas1/luchunlin/project/SJZZK-431-10/")
if (! dir.exists("./05_PPI")){
  dir.create("./05_PPI")
}
setwd("./05_PPI")

###做一个相关性
##选择|cor|>0.3, P<0.05作为筛选阈值，构建共表达网络

# pearson
library(lance)
library(tidyverse)
library(Ipaper)
geneset <- read.delim2('../03_DEAKG/DEAKGs.xls')
dat = read.delim2("../00_rawdata/dat.xls", row.names = 1) %>% lc.tableToNum

#dat<-dat[,colnames(dat1)]
dat.gene<-dat[geneset$symbol,]

library(Hmisc)
## 氧化应激相关基因进行相关性分析（cor>0.4），
corr.dat<-t(rbind(dat.gene,dat.gene))
#基因表达值的相关性分析，以Pearson相关系数为例
gene_cor <- rcorr(corr.dat)$r[1:nrow(dat.gene),(ncol(corr.dat)-length(rownames(dat.gene))+1):ncol(corr.dat)]
#将获得的相关性矩阵转换为两两对应的数据框结构
gene_cor <- reshape2::melt(gene_cor)
#gene_cor <- subset(gene_cor, value != 0)  #去除0值的相关性
head(gene_cor)  #前两列是两个基因名称，第三列为两个基因的相关性
colnames(gene_cor) <- c('symbol','symbol','correlation')
## 计算p值
gene_p <- rcorr(corr.dat)$P[1:nrow(dat.gene),(ncol(corr.dat)-length(rownames(dat.gene))+1):ncol(corr.dat)]
gene_p <- reshape2::melt(gene_p)
head(gene_p)  #前两列是两个基因名称，第三列为
gene_cor$pvalue <- gene_p$value
write.table(gene_cor,file = 'correlation.all.xls',sep = '\t',row.names = F,quote = F)
##|Cor| > 0.4, P< 0.001
cor.final <- subset(gene_cor,gene_cor$pvalu<0.001 & abs(gene_cor$correlation)>0.9)

## 62
write.table(cor.final,file = '04.correlation.final.xls',sep = '\t',row.names = F,quote = F)
