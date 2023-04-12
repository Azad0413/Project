rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/GY0324-12/")
if (! dir.exists("./00_rawdata")){
  dir.create("./00_rawdata")
}
setwd("./00_rawdata")

##转录组数据---------
library(readxl)
library(tidyverse)
library(lance)
exp <- read_xlsx('datall.xlsx')%>%
  select(-'ProbeName')%>%
  na.omit()
dat <- exp%>%lc.tableToNum()%>%
  mutate(rowMean=rowMeans(exp[2:16]))%>%
  arrange(desc(rowMean))%>%
  distinct(GeneSymbol,.keep_all = T)%>%
  select(-rowMean)%>%
  column_to_rownames(colnames(.)[1])

colnames(dat)
dat1 <- dat[,c(1:3,13:15)]
dat2 <- dat[,c(4:6,13:15)]
dat3 <- dat[,c(7:9,13:15)]
dat4 <- dat[,c(10:12,13:15)]



write.table(dat1,file = 'datA1.xls',sep = '\t',row.names = T,quote = F)
write.table(dat2,file = 'datA2.xls',sep = '\t',row.names = T,quote = F)
write.table(dat3,file = 'datA3.xls',sep = '\t',row.names = T,quote = F)
write.table(dat4,file = 'datB4.xls',sep = '\t',row.names = T,quote = F)
#write.table(dat5,file = 'datC5.xls',sep = '\t',row.names = T,quote = F)


###代谢组
meta <- read_xlsx('meta.xlsx')%>%select(-ID)
meta <- meta[!duplicated(meta$Peak),]%>%column_to_rownames(var = 'Peak')
write.table(meta,file = 'dat.meta.xls',sep = '\t',row.names = T,quote = F)

metabo <- t(meta)%>%as.data.frame()%>%rownames_to_column(var = 'sample')
metabo <- metabo[c(1:50),]
metabo$group <- c(rep('A1',10),rep('A2',10),rep('A3',10),rep('B4',10),rep('C5',10))
metabo <- metabo[,c(1,384,2:383)]
write.table(metabo,file = 'metabo.txt',sep = '\t',row.names = F,quote = F)
