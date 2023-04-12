rm(list = ls())
# 01 获取数据集--------------
setwd("/data/nas1/luchunlin/project/BJTC-356/")
if (! dir.exists("./00_rawdata")){
  dir.create("./00_rawdata")
}
setwd("./00_rawdata")
library(tidyverse)
library(lance)
##转录组
dat.all <- read.delim2('core_table.xls')
colnames(dat.all)
dat <- dat.all[,c(2,12:21)]%>%column_to_rownames(var = 'gene_id_alias')%>%lc.tableToNum()
dat <- na.omit(dat)
write.table(dat,file = 'dat.fpkm.xls',sep = '\t',row.names = T,quote = F)

##蛋白组

dat.protein.all <- read_xlsx('BPI_16763_renzuzhi_TMT_uniprot-human_9606_20180424.fasta_psms.xlsx')
colnames(dat.protein.all)
dat.protein <- dat.protein.all[,c(10,17:22)]
length(unique(dat.protein$Protein.Group.Accessions))
#colnames(dat.protein) <- c('ID','sample1_1','sample1_2','sample1_3','sample2_1','sample2_2','sample2_3')
dat.protein <- dat.protein[!duplicated(dat.protein$Protein.Group.Accessions)]
dat.protein <- na.omit(dat.protein)
dat.protein <- column_to_rownames(dat.protein,var = 'Protein.Group.Accessions')

write.table(dat.protein,file = 'protein.dat.xls',sep = '\t',row.names = T,quote = F)

