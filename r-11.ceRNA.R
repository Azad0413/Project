rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/XA-0214-1/")
if (! dir.exists("./11_ceRNA")){
  dir.create("./11_ceRNA")
}
setwd("./11_ceRNA")

library(multiMiR)
mir1 <- read.csv('Stat3.csv')
mir2 <- read.csv('Col18a1.csv')
mir3 <- read.csv('Hmox1.csv')
mir4 <- read.csv('Ptgs2.csv')
mir5 <- read.csv('Egfr.csv')
length(unique(mir1$mirnaid))
##114
length(unique(mir2$mirnaid))
##475
length(unique(mir3$mirnaid))
##363
length(unique(mir4$mirnaid))
##390
mir <- rbind(mir1,mir2,mir3,mir4,mir5)
mir <- mir[!duplicated(mir$mirnaid),]

library(multiMiR)
mir2gene <- get_multimir(org = 'rat',mirna = mir$mirnaid,summary = TRUE)
table(mir2gene@data$target_symbol)

result <- mir2gene@data

library("rtracklayer")
#gtf_data = import('/data/nas1/luchunlin/pipeline/GENEANNO/gencode.v22.annotation.gtf.gz') #gtf的路径
gtf_data = import('Rattus_norvegicus.Rnor_6.0.96.gtf.gz') #gtf的路径
gtf_data = as.data.frame(gtf_data)
head(gtf_data)
table(gtf_data$gene_biotype)
lncRNA=gtf_data%>%
  dplyr::filter(gene_biotype=="lincRNA")%>%
  dplyr::select(gene_id,gene_biotype,gene_name)

result <- result[result$target_ensembl%in%lncRNA$gene_id,]


