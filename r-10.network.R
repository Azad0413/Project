rm(list = ls())
setwd("/data/nas1/luchunlin/project/CD-0601-2/")
if (! dir.exists("./10_network")){
  dir.create("./10_network")
}
setwd("./10_network")

hubgene <- read.delim2('../07_features/features.xls')

library(multiMiR)
gene2mir <- get_multimir(org = 'hsa',target = hubgene$symbol,summary = TRUE)

table(gene2mir@data$mature_mirna_id)

result <- gene2mir@data



gene2mir <- result[,c(3,4)]
colnames(gene2mir) <- c('miRNA','mRNA')

# gene2mir <- gene2mir[duplicated(gene2mir$miRNA),]
table(gene2mir$miRNA)

write.table(gene2mir,file = 'gene2mir.xls',sep = '\t',row.names = F,quote = F)

gene2tf <- read.csv('gene2tf.csv')

gene2tf <- gene2tf[,c(1,3)]

colnames(gene2tf) <- c('TF','mRNA')

network <- merge(gene2tf,gene2mir,by='mRNA')

network$link <- paste0(network$mRNA,network$TF,network$miRNA)
network <- network[!duplicated(network$link),]%>%dplyr::select(-'link')

write.table(network,file = 'network.xls',sep = '\t',row.names = F,quote = F)
