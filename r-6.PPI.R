rm(list = ls())
setwd("/data/nas1/luchunlin/project/NN-0118-2/")
if (! dir.exists("./06_PPI")){
  dir.create("./06_PPI")
}
setwd("./06_PPI")

hubgene <- data.frame(symbol=c('MKI67','NDC80','ANLN','CENPF','CCNB1','MELK','AURKA','KIF4A',
                               'CCNB2','KIF20A','PRC1','BUB1B','BUB1','TPX2','CCNA2','CEP55',
                               'PBK','NUF2','NUSAP1','ASPM','RRM2','CENPE','TOP2A','KIF11'))
write.table(hubgene,file = 'hubgene.xls',sep = '\t',row.names = F,quote = F)
