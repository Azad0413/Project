rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/BJTC-420-1/")
if (! dir.exists("./13_drug")){
  dir.create("./13_drug")
}
setwd("./13_drug")
install.packages('enrichR')
library(enrichR)
dbs <- listEnrichrDbs()

##DsigDB
dbs$libraryName
dbs <- c("DSigDB")
hubgene <- read.delim('../05_PPI/hubgene.xls')
symbol <- hubgene$symbol

enrichr <- enrichr(symbol,dbs)
result <- data.frame(enrichr$DSigDB)
result <- result[which(result$Adjusted.P.value<0.05),]
write.table(result,file = 'result.xls',sep = '\t',row.names = F,quote = F)
