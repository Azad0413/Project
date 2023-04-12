rm(list = ls())
setwd("/data/nas1/luchunlin/project/HF-0106-2/")
if (! dir.exists("./07_features")){
  dir.create("./07_features")
}
setwd("./07_features")
library(tidyverse)
library(lance)
lasso <- read.delim2('../04_lasso/lasso.gene.xls')
rf <- read.delim2('../05_randomforest/rf.gene.xls')

svmrfe <- read.delim2('../06_SVM_RFE/svmrfe_result.txt')
hubgene <- data.frame(symbol=intersect(lasso$x,svmrfe$x))

hubgene <- data.frame(symbol=intersect(hubgene$symbol,rf$symbol))

## 4
# hubgene$symbol <- gsub('_','-',hubgene$symbol)
write.table(hubgene,file = 'features.xls',sep = '\t',row.names = F,quote = F)
library(ggvenn)
mydata<-list('Lasso'=lasso$x,'RandomForest'=rf$symbol,'SVM-RFE'=svmrfe$x)
pdf('01.venn.pdf',w=5,h=5)
ggvenn(mydata,c('Lasso','RandomForest','SVM-RFE'),
       fill_color = c("#ffb2b2","#b2e7cb","#F0E442"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       text_size = 4,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb","#F0E442"),
       text_color = 'black')
dev.off()
png('01.venn.png',w=400,h=400)
ggvenn(mydata,c('Lasso','RandomForest','SVM-RFE'),
       fill_color = c("#ffb2b2","#b2e7cb","#F0E442"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       text_size = 4,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb","#F0E442"),
       text_color = 'black')
dev.off()
