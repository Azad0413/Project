rm(list = ls())
setwd("/data/nas1/luchunlin/project/CD-0601-2/")
if (! dir.exists("./07_features")){
  dir.create("./07_features")
}
setwd("./07_features")
library(tidyverse)
library(lance)
lasso <- read.delim2('../05_lasso/lasso.gene.xls')%>%lc.tableToNum()
svmrfe <- read.delim2('../06_SVM_RFE/svmrfe_result.txt')
hubgene <- data.frame(symbol=intersect(lasso$x,svmrfe$x))
## 7
# hubgene$symbol <- gsub('_','-',hubgene$symbol)
write.table(hubgene,file = 'features.xls',sep = '\t',row.names = F,quote = F)
library(ggvenn)
mydata<-list('Lasso'=lasso$x,'SVM-RFE'=svmrfe$x)
pdf(file = '01.venn.pdf',w=6,h=6)
ggvenn(mydata,c('Lasso','SVM-RFE'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()
png(file = '01.venn.png',w=400,h=400)
ggvenn(mydata,c('Lasso','SVM-RFE'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()
