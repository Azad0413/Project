rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-386-10/")
if (! dir.exists("./07_hubgene")){
  dir.create("./07_hubgene")
}
setwd("./07_hubgene")
xgb <- read.delim2('../05_XGBoost/xgbgene.xls')
brt <- read.delim2('../06_Boruta/02.brtgene.xls')

hubgene <- data.frame(symbol=intersect(xgb$Feature,brt$symbol))
## 8
hubgene$symbol <- gsub('_','-',hubgene$symbol)
write.table(hubgene,file = 'hubgene.xls',sep = '\t',row.names = F,quote = F)
library(ggvenn)
mydata<-list('XGBoost'=xgb$Feature,'Boruta'=brt$symbol)
pdf(file = '01.hubgene.pdf',w=6,h=6)
ggvenn(mydata,c('XGBoost','Boruta'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()
png(file = '01.hubgene.png',w=400,h=400)
ggvenn(mydata,c('XGBoost','Boruta'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()
