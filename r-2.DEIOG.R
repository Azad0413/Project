# DEIOG -----
rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-334/")
if (! dir.exists("./02_DEIOG")){
  dir.create("./02_DEIOG")
}
setwd("./02_DEIOG")
library(readxl)
library(tidyverse)
## 差异基因与免疫基因集、氧化应激基因集取交集
sig_diff<-read.delim2("/data/nas1/luchunlin/project/BJTC-334/01_DEGs/DEG_sig.xls", row.names = 1)
## 免疫基因集
immport<-read_xlsx('Immport_IRGs.xlsx')
innatedb<-read_xls('innatedb.xls')
immport<-data.frame(symbol=immport$Symbol)
innatedb<-data.frame(symbol=innatedb$`Gene Symbol`)
immune<-rbind(immport,innatedb)
immune<-immune[!duplicated(immune$symbol),]%>%as.data.frame()
## 2533

## 氧化应激基因集
oxstress<-read_xlsx('oxidative+stress.xlsx')

immu.oxstress<-immune[immune$.%in%oxstress$Symbol,]%>%as.data.frame()

DEIOG<-immu.oxstress[immu.oxstress$.%in%rownames(sig_diff),]%>%as.data.frame()
write.table(DEIOG,file = 'DEIOG.xls',row.names = F,quote = F,sep = '\t')

library(ggvenn)
mydata<-list('DEG'=rownames(sig_diff),'immune'=immune$.,'oxidative stress'=oxstress$Symbol)
ggvenn(mydata,c('DEG','immune','oxidative stress'),
       fill_color = c("#ffb2b2","#b2e7cb","#F0E442"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb","#F0E442"),
       text_color = 'black')
