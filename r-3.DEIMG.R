rm(list = ls())
# 04 DEOIG----------
setwd("/data/nas1/luchunlin/project/BJTC-320")
if (! dir.exists("./03_DEIMG")){
  dir.create("./03_DEIMG")
}
setwd("./03_DEIMG")

## wgcna 和 deg取交集------
sig_diff<-read.delim2("/data/nas1/luchunlin/project/BJTC-320/01_DEGs/DEG_sig.xls", row.names = 1)
modGenes<-read.delim2('/data/nas1/luchunlin/project/BJTC-320/02_WGCNA/modGene.xls')
library(ggvenn)
mydata1<-list('DEGs'=rownames(sig_diff),'MOD gene'=modGenes$modgene)
ggvenn(mydata1,c('DEGs','MOD gene'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
DEmod<-modGenes[modGenes$modgene%in%rownames(sig_diff),]
## 761
write.table(DEmod,file = 'DEmod.xls',sep = '\t',quote = F,row.names = F)
## 氧化应激相关基因集、免疫相关基因集取交集
immport<-read_xlsx('Immport_IRGs.xlsx')
innatedb<-read_xls('innatedb.xls')
immport<-data.frame(symbol=immport$Symbol)
innatedb<-data.frame(symbol=innatedb$`Gene Symbol`)
immune<-rbind(immport,innatedb)
immune<-immune[!duplicated(immune$symbol),]%>%as.data.frame()
##2533
oxidative<-read_xlsx('oxidative+stress.xlsx')

mydata2<-list('DEMod gene'=DEmod,'immune'=immune$.,'oxidative stress'=oxidative$Symbol)
ggvenn(mydata2,c('DEMod gene','immune','oxidative stress'),
       fill_color = c("#ffb2b2","#b2e7cb","#F0E442"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb","#F0E442"),
       text_color = 'black')
DEIMG<-immune[immune$.%in%DEmod,]%>%as.data.frame()
DEOSG<-DEIMG[DEIMG$.%in%oxidative$Symbol,]%>%as.data.frame()

write.table(DEOSG,file = 'DEOSG.xls',sep = '\t',row.names = F,quote = F)
