rm(list = ls())
# 03 WGCNA----------
setwd("/data/nas1/luchunlin/project/BJTC-302")
if (! dir.exists("./05_PPI")){
  dir.create("./05_PPI")
}
setwd("./05_PPI")
degree<-read.csv('string_interactions_short.tsv_Degree_top10 default node.csv')
#dmnc<-read.csv('string_interactions_short.tsv_DMNC_top10 default node.csv')
mcc<-read.csv('string_interactions_short.tsv_MCC_top10 default node.csv')
mnc<-read.csv('string_interactions_short.tsv_MNC_top10 default node.csv')
hub<-degree[degree$name%in%mcc$name,]
hub<-hub[hub$name%in%mnc$name,]
# hub<-degree[degree$name%in%dmnc$name,]
# hub<-hub[hub$name%in%mcc$name,]
# hub<-hub[hub$name%in%mnc$name,]
library(ggvenn)
mydata<-list('Degree'=degree$name,'MCC'=mcc$name,'MNC'=mnc$name)
pdf('hub_venn.pdf',w=8,h=6)
ggvenn(mydata,c('Degree','MCC','MNC'),
       fill_color = c("#ffb2b2","#b2e7cb","#F0E442"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb","#F0E442"),
       text_color = 'black')
dev.off()
png('hub_venn.png',w=700,h=500)
ggvenn(mydata,c('Degree','DMNC','MCC','MNC'),
       fill_color = c("#ffb2b2","#b2e7cb","#F0E442", "#0072B2"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb","#F0E442", "#0072B2"),
       text_color = 'black')
dev.off()
hub<-data.frame(hubgene=hub$name)
write.table(hub,file = 'hubgene.xls',sep = '\t',row.names = F,quote = F)
