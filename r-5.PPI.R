rm(list = ls())
# 02 差异分析----------
setwd("/data/nas1/luchunlin/project/BJTC-321")
if (! dir.exists("./05_PPI")){
  dir.create("./05_PPI")
}
setwd("./05_PPI")
Degree<-read.csv('string_interactions_short.tsv_Degree_top20 default node.csv')
EPC<-read.csv('string_interactions_short.tsv_EPC_top20 default node.csv')
MCC<-read.csv('string_interactions_short.tsv_MCC_top20 default node.csv')
MNC<-read.csv('string_interactions_short.tsv_MNC_top20 default node.csv')

hubgene<-Degree[Degree$name%in%EPC$name,]
hubgene<-hubgene[hubgene$name%in%MCC$name,]
hubgene<-hubgene[hubgene$name%in%MNC$name,]

hubgene<-data.frame(hubgene=hubgene$name)
write.table(hubgene,file = 'hubgene.xls',sep = '\t',row.names = F,quote = F)

mydata<-list('Degree'=Degree$name,'MCC'=MCC$name,'MNC'=MNC$name,'EPC'=EPC$name)
pdf('hub_venn.pdf',w=8,h=6)
ggvenn(mydata,c('Degree','MCC','MNC','EPC'),
       fill_color = c("#ffb2b2","#b2e7cb","#F0E442", "#0072B2"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb","#F0E442", "#0072B2"),
       text_color = 'black')
dev.off()
png('hub_venn.png',w=700,h=500)
ggvenn(mydata,c('Degree','MCC','MNC','EPC'),
       fill_color = c("#ffb2b2","#b2e7cb","#F0E442", "#0072B2"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb","#F0E442", "#0072B2"),
       text_color = 'black')
dev.off()
