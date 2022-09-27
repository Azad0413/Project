rm(list = ls())
setwd("/data/nas1/luchunlin/project/JNZK-218-8/")
if (! dir.exists("./04_hubgene")){
  dir.create("./04_hubgene")
}
setwd("./04_hubgene")

POAG_MOD<-read.delim2('../02_WGCNA(POAG)/modGene(POAG).xls')
AS_MOD<-read.delim2('../03_WGCNA(AS)/modGene(AS).xls')
POAG_AS_MOD<-POAG_MOD[POAG_MOD$modgene%in%AS_MOD$modgene,]%>%as.data.frame()
### 865
write.table(POAG_AS_MOD,file = 'POAG_AS_MODgene.xls',sep = '\t',row.names = F,quote = F)

library(ggvenn)
mydata1<-list('POAG modgene'=POAG_MOD$modgene,'AS modgene'=AS_MOD$modgene)
pdf('01.Mod_venn.pdf',w=5,h=5)
ggvenn(mydata1,c('POAG modgene','AS modgene'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()

png('01.Mod_venn.png',w=400,h=400)
ggvenn(mydata1,c('POAG modgene','AS modgene'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()


DEGs<-read.delim2('../01_DEGs/interDEGs.xls')

hubgene<-DEGs[DEGs$.%in%POAG_AS_MOD$.,]%>%as.data.frame()
colnames(hubgene)<-'symbol'
### 93
write.table(hubgene,file = 'hubgene.xls',sep = '\t',row.names = F,quote = F)
library(ggvenn)
mydata2<-list('DEGs'=DEGs$.,'MOD gene'=POAG_AS_MOD$.)
pdf('02.hub_venn.pdf',w=5,h=5)
ggvenn(mydata2,c('DEGs','MOD gene'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()
png('02.hub_venn.png',w=400,h=400)
ggvenn(mydata2,c('DEGs','MOD gene'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()
