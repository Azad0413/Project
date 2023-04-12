rm(list = ls())
# 04 DEOIG----------
setwd("/data/nas1/luchunlin/project/BJTC-317")
if (! dir.exists("./03_DEMDG")){
  dir.create("./03_DEMDG")
}
setwd("./03_DEMDG")
library(readxl)
## wgcna 和 deg取交集------
sig_diff<-read.delim2("/data/nas1/luchunlin/project/BJTC-317/01_DEGs/DEG_sig.xls", row.names = 1)
modGenes<-read.delim2('/data/nas1/luchunlin/project/BJTC-317/02_WGCNA/modGene.xls')
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
DEmod<-modGenes[modGenes$modgene%in%rownames(sig_diff),]%>%as.data.frame()
## 444
write.table(DEmod,file = 'DEmod.xls',sep = '\t',quote = F,row.names = F)
#核线粒体基因取交集
mito<-read_csv(file = 'GeneCards-SearchResults.csv')
mito<-mito[which(mito$`Relevance score`>=4),]
##1870
demids<-DEmod[DEmod$.%in%mito$`Gene Symbol`,]%>%as.data.frame()
write.table(demids,file = 'DEMDG.xls',sep = '\t',row.names = F,quote = F)
mydata2<-list('DEmod'=DEmod$.,'MIDS'=mito$`Gene Symbol`)
ggvenn(mydata2,c('DEmod','MIDS'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')


