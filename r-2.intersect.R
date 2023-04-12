rm(list = ls())
setwd("/data/nas1/luchunlin/project/HF-0106-2/")
if (! dir.exists("./02_intersect")){
  dir.create("./02_intersect")
}
setwd("./02_intersect")

diff <- read.delim2('../01_DEGs/DEG_sig(GSE169568).xls')
geneset <- read.csv('GeneCards-SearchResults.csv')
geneset <- geneset[which(geneset$Relevance.score>5),]

DECFRGs <- data.frame(symbol=intersect(rownames(diff),geneset$Gene.Symbol))
##10
write.table(DECFRGs,file = 'DECFRGs.xls',sep = '\t',row.names = F,quote = F)
BiocManager::install('ggvenn')
library(ggvenn)
mydata1<-list('DEGs'=rownames(diff),'CFRGs'=geneset$Gene.Symbol)
pdf(file = '01.venn.pdf',w=6,h=6)
ggvenn(mydata1,c('DEGs','CFRGs'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()
png(file = '01.venn.png',w=450,h=450)
ggvenn(mydata1,c('DEGs','CFRGs'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()

