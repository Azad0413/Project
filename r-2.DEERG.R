rm(list = ls())
setwd("/data/nas1/luchunlin/project/TY0307-11")
if (! dir.exists("./02_DEERS")){
  dir.create("./02_DEERS")
}
setwd("./02_DEERS")
geneset <- read.csv('GeneCards-SearchResults.csv')

ERS <- geneset[which(geneset$Relevance.score>=10),]
##439
diff <- read.delim2('../01_DEGs/DEG_sig.xls')

DEERS <- data.frame(symbol=intersect(rownames(diff),ERS$Gene.Symbol))
## 13
write.table(DEERS,file = 'DEERS.xls',sep = '\t',row.names = F,quote = F)

library(ggvenn)
mydata<-list('DEGs'=rownames(diff),'ERS'=ERS$Gene.Symbol)
pdf('01.DEERS.pdf',w=6,h=6)
ggvenn(mydata,c('DEGs','ERS'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()
png('01.DEERS.png',w=500,h=500)
ggvenn(mydata,c('DEGs','ERS'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()

