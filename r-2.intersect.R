rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-441-3/")
if (! dir.exists("./02_intersect")){
  dir.create("./02_intersect")
}
setwd("./02_intersect")
diff <- read.csv('../01_DEGs/DEG_sig.xls',sep = '\t',row.names = 1)
library(Ipaper)
smgs <- read_xlsx('SMGs.xlsx')
length(unique(smgs$symbol))
intersect <- data.frame(symbol=intersect(rownames(diff),smgs$symbol)%>%as.data.frame()) 
##13
write.table(intersect,file = 'DESMGs.xls',sep = '\t',row.names = F,quote = F)
library(ggvenn)
mydata<-list('DEGs'=rownames(diff),'SMGs'=smgs$symbol)
pdf('01.venn.pdf',w=5,h=5)
ggvenn(mydata,c('DEGs','SMGs'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       text_size = 4,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()
png('01.venn.png',w=400,h=400)
ggvenn(mydata,c('DEGs','SMGs'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       text_size = 4,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()
