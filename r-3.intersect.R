rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-385-10/")
if (! dir.exists("./03_intersect")){
  dir.create("./03_intersect")
}
setwd("./03_intersect")
sig_diff <- read.delim2('../01_DEGs/DEG_sig.xls')
modgene <- read.delim2('../02_WGCNA/modgene.xls')
intersect <- intersect(rownames(sig_diff),modgene$modgene)%>%as.data.frame()
colnames(intersect) <- 'symbol'
write.table(intersect,file = 'intersect.xls',sep = '\t',row.names = F,quote = F)
library(ggvenn)
mydata<-list('DEGs'=rownames(sig_diff),'modGene'=modgene$modgene)
pdf('01.venn.pdf',w=5,h=5)
ggvenn(mydata,c('DEGs','modGene'),
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
ggvenn(mydata,c('DEGs','modGene'),
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
