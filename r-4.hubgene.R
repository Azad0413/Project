rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-399-11/")
if (! dir.exists("./04_hubgene")){
  dir.create("./04_hubgene")
}
setwd("./04_hubgene")

sig_diff <- read.delim2('../01_DEGs/DEG_sig.xls',row.names = 1)%>%lc.tableToNum()
modgene <- read.delim2('../03_WGCNA/modGene.xls')
hubgene <- data.frame(symbol=intersect(rownames(sig_diff),modgene$modgene))
write.table(hubgene,file = 'hubgene.xls',sep = '\t',row.names = F,quote = F)
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

