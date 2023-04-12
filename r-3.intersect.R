rm(list = ls())
setwd("/data/nas1/luchunlin/project/CD-0601-2/")
if (! dir.exists("./03_intersect")){
  dir.create("./03_intersect")
}
setwd("./03_intersect")
diff <- read.delim2('../02_DEGs/DEG_sig(GSE65682).xls',row.names = 1)
modgene <- read.delim2('../01_WGCNA/modgene.xls')
intersect <- intersect(rownames(diff),modgene$symbol)%>%as.data.frame()
##719
geneset <- read_xlsx('geneset.xlsx')
geneset <- geneset[!duplicated(geneset$symbol),]
##471

intersect <- data.frame(symbol=intersect(intersect$.,geneset$symbol))
##21
write.table(intersect,file = 'intersect.xls',sep = '\t',row.names = F,quote = F)
library(ggvenn)
mydata<-list('DEGs'=rownames(diff),'modGenes'=modgene$symbol,'AAMRGs'=geneset$symbol)
pdf('01.venn.pdf',w=5,h=5)
ggvenn(mydata,c('DEGs','modGenes','AAMRGs'),
       fill_color = c("#ffb2b2","#b2e7cb","#F0E442"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       text_size = 4,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb","#F0E442"),
       text_color = 'black')
dev.off()
png('01.venn.png',w=400,h=400)
ggvenn(mydata,c('DEGs','modGenes','AAMRGs'),
       fill_color = c("#ffb2b2","#b2e7cb","#F0E442"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       text_size = 4,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb","#F0E442"),
       text_color = 'black')
dev.off()
