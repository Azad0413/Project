rm(list = ls())
setwd("/data/nas1/luchunlin/project/JNZK-214-8/")
if (! dir.exists("./06_DEIRG")){
  dir.create("./06_DEIRG")
}
setwd("./06_DEIRG")

scDEGs <- read.delim2('../04_scRNA/06_DEGs/scDEGs.xls')

#scDEGs <- scDEGs[which(scDEGs$avg_log2FC>1),]
DEGs <- read.delim2('../05_DEGs/DEG_sig.xls',row.names = 1)
modGene <- read.delim2('../03_WGCNA/modgene.xls')
intersect <- intersect(scDEGs$symbol,rownames(DEGs))%>%as.data.frame()
intersect <- intersect(intersect$.,modGene$modgene)%>%as.data.frame()
colnames(intersect) <- 'symbol'
write.table(intersect,file = 'DEIRG.xls',sep = '\t',row.names = F,quote = F)
library(ggvenn)
mydata1<-list('DEG(scRNA)'=scDEGs$symbol,'DEG(TCGA)'=rownames(DEGs),'modGene'=modGene$modgene)
pdf('01.venn.pdf',w=5,h=5)
ggvenn(mydata1,c('DEG(scRNA)','DEG(TCGA)','modGene'),
       fill_color = c("#ffb2b2","#b2e7cb","#F0E442"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb","#F0E442"),
       text_color = 'black')
dev.off()
png('01.venn.png',w=400,h=400)
ggvenn(mydata1,c('DEG(scRNA)','DEG(TCGA)','modGene'),
       fill_color = c("#ffb2b2","#b2e7cb","#F0E442"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb","#F0E442"),
       text_color = 'black')
dev.off()
