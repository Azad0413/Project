rm(list = ls())
setwd("/data/nas1/luchunlin/project/HZ0301-3/")
if (! dir.exists("./03_intersect")){
  dir.create("./03_intersect")
}
setwd("./03_intersect")

diff <- read.delim2('../01_DEGs/DEG_sig.xls',row.names = 1)
modgene <- read.delim2('../02_WGCNA/modGene.xls')
intersect <- intersect(rownames(diff),modgene$modgene)%>%as.data.frame()
##478
geneset1 <- read.csv('Table_1.csv')
##118
geneset2 <- read.csv('GeneCards-SearchResults.csv')
geneset2 <- geneset2[which(geneset2$Relevance.score>5),]
geneset2 <- data.frame(GeneID=geneset2$Gene.Symbol)

geneset <- rbind(geneset1,geneset2)
length(unique(geneset$GeneID))
geneset <- geneset[!duplicated(geneset$GeneID),]%>%as.data.frame()


intersect <- data.frame(symbol=intersect(intersect$.,geneset$.))
##33
write.table(intersect,file = 'intersect.xls',sep = '\t',row.names = F,quote = F)
library(ggvenn)
mydata<-list('DEGs'=rownames(diff),'modGene'=modgene$modgene,'Glutamine metabolism'=geneset$.)
pdf('01.venn.pdf',w=5,h=5)
ggvenn(mydata,c('DEGs','modGene','Glutamine metabolism'),
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
ggvenn(mydata,c('DEGs','modGene','Glutamine metabolism'),
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

