rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/XA-0214-1/")
if (! dir.exists("./02_intersect")){
  dir.create("./02_intersect")
}
setwd("./02_intersect")

geneset <- read_xlsx('angiogenesis.xlsx')

genecard <- read.csv('GeneCards-SearchResults.csv')
gene <- genecard[which(genecard$Relevance.score>5),]
##155
gene <- data.frame(Symbol=gene$Gene.Symbol)
gene.all <- rbind(geneset,gene)
gene.all <- gene.all[!duplicated(gene.all$Symbol),]

##同源基因转化
## 将基因转换为大鼠的同源基因
# install.packages('homologene')
library(homologene)
homologene::taxData
homologene <- homologene(gene.all$Symbol,inTax = 9606,outTax = 10116)
dim(homologene)
# 169   4
write.table(homologene,file = "homologene.xls",
            quote = F,
            sep = "\t",
            row.names = T)

diff <- read.delim2('../01_DEGs/DEG_sig(GSE97537).xls')%>%lc.tableToNum()

DEAOG <- data.frame(symbol=intersect(rownames(diff),homologene$`10116`))
## 38
write.table(DEAOG,file = 'DEAOG.xls',sep = '\t',row.names = F,quote = F)

library(ggvenn)
mydata<-list('DEGs'=rownames(diff),'Angiogenesis'=homologene$`10116`)
pdf('01.DEERS.pdf',w=6,h=6)
ggvenn(mydata,c('DEGs','Angiogenesis'),
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
ggvenn(mydata,c('DEGs','Angiogenesis'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()

