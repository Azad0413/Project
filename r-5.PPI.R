rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/BJTC-420-1/")
if (! dir.exists("./05_PPI")){
  dir.create("./05_PPI")
}
setwd("./05_PPI")

hubba <- read.csv('HubbaTable.csv')
colnames(hubba)

degree <- hubba[,c('Symbol','Degree')]
closeness <- hubba[,c('Symbol','Closeness')]
betweenness <- hubba[,c('Symbol','Betweenness')]

degree <- degree[order(degree$Degree,decreasing = T),]
degree <- degree[which(degree$Degree>=11),]
write.table(degree,file = 'topdegree.xls',sep = '\t',row.names = F,quote = F)

closeness <- closeness[order(closeness$Closeness,decreasing = T),]
closeness <- closeness[c(1:21),]
write.table(closeness,file = 'topdcloseness.xls',sep = '\t',row.names = F,quote = F)

betweenness <- betweenness[order(betweenness$Betweenness,decreasing = T),]
betweenness <- betweenness[c(1:21),]
write.table(betweenness,file = 'topdbetweenness.xls',sep = '\t',row.names = F,quote = F)


hubgene <- data.frame(symbol=intersect(degree$Symbol,closeness$Symbol))
hubgene <- data.frame(symbol=intersect(hubgene$symbol,betweenness$Symbol))
##9

write.table(hubgene,file = 'hubgene.xls',sep = '\t',row.names = F,quote = F)

library(ggvenn)
mydata<-list('Degree'=degree$Symbol,'Closeness'=closeness$Symbol,'Betweenness'=betweenness$Symbol)
pdf('01.venn.pdf',w=5,h=5)
ggvenn(mydata,c('Degree','Closeness','Betweenness'),
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
ggvenn(mydata,c('Degree','Closeness','Betweenness'),
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
