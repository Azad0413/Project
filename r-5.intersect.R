rm(list = ls())
setwd("/data/nas1/luchunlin/project/SJZZK-428-10/")
if (! dir.exists("./05_intersect")){
  dir.create("./05_intersect")
}
setwd("./05_intersect")

# PART 1 GEO差异取交集-----
## 
diff1 <- read.delim2('../05_DEGs(GSE40595)/DEG_sig(GSE40595).xls')
up1 <- data.frame(symbol=rownames(diff1[which(diff1$change=='UP'),]))
down1 <- data.frame(symbol=rownames(diff1[which(diff1$change=='DOWN'),]))

diff2 <- read.delim2('../04_DEGs(GSE66957)/DEG_sig(GSE66957).xls')
up2 <- data.frame(symbol=rownames(diff2[which(diff2$change=='UP'),]))
down2 <- data.frame(symbol=rownames(diff2[which(diff2$change=='DOWN'),]))

## UP
up <- data.frame(symbol=intersect(up1$symbol,up2$symbol))
##3302
library(ggvenn)
mydata1<-list('UP(GSE40595)'=up1$symbol,'UP(GSE66957)'=up2$symbol)
pdf('01.venn(up).pdf',w=6,h=6)
ggvenn(mydata1,c('UP(GSE40595)','UP(GSE66957)'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()
png('01.venn(up).png',w=500,h=500)
ggvenn(mydata1,c('UP(GSE40595)','UP(GSE66957)'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()
down <- data.frame(symbol=intersect(down1$symbol,down2$symbol))
##1766
mydata2<-list('DOWN(GSE40595)'=down1$symbol,'DOWN(GSE66957)'=down2$symbol)
pdf('02.venn(down).pdf',w=6,h=6)
ggvenn(mydata2,c('DOWN(GSE40595)','DOWN(GSE66957)'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()
png('02.venn(down).png',w=500,h=500)
ggvenn(mydata2,c('DOWN(GSE40595)','DOWN(GSE66957)'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()

diff <- rbind(up,down)
##5068
write.table(diff,file = 'DEG(intersect).xls',sep = '\t',row.names = F,quote = F)

# PART 2 GEO-mRNAsi-氧化应激取交集-----
## DEG mRNAsi
diffsi <- read.delim2('../02_DEmRNAsi/DEG_sig(mRNAsi).xls',row.names = 1)
DESIG <- data.frame(symbol=intersect(rownames(diffsi),diff$symbol))
## 202
## 氧化应激
osgene <- read.csv('GeneCards-SearchResults.csv') 
osgene <- osgene[which(osgene$Relevance.score>2),]
##3840
intersect <- data.frame(symbol=intersect(osgene$Gene.Symbol,DESIG$symbol))
##43
write.table(intersect,file = 'intersect_final.xls',sep = '\t',row.names = F,quote = F)
mydata3<-list('DEG(OV)'=diff$symbol,'DEG(mRNAsi)'=rownames(diffsi),'OS'=osgene$Gene.Symbol)
pdf('03.intersect.pdf',w=5,h=5)
ggvenn(mydata3,c('DEG(OV)','DEG(mRNAsi)','OS'),
       fill_color = c("#ffb2b2","#b2e7cb","#F0E442"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb","#F0E442"),
       text_color = 'black')
dev.off()
png('03.intersect.png',w=400,h=400)
ggvenn(mydata3,c('DEG(OV)','DEG(mRNAsi)','OS'),
       fill_color = c("#ffb2b2","#b2e7cb","#F0E442"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb","#F0E442"),
       text_color = 'black')
dev.off()
