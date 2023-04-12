rm(list = ls())
setwd("/data/nas1/luchunlin/project/TY0307-11/")
if (! dir.exists("./10_immucorrelation")){
  dir.create("./10_immucorrelation")
}
setwd("./10_immucorrelation")

##PART A 趋化因子--------
chemokine <- c('CCL27','CXCL5','CCL16','CCL7','CCL19','CCL11','XCL1','CXCL11',
               'CCL21','CCL13','CCL17','CCL5','CXCL1','CCL20','CCL23','CCL8',
               'CCL18','CXCL8','CCL25','CCL22','CXCL10','CCL1','CCL4','CXCL2',
               'CXCL9','CXCL13','CCL2','CXCL3','CXCL6','CXCL12','CX3CL1')

dat <- read.delim2('../00_rawdata/dat(GSE113079).xls',row.names = 1)%>%lc.tableToNum()
hubgene <- read.delim2('../04_model/hubgene.xls')
hub.dat <- dat[hubgene$x,]
cor.dat <- dat[chemokine,]
cor.dat <- na.omit(cor.dat)
library(Hmisc)
nc<-t(rbind(cor.dat,hub.dat))
m=rcorr(nc)$r[1:nrow(cor.dat),(ncol(nc)-length(rownames(hub.dat))+1):ncol(nc)]
m<-t(m)

p=rcorr(nc)$P[1:nrow(cor.dat),(ncol(nc)-length(rownames(hub.dat))+1):ncol(nc)]
p<-t(p)


tmp = matrix(case_when(p<0.0001~"****",
                       p<0.0001~"***",
                       p<0.01~"**",
                       p<0.05~"*",
                       T~""),nrow = nrow(p))
cor <- m
cor <- round(cor,digits = 3)
cor[abs(cor)<0.15] <- ''

textMatrix = paste(cor,"\n",
                   tmp, sep = "")
dim(textMatrix) = dim(m)
colnames(textMatrix)<-colnames(m)
rownames(textMatrix)<-rownames(m)
write.table(textMatrix,file = 'correlaion(chemokine).xls',sep = '\t',row.names = T,quote = F)

library(WGCNA)
pdf(file = '01.correlation(chemokine).pdf',w=15,h=7)
par(mar = c(9, 8, 3, 1));
labeledHeatmap(Matrix = m, 
               xLabels = colnames(m), 
               yLabels = rownames(m), 
               cex.lab = 1, 
               # ySymbols = colnames(MEs_col)[-20], 
               colorLabels = FALSE, 
               colors = blueWhiteRed(50, 0.8),
               textMatrix = textMatrix, 
               setStdMargins = FALSE, 
               cex.text = 0.8, 
               zlim = c(-1,1),
               main = paste("hubgene-chemokine correlation"))
dev.off()
png(file = '01.correlation(chemokine).png',w=1000,h=450)
par(mar = c(9, 8, 3, 1));
labeledHeatmap(Matrix = m, 
               xLabels = colnames(m), 
               yLabels = rownames(m), 
               cex.lab = 1, 
               # ySymbols = colnames(MEs_col)[-20], 
               colorLabels = FALSE, 
               colors = blueWhiteRed(50, 0.8),
               textMatrix = textMatrix, 
               setStdMargins = FALSE, 
               cex.text = 0.8, 
               zlim = c(-1,1),
               main = paste("hubgene-chemokine correlation"))
dev.off()


##PART B 免疫抑制剂--------
inhibitor <- c('IDO1','TGFB1','PDCD1','LGALS9','KDR','LAG3','IL10','CD96','TGFBR1','CD160','IL10RB')

cor.dat <- dat[inhibitor,]
cor.dat <- na.omit(cor.dat)
library(Hmisc)
nc<-t(rbind(cor.dat,hub.dat))
m=rcorr(nc)$r[1:nrow(cor.dat),(ncol(nc)-length(rownames(hub.dat))+1):ncol(nc)]
m<-t(m)

p=rcorr(nc)$P[1:nrow(cor.dat),(ncol(nc)-length(rownames(hub.dat))+1):ncol(nc)]
p<-t(p)


tmp = matrix(case_when(p<0.0001~"****",
                       p<0.0001~"***",
                       p<0.01~"**",
                       p<0.05~"*",
                       T~""),nrow = nrow(p))
cor <- m
cor <- round(cor,digits = 3)
cor[abs(cor)<0.15] <- ''

textMatrix = paste(cor,"\n",
                   tmp, sep = "")
dim(textMatrix) = dim(m)
colnames(textMatrix)<-colnames(m)
rownames(textMatrix)<-rownames(m)
write.table(textMatrix,file = 'correlaion(Immunoinhibitor).xls',sep = '\t',row.names = T,quote = F)


library(WGCNA)
pdf(file = '02.correlation(Immunoinhibitor).pdf',w=8,h=6)
par(mar = c(9, 8, 3, 1));
labeledHeatmap(Matrix = m, 
               xLabels = colnames(m), 
               yLabels = rownames(m), 
               cex.lab = 1, 
               # ySymbols = colnames(MEs_col)[-20], 
               colorLabels = FALSE, 
               colors = blueWhiteRed(50, 0.8),
               textMatrix = textMatrix, 
               setStdMargins = FALSE, 
               cex.text = 0.8, 
               zlim = c(-1,1),
               main = paste("hubgene-Immunoinhibitor correlation"))
dev.off()
png(file = '02.correlation(Immunoinhibitor).png',w=600,h=400)
par(mar = c(9, 8, 3, 1));
labeledHeatmap(Matrix = m, 
               xLabels = colnames(m), 
               yLabels = rownames(m), 
               cex.lab = 1, 
               # ySymbols = colnames(MEs_col)[-20], 
               colorLabels = FALSE, 
               colors = blueWhiteRed(50, 0.8),
               textMatrix = textMatrix, 
               setStdMargins = FALSE, 
               cex.text = 0.8, 
               zlim = c(-1,1),
               main = paste("hubgene-Immunoinhibitor correlation"))
dev.off()

##PART C 免疫激活剂--------
stimulator <- c('TNFRSF9','IL6','TNFRSF17','CD40','TNFRSF25','TNFRSF13B','ICOSLG',
                'CD48','CD70','TNFSF14','MICB','IL2RA','PVR','TNFRSF8','TNFSF4',
                'TNFSF9','ICOS','NT5E','CD86','CD28','LTA','CD27','CXCR4','TNFRSF14','CD80','CD40LG','ENTPD1')

cor.dat <- dat[stimulator,]
cor.dat <- na.omit(cor.dat)
library(Hmisc)
nc<-t(rbind(cor.dat,hub.dat))
m=rcorr(nc)$r[1:nrow(cor.dat),(ncol(nc)-length(rownames(hub.dat))+1):ncol(nc)]
m<-t(m)

p=rcorr(nc)$P[1:nrow(cor.dat),(ncol(nc)-length(rownames(hub.dat))+1):ncol(nc)]
p<-t(p)


tmp = matrix(case_when(p<0.0001~"****",
                       p<0.0001~"***",
                       p<0.01~"**",
                       p<0.05~"*",
                       T~""),nrow = nrow(p))
cor <- m
cor <- round(cor,digits = 3)
cor[abs(cor)<0.15] <- ''

textMatrix = paste(cor,"\n",
                   tmp, sep = "")
dim(textMatrix) = dim(m)
colnames(textMatrix)<-colnames(m)
rownames(textMatrix)<-rownames(m)
write.table(textMatrix,file = 'correlaion(Immunostimulator).xls',sep = '\t',row.names = T,quote = F)

library(WGCNA)
pdf(file = '03.correlation(Immunostimulator).pdf',w=14,h=7)
par(mar = c(9, 8, 3, 1));
labeledHeatmap(Matrix = m, 
               xLabels = colnames(m), 
               yLabels = rownames(m), 
               cex.lab = 1, 
               # ySymbols = colnames(MEs_col)[-20], 
               colorLabels = FALSE, 
               colors = blueWhiteRed(50, 0.8),
               textMatrix = textMatrix, 
               setStdMargins = FALSE, 
               cex.text = 0.8, 
               zlim = c(-1,1),
               main = paste("hubgene-Immunostimulator correlation"))
dev.off()
png(file = '03.correlation(Immunostimulator).png',w=900,h=450)
par(mar = c(9, 8, 3, 1));
labeledHeatmap(Matrix = m, 
               xLabels = colnames(m), 
               yLabels = rownames(m), 
               cex.lab = 1, 
               # ySymbols = colnames(MEs_col)[-20], 
               colorLabels = FALSE, 
               colors = blueWhiteRed(50, 0.8),
               textMatrix = textMatrix, 
               setStdMargins = FALSE, 
               cex.text = 0.8, 
               zlim = c(-1,1),
               main = paste("hubgene-Immunostimulator correlation"))
dev.off()

