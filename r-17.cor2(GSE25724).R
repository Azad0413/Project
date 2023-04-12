rm(list = ls())
setwd("/data/nas1/luchunlin/project/LZZK-519-10/")
if (! dir.exists("./17_cor2(GSE25724)")){
  dir.create("./17_cor2(GSE25724)")
}
setwd("./17_cor2(GSE25724)")

##PART A 炎性因子-------
geneset <- c('TNF','IL1B','IFNG','VCAM1','ICAM1','CXCL10','IL2','IL4',
             'IL6','IL10','IL12A','IL12B','CD163','CCL5','NFKB1','TLR4')
dat <- read.delim2('../00_rawdata/dat(GSE25724).xls',row.names = 1)%>%lc.tableToNum()
hubgene <- read.delim2('../06_PPI/hubgene.xls')
hub.dat <- dat[hubgene$Symbol,]
cor.dat <- dat[geneset,]
cor.dat <- na.omit(cor.dat)
library(Hmisc)
nc<-t(rbind(cor.dat,hub.dat))
m=rcorr(nc,type = 'spearman')$r[1:nrow(cor.dat),(ncol(nc)-length(rownames(hub.dat))+1):ncol(nc)]
m<-t(m)

p=rcorr(nc,type = 'spearman')$P[1:nrow(cor.dat),(ncol(nc)-length(rownames(hub.dat))+1):ncol(nc)]
p<-t(p)

tmp = matrix(case_when(p<0.0001~"****",
                       p<0.0001~"***",
                       p<0.01~"**",
                       p<0.05~"*",
                       T~""),nrow = nrow(p))
cor <- m
cor <- round(cor,digits = 3)
cor[abs(cor)<0.55] <- ''

textMatrix = paste(cor,"\n",
                   tmp, sep = "")
dim(textMatrix) = dim(m)
colnames(textMatrix)<-colnames(m)
rownames(textMatrix)<-rownames(m)
write.table(textMatrix,file = 'correlaion.xls',sep = '\t',row.names = T,quote = F)

library(WGCNA)
pdf(file = '01.correlation.pdf',w=10,h=5)
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
               main = paste("hubgene-Inflammatory factor correlation"))
dev.off()
png(file = '01.correlation.png',w=700,h=350)
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
               main = paste("hubgene-Inflammatory factor correlation"))
dev.off()
