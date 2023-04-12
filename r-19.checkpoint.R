rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-386-10/")
if (! dir.exists("./19_checkpoint")){
  dir.create("./19_checkpoint")
}
setwd("./19_checkpoint")

hubgene <- read.delim2('../07_hubgene/hubgene.xls')
dat<-read.delim2('../00_rawdata/dat.all(fpkm).xls',row.names = 1)%>%lc.tableToNum()
colnames(dat)<-gsub('.','-',colnames(dat),fixed = T)
dat <- log2(dat+1)
hub_exp <- dat[hubgene$symbol,]

checkpoint <- data.frame(V1=c('PDCD1','CD274','PDCD1LG2','CTLA4','HAVCR2','LAG3'))
checkpoint_dat <- dat[checkpoint$V1,]

## 相关性----------
library(Hmisc)
gene.exp <- hub_exp

nc<-t(rbind(checkpoint_dat,gene.exp))
m=rcorr(nc)$r[1:nrow(checkpoint_dat),(ncol(nc)-length(hubgene$symbol)+1):ncol(nc)]
m<-t(m)
p=rcorr(nc)$P[1:nrow(checkpoint_dat),(ncol(nc)-length(hubgene$symbol)+1):ncol(nc)]
p<-t(p)
library(dplyr)
library(dplyr)

tmp = matrix(case_when(p<0.0001~"****",
                       p<0.0001~"***",
                       p<0.01~"**",
                       p<0.05~"*",
                       T~""),nrow = nrow(p))
cor <- m
cor <- signif(cor,3)
cor[abs(cor)<0.1] <- ''

textMatrix = paste(cor,"\n",
                   tmp, sep = "")
dim(textMatrix) = dim(m)
colnames(textMatrix)<-colnames(m)
rownames(textMatrix)<-rownames(m)
library(WGCNA)
pdf(file = '01.correlation.pdf',w=7,h=5.5)
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
               main = paste("lncRNA-checkpoint correlation"))
dev.off()
png(file = '01.correlation.png',w=500,h=400)
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
               main = paste("lncRNA-checkpoint correlation"))
dev.off()



##最负相关的

##HAVCR2-RP11-321G12.1  HAVCR2-RP11-195B3.1

cor_plot_dat <- data.frame(sample=colnames(gene.exp),t(checkpoint_dat[c('HAVCR2'),]),
                           t(gene.exp[c('RP11-321G12.1','RP11-195B3.1'),]))


dim(cor_plot_dat)
library(ggplot2)
library(ggpubr)
##HAVCR2-RP11-321G12.1 
colnames(cor_plot_dat)
library(ggstatsplot)
cor1 <- ggscatterstats(data = cor_plot_dat,
                       x = RP11.321G12.1,
                       y = HAVCR2,
                       centrality.para = "mean",
                       margins = "both",
                       xfill = "#A73030FF",
                       yfill = "#0073C2FF",
                       type = "pearson",
                       ylab = "HAVCR2",
                       xlab = "RP11-321G12.1",
                       marginal.type = "histogram",
                       title = "Relationship between RP11-321G12.1 and HAVCR2"
)+theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(size = 12))
cor1
ggsave(filename = "01.HAVCR2-RP11-321G12.1.png", height = 7, width = 8,cor1)
ggsave(filename = "01.HAVCR2-RP11-321G12.1.pdf", height = 7, width = 8,cor1)


##HAVCR2-RP11-195B3.1
colnames(cor_plot_dat)
library(ggstatsplot)
cor2 <- ggscatterstats(data = cor_plot_dat,
                       x = RP11.195B3.1,
                       y = HAVCR2,
                       centrality.para = "mean",
                       margins = "both",
                       xfill = "#A73030FF",
                       yfill = "#0073C2FF",
                       type = "pearson",
                       ylab = "HAVCR2",
                       xlab = "RP11-195B3.1",
                       marginal.type = "histogram",
                       title = "Relationship between RP11-195B3.1 and HAVCR2"
)+theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(size = 12))
cor2
ggsave(filename = "02.HAVCR2-RP11-195B3.1.png", height = 7, width = 8,cor2)
ggsave(filename = "02.HAVCR2-RP11-195B3.1.pdf", height = 7, width = 8,cor2)


