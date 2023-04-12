rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-406-12/")
if (! dir.exists("./13_TMB")){
  dir.create("./13_TMB") 
}
setwd("./13_TMB")
## PART A TMB------

maf <- TCGAmutations::tcga_load(study = "LUAD")
library(maftools)
riskscore <- read.table("../09_risk/risk.xls", header = T)
riskscore$risk <- ifelse(riskscore$riskScore>median(riskscore$riskScore), "High", "Low")
model_gene <- read.csv("../08_Lasso/lasso_genes.csv",header = F)
maf_sample <- data.frame(barcode = maf@clinical.data$Tumor_Sample_Barcode)
maf_sample$sample <- stringr::str_sub(maf_sample$barcode, 1, 16)
maf_sample <- merge(maf_sample, riskscore, by.x = "sample",by.y = 'id')
High.sample <- riskscore$id[which(riskscore$risk=='High')]
Low.sample <- riskscore$id[which(riskscore$risk=='Low')]
sample<-subset(maf_sample)$barcode
maf<-subsetMaf(maf,tsb = sample)
#maf<-subsetMaf(maf,tsb = sample,genes = model_gene$V1)

##HIGH LOW分开
high <- maf_sample[(maf_sample$sample)%in%High.sample,]
high <- subset(high)$barcode
maf.high <- subsetMaf(maf,tsb = high)
pdf(file = "01.oncoplot.high.pdf", height = 8, width = 10)
oncoplot(maf = maf.high, top = 20)
dev.off()

png(file = "01.oncoplot.high.png", family = "Times", height = 8, width = 10, units = "in", res = 600)
oncoplot(maf = maf.high, top = 20)
dev.off()

low <- maf_sample[(maf_sample$sample)%in%Low.sample,]
low <- subset(low)$barcode
maf.low <- subsetMaf(maf,tsb = low)
pdf(file = "02.oncoplot.low.pdf", height = 8, width = 10)
oncoplot(maf = maf.low, top = 20)
dev.off()

png(file = "02.oncoplot.low.png", family = "Times", height = 8, width = 10, units = "in", res = 600)
oncoplot(maf = maf.low, top = 20)
dev.off()


## 与显著突变基因的相关性------
dat <- read.delim2('../00_rawdata/dat.fpkm.xls',row.names = 1)%>%lc.tableToNum()
colnames(dat)<-gsub('.','-',colnames(dat),fixed = T)

### high---
high.gene <- maf.high@gene.summary$Hugo_Symbol[c(1:20)]
high.dat <- dat[high.gene,]
hub.dat <- dat[model_gene$V1,]
library(Hmisc)
gene.exp <- hub.dat
nc<-t(rbind(high.dat,gene.exp))
m=rcorr(nc,type = 'spearman')$r[1:nrow(high.dat),(ncol(nc)-length(rownames(gene.exp))+1):ncol(nc)]
m<-t(m)
p=rcorr(nc,type = 'spearman')$P[1:nrow(high.dat),(ncol(nc)-length(rownames(gene.exp))+1):ncol(nc)]
p<-t(p)

library(dplyr)
library(dplyr)

tmp = matrix(case_when(p<0.0001~"****",
                       p<0.0001~"***",
                       p<0.01~"**",
                       p<0.05~"*",
                       T~""),nrow = nrow(p))
cor <- m
cor <- round(cor,digits = 3)
#cor <- signif(cor,3)
cor[abs(cor)<0.083] <- ''
textMatrix = paste(cor,"\n",
                   tmp, sep = "")

dim(textMatrix) = dim(m)
colnames(textMatrix)<-colnames(m)
rownames(textMatrix)<-rownames(m)

library(WGCNA)
pdf(file = '03.correlation(high).pdf',w=12,h=5)
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
               main = paste("Model gene-TMB gene Correlation(High risk)"))
dev.off()
png(file = '03.correlation(high).png',w=800,h=400)
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
               main = paste("Model gene-TMB gene Correlation(High risk)"))
dev.off()


### low---
low.gene <- maf.low@gene.summary$Hugo_Symbol[c(1:20)]
low.dat <- dat[low.gene,]
hub.dat <- dat[model_gene$V1,]
library(Hmisc)
gene.exp <- hub.dat
nc<-t(rbind(low.dat,gene.exp))
m=rcorr(nc,type = 'spearman')$r[1:nrow(low.dat),(ncol(nc)-length(rownames(gene.exp))+1):ncol(nc)]
m<-t(m)
p=rcorr(nc,type = 'spearman')$P[1:nrow(low.dat),(ncol(nc)-length(rownames(gene.exp))+1):ncol(nc)]
p<-t(p)

library(dplyr)
library(dplyr)

tmp = matrix(case_when(p<0.0001~"****",
                       p<0.0001~"***",
                       p<0.01~"**",
                       p<0.05~"*",
                       T~""),nrow = nrow(p))
cor <- m
cor <- round(cor,digits = 3)
#cor <- signif(cor,3)
cor[abs(cor)<0.083] <- ''
textMatrix = paste(cor,"\n",
                   tmp, sep = "")

dim(textMatrix) = dim(m)
colnames(textMatrix)<-colnames(m)
rownames(textMatrix)<-rownames(m)


pdf(file = '04.correlation(low).pdf',w=12,h=5)
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
               main = paste("Model gene-TMB gene Correlation(Low risk)"))
dev.off()
png(file = '04.correlation(low).png',w=800,h=400)
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
               main = paste("Model gene-TMB gene Correlation(Low risk)"))
dev.off()

