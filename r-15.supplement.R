rm(list = ls())
setwd("/data/nas1/luchunlin/project/JNZK-214-8/")
if (! dir.exists("./16_supplement")){
  dir.create("./16_supplement")
}
setwd("./16_supplement")

hubgene <- read.delim2('../08_Lasso/lasso_genes.csv',header = F)
library(lance)
library(tidyverse)

## part 1 与免疫因子的相关性---------
##PART A 趋化因子--------
chemokine <- readxl::read_xlsx('immugeneset.xlsx')

dat <- read.delim2('../00_rawdata/dat.fpkm.xls',row.names = 1)%>%lc.tableToNum()
colnames(dat)<-gsub('.','-',colnames(dat),fixed = T)
group<-read.delim2("../01_CIBERSORT/group.xls")
tumor.sample<-group$sample[which(group$group=='Tumor')]
hub.dat <- dat[hubgene$V1,]
cor.dat <- dat[chemokine$symbol,]
cor.dat <- na.omit(cor.dat)
library(Hmisc)
nc<-t(rbind(cor.dat,hub.dat))
m=rcorr(nc,type = 'spearman')$r[1:nrow(cor.dat),(ncol(nc)-length(rownames(hub.dat))+1):ncol(nc)]
m<-t(m)

p=rcorr(nc,type = 'spearman')$P[1:nrow(cor.dat),(ncol(nc)-length(rownames(hub.dat))+1):ncol(nc)]
p<-t(p)
class(p)
library(dplyr)
# tmp = matrix(case_when(p<0.0001~"****",
#                        p<0.001~"***",
#                        p<0.01~"**",
#                        p<0.05~"*",
#                        T~""),nrow = nrow(p))
tmp <- p
tmp[tmp<0.0001] <- '****'
tmp[tmp<0.001] <- '***'
tmp[tmp<0.01] <- '**'
tmp[tmp<0.05] <- '*'
tmp[tmp>0.05] <- 'ns'


cor <- m
cor <- round(cor,digits = 2)
# cor[abs(cor)<0.15] <- ''

textMatrix = paste(cor,"\n",
                   tmp, sep = "")
dim(textMatrix) = dim(m)
colnames(textMatrix)<-colnames(m)
rownames(textMatrix)<-rownames(m)
write.table(textMatrix,file = 'correlaion(chemokine).xls',sep = '\t',row.names = T,quote = F)

library(WGCNA)
pdf(file = '01.correlation(chemokine).pdf',w=16,h=7)
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
png(file = '01.correlation(chemokine).png',w=1100,h=450)
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

# part2 TMB MSI----
###TMB-------
## PART A TMB------
library(TCGAmutations)
maf <- TCGAmutations::tcga_load(study = "COAD")
library(maftools)

maf_sample <- data.frame(barcode = maf@clinical.data$Tumor_Sample_Barcode)
maf_sample$sample <- stringr::str_sub(maf_sample$barcode, 1, 16)
sample<-subset(maf_sample)$barcode
maf<-subsetMaf(maf,tsb = sample)
#maf<-subsetMaf(maf,tsb = sample,genes = model_gene$V1)

pdf(file = "oncoplot.pdf", height = 8, width = 10)
oncoplot(maf = maf, top = 20)
dev.off()

png(file = "oncoplot.png", family = "Times", height = 8, width = 10, units = "in", res = 600)
oncoplot(maf = maf, top = 20)
dev.off()


##预后基因突变率--------
maf.mod<-subsetMaf(maf,tsb = sample,genes = hubgene$V1)

pdf(file = "mut_modgene.pdf", family = "Times", height = 5, width = 6)
plotmafSummary(maf = maf.mod, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

png(file = "mut_modgene.png", family = "Times", height = 5, width = 6, units = "in", res = 600)
plotmafSummary(maf = maf.mod, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

###MSI---------
setwd("/data/nas1/luchunlin/project/JNZK-214-8/16_supplement/")
if (! dir.exists("./MSI")){
  dir.create("./MSI")
}
setwd("./MSI")

msi <- readxl::read_xlsx('../msi.xlsx')

dat.hub <- t(hub.dat)%>%as.data.frame()
dat.hub <- log2(dat.hub+1)
dat.hub$sample <- substr(rownames(dat.hub),1,15)
dat.hub$rownames <- rownames(dat.hub)
dat.cor <- merge(msi,dat.hub,by='sample')

colnames(dat.cor)
dat.gene <- dat.cor[,c(3:9)]%>%column_to_rownames(var = 'rownames')
dat.msi <- dat.cor[,c(2,9)]%>%column_to_rownames(var = 'rownames')
dat.gene <- t(dat.gene)%>%as.data.frame()

library(ggplot2)
library(ggpubr)
library(ggExtra)
cnt <- 1
while (cnt<7) {
  x=as.numeric(dat.gene[cnt,])
  y=dat.msi
  
  library(psych)
  d <- psych::corr.test(y,x,use="complete",method = 'spearman')
  r <- data.frame(d$r)
  p <- data.frame(d$p)
  # write.table(p,file=paste(rownames(dat.hub)[cnt], "p.xls",sep='_'),sep="\t",quote=F)
  #  write.table(r,file=paste(rownames(dat.hub)[cnt], "r.xls",sep='_'),sep="\t",quote=F)
  correlation<-data.frame(rownames(p),r$d.r,p$d.p)
  colnames(correlation) <- c('gene','Correlation','p.value')
  write.table(correlation,file=paste(rownames(dat.gene)[cnt], "correlaion.xls",sep='_'),sep="\t",quote=F,row.names = F)
  
  for (i in 1:length(colnames(y))) {
    marker <- as.numeric(t(y)[i,])
    hubgene <- as.numeric(dat.gene[cnt,])
    df1 <- cbind(hubgene,marker)%>%as.data.frame()
    p1 <- ggplot(df1,aes(marker,hubgene))+
      xlab(colnames(dat.msi[i]))+ylab(rownames(dat.gene[cnt,]))+
      geom_point()+geom_smooth(method = 'lm',formula = y~x)+
      theme_bw()+
      stat_cor(method = 'spearman',aes(x=marker,y=hubgene))
    p1  
    p2 <- ggMarginal(p1,type = 'density',xparams = list(fill='orange'),
                     yparams = list(fill='blue'))
    
    p2
    ggsave(p2,filename = paste0(cnt,'.',rownames(dat.gene[cnt,]),'-',colnames(dat.msi[i]),'.pdf'),w=5,h=5)
    ggsave(p2,filename = paste0(cnt,'.',rownames(dat.gene[cnt,]),'-',colnames(dat.msi[i]),'.png'),w=5,h=5)
  }  
  cnt <- cnt+1
}




##NEO--------
setwd("/data/nas1/luchunlin/project/JNZK-214-8/16_supplement/")
if (! dir.exists("./NEO")){
  dir.create("./NEO")
}
setwd("./NEO")

NEO <- readxl::read_xlsx('../NEO.xlsx')

dat.hub <- t(hub.dat)%>%as.data.frame()
dat.hub <- log2(dat.hub+1)
dat.hub$sample <- substr(rownames(dat.hub),1,15)
dat.hub$rownames <- rownames(dat.hub)
colnames(NEO) <- c('sample','Neoantigen')
dat.cor <- merge(NEO,dat.hub,by='sample')

colnames(dat.cor)
dat.gene <- dat.cor[,c(3:9)]%>%column_to_rownames(var = 'rownames')
dat.NEO <- dat.cor[,c(2,9)]%>%column_to_rownames(var = 'rownames')
dat.gene <- t(dat.gene)%>%as.data.frame()

library(ggplot2)
library(ggpubr)
library(ggExtra)
cnt <- 1
while (cnt<7) {
  x=as.numeric(dat.gene[cnt,])
  y=dat.NEO
  
  library(psych)
  d <- psych::corr.test(y,x,use="complete",method = 'spearman')
  r <- data.frame(d$r)
  p <- data.frame(d$p)
  # write.table(p,file=paste(rownames(dat.hub)[cnt], "p.xls",sep='_'),sep="\t",quote=F)
  #  write.table(r,file=paste(rownames(dat.hub)[cnt], "r.xls",sep='_'),sep="\t",quote=F)
  correlation<-data.frame(rownames(p),r$d.r,p$d.p)
  colnames(correlation) <- c('gene','Correlation','p.value')
  write.table(correlation,file=paste(rownames(dat.gene)[cnt], "correlaion.xls",sep='_'),sep="\t",quote=F,row.names = F)
  
  for (i in 1:length(colnames(y))) {
    marker <- as.numeric(t(y)[i,])
    hubgene <- as.numeric(dat.gene[cnt,])
    df1 <- cbind(hubgene,marker)%>%as.data.frame()
    p1 <- ggplot(df1,aes(marker,hubgene))+
      xlab(colnames(dat.NEO[i]))+ylab(rownames(dat.gene[cnt,]))+
      geom_point()+geom_smooth(method = 'lm',formula = y~x)+
      theme_bw()+
      stat_cor(method = 'spearman',aes(x=marker,y=hubgene))
    p1  
    p2 <- ggMarginal(p1,type = 'density',xparams = list(fill='orange'),
                     yparams = list(fill='blue'))
    
    p2
    ggsave(p2,filename = paste0(cnt,'.',rownames(dat.gene[cnt,]),'-',colnames(dat.NEO[i]),'.pdf'),w=5,h=5)
    ggsave(p2,filename = paste0(cnt,'.',rownames(dat.gene[cnt,]),'-',colnames(dat.NEO[i]),'.png'),w=5,h=5)
  }  
  cnt <- cnt+1
}


#part3 与差异免疫细胞的相关性--------
## spearman 生物标志物与差异免疫浸润细胞的相关性。
setwd("/data/nas1/luchunlin/project/JNZK-214-8/16_supplement/")
if (! dir.exists("./DEtiic_cor")){
  dir.create("./DEtiic_cor")
}
setwd("./DEtiic_cor")
res.cibersort <- read.delim2('/data/nas1/luchunlin/project/JNZK-214-8/01_CIBERSORT/cibersort.txt',row.names = 1)%>%lc.tableToNum()
colnames(res.cibersort) <- gsub('.','-',colnames(res.cibersort),fixed = T)
cell <- read.delim2('/data/nas1/luchunlin/project/JNZK-214-8/02_CIBERSORT(GEO)/DEcell.xls')

hub_exp <- log2(hub.dat+1)

DEtiic <- res.cibersort[cell$ImmuneCell,]

DEtiic <- t(DEtiic)%>%as.data.frame()
cnt<-1
while(cnt<7){
  library(psych)
  x <-as.numeric(hub_exp[cnt,])
  y <-DEtiic
  library(psych)
  d <- corr.test(y,x,use="complete",method = 'pearson')
  r <- data.frame(d$r)
  p <- data.frame(d$p)
  write.table(p,file=paste(rownames(hub_exp)[cnt], "p.xls",sep='_'),sep="\t",quote=F)
  write.table(r,file=paste(rownames(hub_exp)[cnt], "r.xls",sep='_'),sep="\t",quote=F)
  correlation<-data.frame(rownames(p),r$d.r,p$d.p)
  colnames(correlation) <- c("cell","Correlation","p.value")
  correlation<-correlation[correlation$p.value<0.05,]
  correlation$sig[correlation$p.value>0.05] <- "ns"   
  correlation$sig[correlation$p.value<0.05&correlation$p.value>0.01] <- "*"   
  correlation$sig[correlation$p.value<0.01&correlation$p.value>0.001] <- "**"  
  correlation$sig[correlation$p.value<0.001&correlation$p.value>0.0001] <- "***"   
  correlation$sig[correlation$p.value<0.0001] <- "****"
  correlation$'correlation'<-abs(correlation$Correlation)
  #correlation_sig<-correlation[c(3,9,10,11,23,29,7,12,18,31,1,14,15,25,28),]
  library("viridis")
  #trace(ggdotchart,edit=T)
  p0<-ggdotchart(correlation, x = "cell", y = "Correlation",
                 
                 dot.size ='correlation',
                 # filled.contour(p.value,col=Lab.palette(20)),
                 color ='p.value',
                 # label='sig',
                 
                 #  font.label = list( size = 9),
                 
                 #ylab="mgp = c(3, 2, 0)",
                 # font.label = list(size=10,face="bold",color='black',position="dodge",
                 #                   vjust=0.5),
                 #  y.label=0.7,
                 # palette = palette,
                 # 按照cyl填充颜色
                 # palette = palette, # 修改颜色
                 sorting = "descending",
                 add = "segments",                             # 添加棒子
                 rotate = TRUE,
                 ggtheme = theme_pubr(),                        # 改变主题
                 ylab="Correlation Coefficient (r)",
                 xlab='',
                 # dot.size = 6 ,
                 title=rownames(hub_exp)[cnt]
                 
  )+scale_colour_gradient( high = "#4682B4",low = "#66CDAA")
  
  p10<-p0+theme(legend.position = "right",
                panel.background = element_blank())+geom_hline(aes(yintercept=0),linetype="dashed",lwd = 0.2)+
    theme(axis.title.x =element_text(size=15,family = "Times", face = "bold"),
          axis.text.x =element_text(size=12,family = "Times", face = "bold"),
          axis.title.y =element_text(size=20,family = "Times", face = "bold"),
          axis.text.y=element_text(size=12,family = "Times", face = "bold"),
          plot.title=element_text(size=20,family = "Times", face = "bold",hjust=0.5),
          legend.text = element_text(size = 16, family = "Times"),
          legend.title = element_text(size = 18, family = "Times",face = "bold"))+
    theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
  
  p10
  ggsave(paste(cnt,rownames(hub_exp)[cnt],'pdf',sep='.'),w=8,h=6)
  ggsave(paste(cnt,rownames(hub_exp)[cnt],'png',sep='.'),w=8,h=6)
  
  cnt<-cnt+1
}





