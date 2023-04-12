## 表达分析----------
rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-334")
if (! dir.exists("./07_circos")){
  dir.create("./07_circos")
}
setwd("./07_circos")

library(lance)
library(tidyverse)
library(circlize)
dat = read.delim2("/data/nas1/luchunlin/project/BJTC-334/00_rawdata/dat.final.xls", row.names = 1) %>% lc.tableToNum
dat<-log2(dat+1)
group = read.delim2("/data/nas1/luchunlin/project/BJTC-334/00_rawdata/group.xls")
control.sample<-group$sample[which(group$group=='Control')]
hubgene<-read.delim2('/data/nas1/luchunlin/project/BJTC-334/05_expression/hubgene.xls')
dat.final<-dat[hubgene$hubgene,]
dat.final$gene<-rownames(dat.final)
dat.final<-dat.final[,c(16,1:15)]
control.dat<-dat.final[,c(1:10)]
CAD.dat<-dat.final[,c(1,11:16)]
circos.dat1<-gather(control.dat,key = sample,value = 'expr',-c("gene"))
circos.dat2<-gather(CAD.dat,key = sample,value = 'expr',-c("gene"))

circos.dat1$expr<-as.numeric(circos.dat1$expr)
circos.dat1$gene<-as.character(circos.dat1$gene)
circos.dat2$expr<-as.numeric(circos.dat2$expr)
circos.dat2$gene<-as.character(circos.dat2$gene)
## 添加相关性
#表达值进行log(1+)转化，使数据更服从正态分布，减少离散度极大值影响
gene <- t(dat.final[,-1])
#基因表达值的相关性分析，以Pearson相关系数为例
gene_cor <- cor(gene, method = 'pearson')
#去除基因的自相关，也就是对角线的值
diag(gene_cor) <- 0
gene_cor  #最终的基因间表达值Pearson相关性矩阵
#将获得的相关性矩阵转换为两两对应的数据框结构
gene_cor <- reshape2::melt(gene_cor)
gene_cor <- subset(gene_cor, value != 0)  #去除0值的相关性
head(gene_cor)  #前两列是两个基因名称，第三列为两个基因的相关性

circos.par("track.height" = 0.05)
circos.initialize(circos.dat1$gene, x = circos.dat1$expr)
bgcol = c("#FF69B4", "#9370DB","#FFA07A","#20B2AA","#87CEFA","#48D1CC","#FFC125","#3CB371","#BC8F8F","#DDA0DD")
circos.track(circos.dat1$gene, y = circos.dat1$expr,
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, 
                           CELL_META$cell.ylim[2] + mm_y(10), 
                           CELL_META$sector.index)
               circos.axis(labels.cex = 0.6)
             },
             bg.col=bgcol)

circos.par("track.height" = 0.2)
# circos.trackHist(circos.dat1$gene, x = circos.dat1$expr, col = "#999999", 
#                  border = "#999999")
circos.trackHist(circos.dat1$gene, x = circos.dat1$expr, bin.size = 0.004, 
                 col = "#5599FF", border = "#999999")
# circos.trackHist(circos.dat1$gene, x = circos.dat1$expr, draw.density = TRUE, 
#                  col = "#999999", border = "#999999")
circos.trackHist(circos.dat2$gene, x = circos.dat2$expr, bin.size = 0.005, 
                 col = "#FFCCCC", border = "#999999")

circos.par("track.height" = 0.2)
chordDiagram(gene_cor, 
             col = colorRamp2(c(-1, 0, 1), c('green', 'white', 'red'), transparency = 0.5), #根据相关性大小展示连线的颜色范围
)

