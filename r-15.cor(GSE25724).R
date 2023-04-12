rm(list = ls())
setwd("/data/nas1/luchunlin/project/LZZK-519-10/")
if (! dir.exists("./15_cor(GSE25724)")){
  dir.create("./15_cor(GSE25724)")
}
setwd("./15_cor(GSE25724)")

## 相关性------
library(ggcorrplot)
library(corrplot)
dat <- read.delim2('../00_rawdata/dat(GSE25724).xls',row.names = 1)%>%lc.tableToNum()
hubgene <- read.delim2('../06_PPI/hubgene.xls')
hub.dat <- dat[hubgene$Symbol,]
cor.dat <- hub.dat
corr<-round(cor(t(cor.dat),method = 'spearman'),3)
write.table(corr,file = 'correlation.xls',
            sep = '\t',
            row.names = T)
## 检验基因之间的相关性p值
## 计算相关性系数并显示基因之间的相关性。相关性系数大于0为正相关，小于0为负相关。
## p小于0.05认为相关性显著
p.mat<-round(cor_pmat(t(cor.dat)),3)
col1 <- colorRampPalette(c("#4169E1","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF",
                           "#FDDBC7","#F4A582","#D6604D", "#B2182B","#CC0000","#990000"))
pdf(file = '01.correlation.pdf',w=4,h=4)
corr_plot<-corrplot(corr,
                    method = "circle",
                    is.corr = T,
                    type = "lower",
                    p.mat = p.mat,
                    insig = "blank",
                    outline = "white",
                    addCoef.col ="black",
                    col = col1(200),
                    tl.col = 'black',
                    tl.offset = 0.4,
                    number.font = 2,
                    cl.cex = 0.8,
                    number.cex = 0.8
)
dev.off()
png(file = '01.correlation.png',w=300,h=300)
corr_plot<-corrplot(corr,
                    method = "circle",
                    is.corr = T,
                    type = "lower",
                    p.mat = p.mat,
                    insig = "blank",
                    outline = "white",
                    addCoef.col ="black",
                    col = col1(200),
                    tl.col = 'black',
                    tl.offset = 0.4,
                    number.font = 2,
                    cl.cex = 0.8,
                    number.cex = 0.8
)
dev.off()
