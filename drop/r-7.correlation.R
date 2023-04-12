rm(list = ls())
# 01 获取数据集--------------
setwd("/data/nas1/luchunlin/project/BJTC-356/")
if (! dir.exists("./07_correlation")){
  dir.create("./07_correlation")
}
setwd("./07_correlation")

##基因和蛋白的关系，蛋白和代谢物的关系

##基因表达
dat.fpkm <- read.delim2('../00_rawdata/dat.fpkm.xls',row.names = 1)%>%lc.tableToNum()
hubgene <- read_xlsx('../03_PPI/mcode.xlsx')

hub.exp <- dat.fpkm[hubgene$symbol,]
#hub.exp <- data.frame(row.names = rownames(hub.exp),Paroxysmal=rowMeans(hub.exp[,c(1:5)]),Persistent=rowMeans(hub.exp[,c(6:10)]))

hub.exp <- log2(hub.exp+1)

group <- read.delim2('../01_GSVA/group.xls')

colnames(hub.exp) <- group$group
##蛋白表达量表

dat.protein <- read.delim2('../00_rawdata/protein.dat.xls',row.names = 1)%>%lc.tableToNum()
hub.protein <- read_xlsx('../05_DEPs/mcode_protein.xlsx')

pro2gene <- read.delim2('../05_DEPs/deprotein.xls')
pro2gene <- pro2gene[pro2gene$SYMBOL%in%hub.protein$symbol,]

pro.exp <- dat.protein[pro2gene$UNIPROT,]
rownames(pro.exp) <- pro2gene$SYMBOL
pro.exp <- log2(pro.exp+1)

pro.mean <- data.frame(row.names = rownames(pro.exp),Paroxysmal=rowMeans(pro.exp[,c(1:3)]),Persistent=rowMeans(pro.exp[,c(4:6)]))

pro.exp <- cbind(pro.exp[,(1:3)],pro.mean$Paroxysmal,pro.mean$Paroxysmal,pro.exp[,c(4:6)],pro.mean$Persistent,pro.mean$Persistent)

colnames(pro.exp) <- colnames(hub.exp)


## 相关性------
library(ggcorrplot)
library(corrplot)

cor.dat <- rbind(hub.exp,pro.exp)
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
pdf(file = '03.correlation.pdf',w=9,h=9)
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
                    number.cex = 0.6
)
dev.off()
png(file = '03.correlation.png',w=800,h=800)
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



