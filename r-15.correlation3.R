rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/GY0324-12/")
if (! dir.exists("./15_correlation3")){
  dir.create("./15_correlation3")
}
setwd("./15_correlation3")

library(lance)
library(tidyverse)
datmeta <- read.delim2('/data/nas1/luchunlin/project/GY0324-12/00_rawdata/dat.meta.xls')%>%lc.tableToNum()
colnames(datmeta)
### A3 B4------

datA <- read.delim2('/data/nas1/luchunlin/project/GY0324-12/00_rawdata/datA3.xls')%>%lc.tableToNum()
datA<- datA[,c(1:3)]
datB <- read.delim2('/data/nas1/luchunlin/project/GY0324-12/00_rawdata/datB4.xls')%>%lc.tableToNum()
datB <- datB[,c(1:3)]
dat1 <- cbind(datA,datB)

sig.diff1 <- read.delim2('/data/nas1/luchunlin/project/GY0324-12/13_DEGs3/DEG_sig(A3vsB4).xls')
sig.dat1 <- dat1[rownames(sig.diff1),]
datmeta1 <- datmeta[,c(21:40)]
colnames(sig.dat1)
colnames(datmeta1)
datmeta1 <- datmeta1[,c("A3_26_1","A3_58_1","A3_86_1","B_104_1","B_120_1","B_2_1")]
colnames(datmeta1) <- c("A3_26","A3_58","A3_86","B_104","B_120","B_2")

colnames(sig.dat1) <- colnames(datmeta1)
meta.sig1 <- read.delim2('/data/nas1/luchunlin/project/GY0324-12/01_DEMs/DEM_sig(B4vs.C5).xls')
sig.meta1 <- datmeta1[rownames(meta.sig1),]

library(Hmisc)
corr.dat<-t(rbind(sig.dat1,sig.meta1))
#基因表达值的相关性分析
gene_cor <- rcorr(corr.dat,type = 'spearman')$r[1:nrow(sig.dat1),(ncol(corr.dat)-length(rownames(sig.meta1))+1):ncol(corr.dat)]

#将获得的相关性矩阵转换为两两对应的数据框结构
gene_cor <- reshape2::melt(gene_cor)
#gene_cor <- subset(gene_cor, value != 0)  #去除0值的相关性
head(gene_cor)  #前两列是两个基因名称，第三列为两个基因的相关性
colnames(gene_cor) <- c('DEG','DEM','correlation')
## 计算p值
gene_p <- rcorr(corr.dat,type = 'spearman')$P[1:nrow(sig.dat1),(ncol(corr.dat)-length(rownames(sig.meta1))+1):ncol(corr.dat)]
gene_p <- reshape2::melt(gene_p)
head(gene_p)  #前两列是两个基因名称，第三列为
gene_cor$pvalue <- gene_p$value

write.table(gene_cor,file = '01.correlation.all.xls',sep = '\t',row.names = F,quote = F)
##|Cor| > 0.8, P< 0.05
cor.final <- subset(gene_cor,gene_cor$pvalue<0.05 & abs(gene_cor$correlation)>0.8)

write.table(cor.final,file = '01.correlation.final.xls',sep = '\t',row.names = F,quote = F)

### 九向限------
genefc1 <- data.frame(symbol=rownames(sig.diff1),logFC_DEGs=sig.diff1$logFC)
data1 <- merge(cor.final,genefc1,by.x ='DEG',by.y = 'symbol')%>%select(-c('correlation','pvalue'))
genefc2 <- data.frame(DEM=rownames(meta.sig1),logFC_DEMs=meta.sig1$log2.FC.)

data.merge1 <- merge(data1,genefc2,by='DEM')

#保存合并后的数据；
# write.csv(data.merge1,"RNA_META.csv",row.names=FALSE)

library(dplyr)
library(ggplot2)
library(ggrepel)

data <- data.merge1
data$logFC_DEGs <- as.numeric(data$logFC_DEGs)
data$logFC_DEMs <- as.numeric(data$logFC_DEMs)
#对数据进行分组；
#生成显著上下调数据标签；
data$part <- case_when(abs(data$logFC_DEGs) >= 1 & abs(data$logFC_DEMs) >= 1 ~ "part1379",
                       abs(data$logFC_DEGs) < 1 & abs(data$logFC_DEMs) > 1 ~ "part28",
                       abs(data$logFC_DEGs) > 1 & abs(data$logFC_DEMs) < 1 ~ "part46",
                       abs(data$logFC_DEGs) < 1 & abs(data$logFC_DEMs) < 1 ~ "part5")

data$part <- case_when(data$logFC_DEGs <= 1 & data$logFC_DEMs >= 1 ~ "part1",
                       data$logFC_DEGs >= 1 & data$logFC_DEMs >= 1 ~ "part3",
                       data$logFC_DEGs <= -1 & data$logFC_DEMs <= 1 ~ "part7",
                       data$logFC_DEGs < -1 & data$logFC_DEMs > 1 ~ "part2",
                       data$logFC_DEGs < -1 & data$logFC_DEMs < 1 ~ "part8",
                       abs(data$logFC_DEGs) > 1 & abs(data$logFC_DEMs) < 1 ~ "part46",
                       abs(data$logFC_DEGs) < -1 & abs(data$logFC_DEMs) < 1 ~ "part5",
                       data$logFC_DEGs > 1 & data$logFC_DEMs < 1 ~ "part9")



table(data$part) 
# part1 part3 part7 part9 
#   10    11  2134  4534  
head(data)
write.csv(data,"RNA_META.csv",row.names=FALSE)

#开始尝试绘图；
p0 <-ggplot(data,aes(logFC_DEGs,logFC_DEMs,color=part))
#添加散点；
p1 <- p0+geom_point(size=1)+guides(color="none")
p1
#改变点颜色
mycolor <- c("red","darkgreen","purple","blue")
p2 <- p1 + scale_colour_manual(name="",values=alpha(mycolor,0.7))
p2
#添加辅助线；
p3 <- p2+geom_hline(yintercept = c(-1,1),
                    size = 1,
                    color = "grey40",
                    lty = "dashed")+
  geom_vline(xintercept = c(-1,1),
             size = 1,
             color = "grey40",
             lty = "dashed")
p3
#调整横轴和纵轴绘图区域的范围；
#设置y轴范围（上下两端的空白区域设为1），修改刻度标签；
#expansion函数设置坐标轴范围两端空白区域的大小；mult为“倍数”模式，add为“加性”模式；
p4<-p3+
  scale_y_continuous(expand=expansion(add = c(0.1, 0.1)),
                     limits = c(-20, 20)
                     # breaks = c(-0.3,-3,0,3,6),
                     # label = c("-6","-3","0","3","6")
  )+
  scale_x_continuous(expand=expansion(add = c(0.5, 0.5)),
                     limits = c(-10, 10)
                     # breaks = c(-6,-3,0,3,6),
                     # label = c("-6","-3","0","3","6")
  )
p4
#自定义图表主题，对图表做精细调整；
top.mar=0.2
right.mar=0.2
bottom.mar=0.2
left.mar=0.2
#隐藏纵轴，并对字体样式、坐标轴的粗细、颜色、刻度长度进行限定；
mytheme<-theme_bw()+
  theme(text=element_text(family = "sans",colour ="gray30",size = 17),
        panel.grid = element_blank(),
        panel.border = element_rect(size = 0.8,colour = "gray30"),
        axis.line = element_blank(),
        axis.ticks = element_line(size = 0.6,colour = "gray30"),
        axis.ticks.length = unit(1.5,units = "mm"),
        plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
                         units="inches"),
  )
#应用自定义主题；
p5 <- p4+mytheme
p5

p6 <- p5+labs(x='log2 ratio of gene',y='log2 ratio of metabolin')
p6
ggsave(filename = "01.nqd.png", width = 7, height = 6)
ggsave(filename = "01.nqd.pdf", width = 7, height = 6)
dev.off()
