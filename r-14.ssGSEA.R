rm(list = ls())
setwd("/data/nas1/luchunlin/project/HZ0301-3/")
if (! dir.exists("./14_ssGSEA")){
  dir.create("./14_ssGSEA")
}
setwd("./14_ssGSEA")

library(tidyverse)
dat<-read.csv('../00_rawdata/dat.fpkm.xls',sep = '\t')
colnames(dat)<-gsub('.','-',colnames(dat),fixed = T)
risk<-read.delim2('../07_risk/risk.xls')
high.sample<-risk$id[which(risk$risk==0)]
low.sample<-risk$id[which(risk$risk==1)]
dat<-dat[,risk$id]
group<-data.frame(sample=risk$id,group=ifelse(risk$risk=='0','High risk','Low risk'))

## 03-1 ssGSEA--------
library(GSVA)
gene_set <- read.table("/data/nas1/luchunlin/pipeline/ssGSEA/mmc3.txt",
                       header = T,
                       sep ="\t")
dat.final2 <- as.matrix(dat)
gene_list <- split(as.matrix(gene_set)[,1],
                   gene_set[,2])

ssgsea_score = gsva(dat.final2, gene_list, 
                    method = "ssgsea", 
                    ssgsea.norm = TRUE, 
                    verbose = TRUE)
write.table(ssgsea_score,
            file = "ssgsea_result.xls",
            sep = "\t",
            quote = F)
## 富集分数画热图
group<-group[order(group$group),]
ssgsea_score<-ssgsea_score[,group$sample]
annotation_col<-as.data.frame(group$group)
colnames(annotation_col)='Group'
rownames(annotation_col)=colnames(ssgsea_score)

color.key<-c("#3300CC", "#3399FF", "white", "#FF3333", "#CC0000")
ann_colors<-list(
  Group = c('Low risk'="#00FFFF",'High risk'="#FFAEB9"))
pheatmap(
  ssgsea_score,
  color = colorRampPalette(color.key)(50),
  border_color = 'darkgrey',
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  labels_row = NULL,
  clustering_method = 'ward.D2',
  show_rownames = T,
  show_colnames = F,
  fontsize_col = 5,
  cluster_cols = F,
  cluster_rows = T)
## 03-2 差异免疫细胞鉴定-------
tiics_result <- ssgsea_score
pvalue = padj = log2FoldChange <- matrix(0, nrow(tiics_result), 1)
group_High<-group[group$group=='High risk',]
group_High<-as.character(group_High$sample)
group_Low<-group[group$group=='Low risk',]
group_Low<-as.character(group_Low$sample)

for (i in 1:nrow(tiics_result)){
  pvalue[i, 1] = p.value = wilcox.test(tiics_result[i, group_High],
                                       tiics_result[i, group_Low])$p.value
  log2FoldChange[i, 1] = mean(tiics_result[i, group_High]) - 
    mean(tiics_result[i, group_Low])
}
padj <- p.adjust(as.vector(pvalue), "fdr", n = length(pvalue))
rTable <- data.frame(log2FoldChange, 
                     pvalue, 
                     padj,
                     row.names = rownames(tiics_result))
Low <- signif(apply(tiics_result[rownames(rTable), group_Low], 
                    1,
                    mean), 4)
High <- signif(apply(tiics_result[rownames(rTable), group_High], 
                     1, 
                     mean), 4)
rTable <- data.frame(Low, 
                     High,
                     rTable[, c("padj", "pvalue", "log2FoldChange")])
rTable$immune_cell <- rownames(rTable)
rTable$sig <- ifelse(rTable$pvalue < 0.05,
                     ifelse(rTable$pvalue < 0.01, 
                            ifelse(rTable$pvalue < 0.001,
                                   ifelse(rTable$pvalue < 0.0001,
                                          paste(rTable$immune_cell, "****",  sep = ""),
                                          paste(rTable$immune_cell, "***", sep = "")),
                                   paste(rTable$immune_cell, "**", sep = "")),
                            paste(rTable$immune_cell, "*",  sep = "")), 
                     rTable$immune_cell)

diff_Table<-rTable[which(rTable$pvalue<0.05),]
## 19
write.table(rTable,
            file = "tiics_wilcox_test.xls",
            quote = F,
            row.names = F,
            sep = '\t')
write.table(diff_Table,
            file = "diff_tiics_wilcox_test.xls",
            quote = F,
            row.names = F,
            sep = '\t')

### 箱线图----
library(tidyr)
library(ggplot2)
library(ggpubr)
library(Ipaper)

xCell2 <- data.frame(Immune_Cell=rownames(tiics_result), 
                     tiics_result, 
                     pvalue=rTable$pvalue)
# xCell3 <- xCell2[which(xCell2$pvalue<0.05),]
xCell3<-xCell2
diff_tiics <- rownames(xCell3)
violin_dat <- gather(xCell3, key=Group, value=score, -c("Immune_Cell","pvalue"))

violin_dat$Group <- ifelse(gsub("\\.","-",violin_dat$Group) %in% group_Low,
                           "Low", "High") 
head(violin_dat)
boxplot_diff_TIICs <- ggplot(violin_dat, aes(x=Immune_Cell, 
                                             y=score,
                                             fill=Group)) +
  # geom_violin(trim=T,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)#"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar",
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  geom_point(aes(fill = Group),
             size = 0.05,
             position = position_dodge(0.9))+
  scale_fill_manual(values= c("#48D1CC", "#FF7256"))+ #设置填充的颜色
  labs(title="", x="", y = "Score",size=20) +
  stat_compare_means(data = violin_dat,
                     mapping = aes(group = Group),
                     label ="p.signif",
                     hide.ns = F) +
  theme_bw()+#把背景设置为白底
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18), # 将图表标题居中
        axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=10), #设置x轴刻度标签的字体显示倾斜角度为45度，并向下调整1(hjust = 1)，字体大小为14
        axis.text.y=element_text(hjust=0.5,colour="black",size=12), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.x=element_text(size=16,face="bold"),#设置x轴标题的字体属性
        axis.title.y=element_text(size=14,face="bold"), #设置y轴标题的字体属性
        legend.text=element_text(face="bold", hjust = 0.5,colour="black", size=11), #设置图例的子标题的字体属性
        legend.title=element_text(face="bold", colour="black", size=11),#设置图例的总标题的字体属性
        #legend.justification=c(-0.1,1.2), #可调整图例的位置。##(1,1)第一个1是调整图例在图的内外(左右移动)，第二个1是在图中上下移动。
        #legend.position=c(0, 1.04), #legend.position=c(0,1)左上角，(1,1)是在右上角。
        panel.grid.major = element_blank(), #不显示网格线
        panel.grid.minor = element_blank()) #不显示网格线
boxplot_diff_TIICs

ggsave(filename = '01.DEcells.pdf',boxplot_diff_TIICs,w=10,h=5.5)
ggsave(filename = '01.DEcells.png',boxplot_diff_TIICs,w=11,h=6)

###相关性1---------
hubgene <- read.delim2('../06_Lasso/lasso_genes.csv',header = F)
hub.exp<-dat[hubgene$V1,]
library(Hmisc)
hub_exp <- log2(hub.exp+1)
cortiic <- ssgsea_score[rownames(diff_Table),]
cortiic <- cortiic[,colnames(hub_exp)]

corr.dat<-t(rbind(cortiic,hub_exp))

gene_cor <- rcorr(corr.dat,type = 'spearman')$r[1:nrow(cortiic),(ncol(corr.dat)-length(rownames(hub_exp))+1):ncol(corr.dat)]
m <- t(gene_cor)
#将获得的相关性矩阵转换为两两对应的数据框结构
gene_cor <- reshape2::melt(gene_cor)
#gene_cor <- subset(gene_cor, value != 0)  #去除0值的相关性
head(gene_cor)  #前两列是两个基因名称，第三列为两个基因的相关性
colnames(gene_cor) <- c('symbol','tiic','correlation')
## 计算p值
gene_p <- rcorr(corr.dat,type = 'spearman')$P[1:nrow(cortiic),(ncol(corr.dat)-length(rownames(hub_exp))+1):ncol(corr.dat)]
p <- t(gene_p)
gene_p <- reshape2::melt(gene_p)
head(gene_p)  #前两列是两个基因名称，第三列为
gene_cor$pvalu <- gene_p$value
write.table(gene_cor,file = 'correlation.gene.xls',sep = '\t',row.names = F,quote = F)
##|Cor| > 0.4, P< 0.001
# cor.final <- subset(gene_cor,gene_cor$pvalu<0.05 & abs(gene_cor$correlation)>0.5)
## 13
# write.table(cor.final,file = 'correlation.final.xls',sep = '\t',row.names = F,quote = F)

## 相关性热图---------
library(dplyr)
tmp <- p
tmp[tmp<0.0001] <- '****'
tmp[tmp<0.001] <- '***'
tmp[tmp<0.01] <- '**'
tmp[tmp<0.05] <- '*'
tmp[tmp>0.05] <- ''

cor <- m
cor <- signif(cor,3)
#cor[abs(cor)<0.1] <- ''

textMatrix = paste(cor,"\n",
                   tmp, sep = "")
dim(textMatrix) = dim(m)
colnames(textMatrix)<-colnames(m)
rownames(textMatrix)<-rownames(m)
library(WGCNA)
pdf(file = '02.correlation(modelgene).pdf',w=12,h=6)
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
               main = paste("modelgene-DEcells correlation"))
dev.off()
png(file = '02.correlation(modelgene).png',w=800,h=400)
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
               main = paste("modelgene-DEcells correlation"))
dev.off()

##最正最负2个
##MYB-MemoryB   LRFN-EffECTOR MEMEORY CD4T
library(ggplot2)
library(ggpubr)
library(ggExtra)

rownames(DEtiic)
cell <- as.numeric(cortiic["Memory B cell",])
hubgene <- as.numeric(hub_exp['MYB',])
df1 <- cbind(hubgene,cell)%>%as.data.frame()
p1 <- ggplot(df1,aes(cell,hubgene))+
  xlab('Memory B cell')+ylab('MYB')+
  geom_point()+geom_smooth(method = 'lm',formula = y~x)+
  theme_bw()+
  stat_cor(method = 'spearman',aes(x=cell,y=hubgene))
p1  

p2 <- ggMarginal(p1,type = 'density',xparams = list(fill='orange'),
                 yparams = list(fill='blue'))

p2

ggsave(filename = '03.cor1.pdf',p2,w=4,h=4)
ggsave(filename = '03.cor1.png',p2,w=4,h=4)

cell <- as.numeric(cortiic["Effector memeory CD4 T cell",])
hubgene <- as.numeric(hub_exp['LRFN4',])
df1 <- cbind(hubgene,cell)%>%as.data.frame()
p3 <- ggplot(df1,aes(cell,hubgene))+
  xlab('Effector memeory CD4 T cell')+ylab('LRFN4')+
  geom_point()+geom_smooth(method = 'lm',formula = y~x)+
  theme_bw()+
  stat_cor(method = 'spearman',aes(x=cell,y=hubgene))
p3  
p4 <- ggMarginal(p3,type = 'density',xparams = list(fill='orange'),
                 yparams = list(fill='blue'))

p4

ggsave(filename = '04.cor2.pdf',p4,w=4,h=4)
ggsave(filename = '04.cor2.png',p4,w=4,h=4)


## 相关性---------
DEtiic <- ssgsea_score[rownames(diff_Table),]
DEtiic <- t(DEtiic)%>%as.data.frame()
DEtiic <- DEtiic[risk$id,]

library(psych)
x <-as.numeric(risk$riskScore)
y <-DEtiic
library(psych)
d <- corr.test(y,x,use="complete",method = 'spearman')
r <- data.frame(d$r)
p <- data.frame(d$p)
correlation<-data.frame(rownames(p),r$d.r,p$d.p)
colnames(correlation) <- c("cell","Correlation","p.value")
write.table(correlation,file = 'correlation.xls',sep = '\t',row.names = F,quote = F)
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
               title='riskScore'
               
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
ggsave(filename = '05.cor.pdf',p10,w=8,h=5)  
ggsave(filename = '05.cor.png',p10,w=8,h=5)


