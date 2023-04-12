rm(list = ls())
setwd("/data/nas1/luchunlin/project/JNZK-255-2/")
if (! dir.exists("./09_CIBERSORT")){
  dir.create("./09_CIBERSORT")
}
setwd("./09_CIBERSORT")

library(lance)
dat<-read.delim2("../00_rawdata/dat.fpkm.xls", row.names = 1)%>% lc.tableToNum
colnames(dat)<-gsub('.','-',colnames(dat),fixed = T)
group<-read.delim2('../06_risk/risk.xls')%>%dplyr::select(c('id','risk'))
colnames(group)<-c('sample','group')
group$group <- ifelse(group$group=='1','Low risk','High risk')
Low.sample<-group$sample[which(group$group=='Low risk')]
High.sample<-group$sample[which(group$group=='High risk')]
dat <- dat[,group$sample]
library(immunedeconv)
set_cibersort_binary("CIBERSORT.R")
set_cibersort_mat("LM22.txt")
res.cibersort<-deconvolute(as.matrix(dat),method = 'cibersort')

write.table(res.cibersort,'cibersort.txt',sep = '\t',col.names = T,row.names = F,quote = F)
save(file = 'cibersort.Rdata',res.cibersort)
library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(7,"Paired"))
##画图
pdf('01.cibersort.box.pdf',w=10,h=8)
res.cibersort %>%
  gather(sample, fraction, -cell_type) %>%
  merge(group,by='sample')%>%
  # 绘制堆积条形图
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
  geom_bar(position = 'stack',stat = 'identity')+
  scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  labs(x='',
       y='Relative Percent',
       fill='')+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'top') +
  scale_fill_manual(values = mypalette(22))+
  facet_grid(~group,scales= "free",space= "free")
dev.off()
png('01.cibersort.box.png',w=700,h=500)
res.cibersort %>%
  gather(sample, fraction, -cell_type) %>%
  merge(group,by='sample')%>%
  # 绘制堆积条形图
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
  geom_bar(position = 'stack',stat = 'identity')+
  scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  labs(x='',
       y='Relative Percent',
       fill='')+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'top') +
  scale_fill_manual(values = mypalette(22))+
  facet_grid(~group,scales= "free",space= "free")
dev.off()

## 差异-------
##去掉在75%的样本中结果为0的细胞
res.cibersort2 <- res.cibersort%>%column_to_rownames(var = 'cell_type')
keep <- rowSums(res.cibersort2>0)>=floor(0.75*ncol(res.cibersort2))
res.cibersort2 <- res.cibersort2[keep,]

dat.cibersort <- res.cibersort2 %>% 
  # tibble::column_to_rownames(var = "cell_type") %>% 
  t %>% as.data.frame() %>% 
  tibble::rownames_to_column(var = "sample")
dat.cibersort <- merge(group, dat.cibersort, by = "sample")
dat.cibersort2 <- tidyr::gather(dat.cibersort, ImmuneCell, Score, -c("sample", "group"))
library(rstatix)
library(ggplot2)
library(ggpubr)
colnames(dat.cibersort2)
stat_cibersort <- dat.cibersort2 %>% 
  group_by(ImmuneCell) %>%
  wilcox_test(Score ~ group) %>% 
  adjust_pvalue(method = "BH") %>%  # method BH == fdr
  add_significance("p")
write.table(stat_cibersort,file = 'stat.cibersort.xls',sep = '\t',row.names = F,quote = F)
DE.cibersort<-stat_cibersort[which(stat_cibersort$p<0.05),]
write.table(DE.cibersort,file = 'DE.cibersort.xls',sep = '\t',row.names = F,quote = F)
colnames(dat.cibersort2)
violin.cibersort<-dat.cibersort2[dat.cibersort2$ImmuneCell%in%stat_cibersort$ImmuneCell[which(stat_cibersort$p<0.05)],]
cibersort_plot <- ggplot(violin.cibersort, aes(x=group,
                                               y=Score,
                                               fill=group)) +
  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar",
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.3,
               position=position_dodge(0.9),
               outlier.shape = NA,fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#FF6A6A", "#20B2AA"), name = "Group")+
  labs(title="Immune Cell", x="", y = "Score",size=20) +
  stat_compare_means(data = violin.cibersort,
                     mapping = aes(group = group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F,label.x = 1.4) +
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18),
        axis.text.x=element_text(angle=0,hjust=,colour="black",face="bold",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+facet_wrap(~ImmuneCell,scales = "free",nrow = 2) 
cibersort_plot
ggsave(filename = '02.cibersort.plot.pdf',cibersort_plot,w=10,h=9)
ggsave(filename = '02.cibersort.plot.png',cibersort_plot,w=10,h=9)

## 相关性---------

library(ggcorrplot)
library(corrplot)
de.dat <- res.cibersort%>%column_to_rownames(var = 'cell_type')
de.dat <- de.dat[DE.cibersort$ImmuneCell,]

cor.dat <- de.dat
library(Hmisc)
nc<-t(rbind(cor.dat,de.dat))
m=rcorr(nc,type = 'spearman')$r[1:nrow(cor.dat),(ncol(nc)-length(rownames(de.dat))+1):ncol(nc)]
m<-t(m)

p=rcorr(nc,type = 'spearman')$P[1:nrow(cor.dat),(ncol(nc)-length(rownames(de.dat))+1):ncol(nc)]
p<-t(p)

tmp <- p
tmp[tmp<0.0001] <- '****'
tmp[tmp<0.001] <- '***'
tmp[tmp<0.01] <- '**'
tmp[tmp<0.05] <- '*'
tmp[tmp>0.05] <- 'ns'

cor <- m
cor <- round(cor,digits = 3)
# cor[abs(cor)<0.15] <- ''

textMatrix = paste(cor,"\n",
                   tmp, sep = "")
dim(textMatrix) = dim(m)
colnames(textMatrix)<-colnames(m)
rownames(textMatrix)<-rownames(m)
write.table(textMatrix,file = 'correlaion(tiic).xls',sep = '\t',row.names = T,quote = F)

library(WGCNA)
pdf(file = '03.correlation.pdf',w=6,h=5)
par(mar = c(10, 10, 3, 1));
labeledHeatmap(Matrix = m, 
               xLabels = colnames(m), 
               yLabels = colnames(m), 
               cex.lab = 0.8, 
               # ySymbols = colnames(MEs_col)[-20], 
               colorLabels = FALSE, 
               colors = blueWhiteRed(50, 0.8),
               textMatrix = textMatrix, 
               setStdMargins = FALSE, 
               cex.text = 0.8, 
               zlim = c(-1,1),
               # main = paste("correlation")
)
dev.off()
png(file = '03.correlation.png',w=450,h=350)
par(mar = c(10, 10, 3, 1));
labeledHeatmap(Matrix = m, 
               xLabels = colnames(m), 
               yLabels = colnames(m), 
               cex.lab = 0.8, 
               # ySymbols = colnames(MEs_col)[-20], 
               colorLabels = FALSE, 
               colors = blueWhiteRed(50, 0.8),
               textMatrix = textMatrix, 
               setStdMargins = FALSE, 
               cex.text = 0.8, 
               zlim = c(-1,1),
               # main = paste("correlation")
)
dev.off()



###相关性1---------
hubgene <- read.delim2('../05_Lasso/lasso_genes.csv',header = F)
hub.exp<-dat[hubgene$V1,]
library(Hmisc)
hub_exp <- log2(hub.exp+1)
cortiic <- res.cibersort%>%column_to_rownames(var = 'cell_type')
cortiic <- cortiic[DE.cibersort$ImmuneCell,]
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
colnames(gene_cor)
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
tmp[tmp>0.05] <- 'ns'

cor <- m
cor <- signif(cor,3)
#cor[abs(cor)<0.1] <- ''

textMatrix = paste(cor,"\n",
                   tmp, sep = "")
dim(textMatrix) = dim(m)
colnames(textMatrix)<-colnames(m)
rownames(textMatrix)<-rownames(m)
library(WGCNA)
pdf(file = '04.correlation(modelgene).pdf',w=7,h=4)
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
               main = paste("modelgenes-DEcells correlation"))
dev.off()
png(file = '04.correlation(modelgene).png',w=500,h=300)
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
               main = paste("modelgenes-DEcells correlation"))
dev.off()


