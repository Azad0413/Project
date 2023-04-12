##  差异基因鉴定-----
rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-370-8/")
if (! dir.exists("./03_DEGs")){
  dir.create("./03_DEGs")
}
setwd("./03_DEGs")
library(tidyverse)
library(lance)
##合并去批次-------
## 样品校正前后PCA,使用PCA检查有无批次效应 ## 样品校正前后（Boxplot）
## 批次信息 一列为样本名字与表达矩阵一致，另一列为批次信息，用1,2表示
dat1 <- read.delim2('../00_rawdata/dat(GSE110811).xls')%>%rownames_to_column(var = 'symbol') %>%lc.tableToNum()
dat2 <- read.delim2('../00_rawdata/dat(GSE97508).xls')%>%rownames_to_column(var = 'symbol')%>%lc.tableToNum()
dat3 <- read.delim2('../00_rawdata/dat(GSE24673).xls')%>%rownames_to_column(var = 'symbol')%>%lc.tableToNum()
group1 <- read.delim2('../00_rawdata/group(GSE110811).xls')
group2 <- read.delim2('../00_rawdata/group(GSE97508).xls')
group3 <- read.delim2('../00_rawdata/group(GSE24673).xls')
dat <- merge(dat1,dat2,by='symbol')
dat <- merge(dat,dat3,by='symbol')%>%column_to_rownames(var = 'symbol')
group <- rbind(group1,group2)
group <- rbind(group,group3)

library(sva)
library(bladderbatch)
## 转化为矩阵
sif<-data.frame(batch=c(rep('1',34),rep('2',9),rep('3',11)))
dat <- t(dat)
dat_pca<-cbind(sif,dat)
dat_pca<-as.data.frame(lapply(dat_pca,as.numeric))
dat_pca<-t(dat_pca)
modcombat=model.matrix(~1,data=sif)
combat=ComBat(dat = dat_pca,
              batch = sif$batch,
              mod = modcombat,
              par.prior = T,
              prior.plots = F)
colnames(combat)<-rownames(dat)

dat_final<-combat[-1,]
dat_final<-t(dat_final)
#绘图风格
library(ggplot2)
pdf('01.expression_before.pdf',w=10,h=5)
par(mar = c(8, 3, 3, 3))
boxplot(t(dat),las=2,
        notch = F, col = rainbow(10))
dev.off()

png('01.expression_before.png',w=800,h=400)
par(mar = c(8, 3, 3, 3))
boxplot(t(dat),las=2,
        notch = F, col = rainbow(10))
dev.off()

pdf('02.expression_after.pdf',w=10,h=5)
par(mar = c(8, 3, 3, 3))
boxplot(t(dat_final),las=2,
        notch = F, col = rainbow(10))
dev.off()

png('02.expression_after.png',w=800,h=400)
par(mar = c(8, 3, 3, 3))
boxplot(t(dat_final),las=2,
        notch = F, col = rainbow(10))
dev.off()
##PCA
library(scatterplot3d)

pca_pre<-prcomp(dat,center = TRUE,scale. = TRUE)
df_pre<-pca_pre$x %>%as.data.frame()
group.pca <- data.frame(sample=group$sample,group=c(rep('GSE110811',34),rep('GSE97508',9),rep('GSE24673',11)))
group.pca$group <- factor(group.pca$group,levels = c('GSE110811','GSE97508','GSE24673'))
shapes <- c(16,17,18)
colors <- c("#32CD32", "#0000FF", "#FF4500")
colors <- colors[as.numeric(group.pca$group)]
colors
shapes <- shapes[(group.pca$group)]
pdf(file = '03.PCA_before.pdf',w=5,h=5)
scatterplot3d(df_pre[,1:3],pch = shapes,color = colors)
dev.off()

png(file = '03.PCA_before.png',w=400,h=400)
scatterplot3d(df_pre[,1:3],pch = shapes,color = colors)
dev.off()

pca_after<-prcomp(dat_final,center = TRUE,scale. = TRUE)
df_after<-pca_after$x %>%as.data.frame()
group.pca <- data.frame(sample=group$sample,group=c(rep('GSE110811',34),rep('GSE97508',9),rep('GSE24673',11)))
group.pca$group <- factor(group.pca$group,levels = c('GSE110811','GSE97508','GSE24673'))
shapes <- c(16,17,18)
colors <- c("#32CD32", "#0000FF", "#FF4500")
colors <- colors[as.numeric(group.pca$group)]
colors
shapes <- shapes[(group.pca$group)]
pdf(file = '04.PCA_after.pdf',w=5,h=5)
scatterplot3d(df_after[,1:3],pch = shapes,color = colors)
dev.off()

png(file = '04.PCA_after.png',w=400,h=400)
scatterplot3d(df_after[,1:3],pch = shapes,color = colors)
dev.off()

library(lance)
library(tidyverse)
library(magrittr)
library(stringr)
library(limma)
df = t(dat_final)
df.group = group[order(group$group),]
write.table(df.group,file = 'group.xls',sep = '\t',row.names = F,quote = F)
df = df[,df.group$sample]
write.table(df,file = 'dat_final.xls',sep = '\t',row.names = T,quote = F)
table(df.group$group)
df.group$group = factor(df.group$group, levels = c("control", "RB"))
design.mat = cbind(control = ifelse(df.group$group == "control", 1, 0), 
                   RB = ifelse(df.group$group == "control", 0, 1))
contrast.mat = makeContrasts(contrasts="RB-control", levels=design.mat)

fit = lmFit(df, design.mat)
fit = contrasts.fit(fit, contrast.mat)
fit = eBayes(fit)
fit = topTable(fit, coef = 1, number = Inf, adjust.method = "fdr")
#fit = fit[c(1,4,5)]
DEG=na.omit(fit)
logFC_cutoff <- 0
DEG$change = as.factor(
  ifelse(DEG$adj.P.Val<0.05 & abs(DEG$logFC) > logFC_cutoff,
         ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
sig_diff <- subset(DEG,
                   DEG$adj.P.Val < 0.05 & abs(DEG$logFC) > logFC_cutoff)

dim(DEG)
dim(sig_diff)
## 1919  6
summary(sig_diff$change)
# DOWN  NOT   UP 
#  1133    0  786 
write.table(DEG,file = "DEG_all.xls",quote = F,sep = "\t",row.names = T)
write.table(sig_diff,file = "DEG_sig.xls",quote = F,sep = "\t",row.names = T)

ferr <- read.table('/data/nas2/database/gene_set/Ferr_genes.txt',header = T)
siginter <- intersect(rownames(sig_diff),ferr$Ferr_gene)%>%as.data.frame()

scDEG <- read.delim2('../02_scRNA/02_PCA/DEGs.xls')
modgene <- read.delim2('../01_WGCNA/modGene.xls')

inter <- intersect(modgene$modgene,scDEG$gene)%>%as.data.frame()
inter <- intersect(inter$.,ferr$Ferr_gene)%>%as.data.frame()
inter <- intersect(inter$.,rownames(sig_diff))%>%as.data.frame()
write.table(inter,file = 'intersect.xls',sep = '\t',row.names = F,quote = F)

library(ggvenn)
mydata<-list('DEG(scRNA)'=scDEG$gene,'DEG(GEO)'=rownames(sig_diff),'modGene'=modgene$modgene,Ferroptosis=ferr$Ferr_gene)
pdf('05.venn.pdf',w=9,h=9)
ggvenn(mydata,c('DEG(scRNA)','DEG(GEO)','modGene','Ferroptosis'),
       fill_color = c("#ffb2b2","#b2e7cb","#F0E442", "#0072B2"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       text_size = 6,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb","#F0E442", "#0072B2"),
       text_color = 'black')
dev.off()
png('05.venn.png',w=700,h=700)
ggvenn(mydata,c('DEG(scRNA)','DEG(GEO)','modGene','Ferroptosis'),
       fill_color = c("#ffb2b2","#b2e7cb","#F0E442", "#0072B2"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       text_size = 6,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb","#F0E442", "#0072B2"),
       text_color = 'black')
dev.off()


## 火山图------
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(Ipaper)
library(scales)
library(ggrepel)
# dat_rep<-DEG[rownames(DEG)%in%
#                rownames(rbind(head(sig_diff[order(sig_diff$logFC,decreasing = T),],10),
#                               head(sig_diff[order(sig_diff$logFC,decreasing = F),],10))),]
dat_rep <- DEG[rownames(DEG)%in%inter$.,]
volcano_plot<- ggplot(data = DEG, 
                      aes(x = logFC,
                          y = -log10(adj.P.Val), 
                          color =change)) +
  scale_color_manual(values = c("blue", "darkgray","red")) +
  scale_x_continuous(breaks = c(-10,-6,-3,-1,0,1,3,6,10)) +
  scale_y_continuous(trans = "log1p",
                     breaks = c(0,1,5,10,20,50, 100,200)) +
  geom_point(size = 1.2, alpha = 0.4, na.rm=T) +
  theme_bw(base_size = 12, base_family = "Times") +
  #    geom_vline(xintercept = c(-1,1),
  #               lty = 4,
  #               col = "darkgray",
  #               lwd = 0.6)+
  geom_hline(yintercept = -log10(0.05),
             lty = 4,
             col = "darkgray",
             lwd = 0.6)+
  theme(legend.position = "right",
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(face="bold",
                                   color="black",
                                   family = "Times",
                                   size=13),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face = "bold",
                                   color = "black",
                                   size = 15),
        axis.text.y = element_text(face = "bold",
                                   color = "black",
                                   size = 15),
        axis.title.x = element_text(face = "bold",
                                    color = "black",
                                    size = 15),
        axis.title.y = element_text(face = "bold",
                                    color = "black",
                                    size = 15)) +
  geom_label_repel(
    data = dat_rep,
    aes(label = rownames(dat_rep)),
    max.overlaps = 20,
    size = 3,
    box.padding = unit(0.5, "lines"),
    min.segment.length = 0,
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )+
  labs(x = "log (Fold Change)",
       y = "-log10 (adj.P.Val)")
volcano_plot
ggsave('06.volcano.png', volcano_plot,width = 8, height = 7)
ggsave('06.volcano.pdf', volcano_plot,width = 8, height = 7)
## 热图-------
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(gplots)
group_rt<-df.group%>%as.data.frame()
group_rt<-group_rt[order(group_rt$group),]
rt<-df[,group_rt$sample]
group_rt<-data.frame(group_rt$group)
colnames(group_rt)<-'group'
rownames(group_rt)<-colnames(rt)
heat<-rt[rownames(dat_rep),]
x<-log2(heat+1)
#x<-t(scale(t(heat)))
ann_colors<-list(
  group = c(control="#00CED1",RB="#F08080"))
pdf('07.heatmap.pdf',w=8,h=6)
pheatmap(mat=x,
         annotation_col = group_rt,
         color=bluered(100),
         scale = "row",
         annotation_colors = ann_colors,
         fontsize = 10,
         show_colnames = F,
         cluster_cols = F,
         cluster_rows = T,
         show_rownames = T,
         annotation_names_row = F)
dev.off()
png('07.heatmap.png',w=650,h=500)
pheatmap(mat=x,
         annotation_col = group_rt,
         color=bluered(100),
         scale = "row",
         annotation_colors = ann_colors,
         fontsize = 10,
         show_colnames = F,
         cluster_cols = F,
         cluster_rows = T,
         show_rownames = T,
         annotation_names_row = F)
dev.off()


