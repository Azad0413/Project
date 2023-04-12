## 02 差异基因鉴定-----
rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-327/")
if (! dir.exists("./01_DEGs")){
  dir.create("./01_DEGs")
}
setwd("./01_DEGs")
library(lance)
library(tidyverse)
dat1<-read.delim2('/data/nas1/luchunlin/project/BJTC-327/00_rawdata/AIH.PBC.count.xls',row.names = 1)%>%lc.tableToNum()
dat2<-read.delim2('/data/nas1/luchunlin/project/BJTC-327/00_rawdata/AIH.HBV.count.xls',row.names = 1)%>%lc.tableToNum()
dat1<-na.omit(dat1)
dat2<-na.omit(dat2)
## 02-1 AIH.PBC------

library(DESeq2)
colData<-data.frame(sample=colnames(dat1),group=c(rep('AIH',3),rep('PBC',3)))
write.table(colData,file = 'group1.xls',sep = '\t',row.names = F,quote = F)
colData$group<-factor(colData$group,levels = c('AIH','PBC'))
dds<-DESeqDataSetFromMatrix(countData = dat1,colData=colData,design = ~group)
dds = dds[rownames(counts(dds)) > 1,]
dds<-estimateSizeFactors(dds)
##提取标准化后的数据 
#normalized_counts <- counts(dds,normalized=T) 
dds<-DESeq(dds)
#write.table(normalized_counts,file = 'normalized.counts.xls',sep = '\t',row.names = T,quote = F)
## 提取差异结果
## deseq2默认的是后一个比前一个
res =results(dds, contrast = c("group","PBC","AIH"))
res =res[order(res$padj),]
head(res)
summary(res)
table(res$pvalue<0.05)
DEG1 <- subset(res, pvalue < 0.05 & abs(log2FoldChange) >0 )
DEG1<-as.data.frame(res)
DEG1<-na.omit(DEG1)
dim(DEG1)
head(DEG1)
## 添加change列
logFC_cutoff<-1
DEG1$change=as.factor(
  ifelse(DEG1$pvalue<0.05&abs(DEG1$log2FoldChange)>logFC_cutoff,
         ifelse(DEG1$log2FoldChange>logFC_cutoff,'UP','DOWN'),'NOT'))
table(DEG1$change)
## DOWN   NOT    UP 
##  237 21407   184 
sig_diff1 <- subset(DEG1,
                   DEG1$pvalue < 0.05 & abs(DEG1$log2FoldChange) >= logFC_cutoff)
## 421
DEG1_write <- cbind(GeneSymbol=rownames(DEG1), DEG1)
write.table(DEG1_write, file = "DEG_all(AIH vs.PBC).xls",
            quote = F,
            sep = "\t",
            row.names = F)
sig_diff_write <- cbind(GeneSymbol=rownames(sig_diff1), sig_diff1)
write.table(sig_diff_write, file = "DEG_sig(AIH vs.PBC).xls",
            quote = F,
            sep = "\t",
            row.names = F)
### 火山图---------
#devtools::install_github("kongdd/Ipaper")
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(Ipaper)
library(scales)
library(ggrepel)
dat_rep<-DEG1[rownames(DEG1)%in%
               c(head(rownames(subset(sig_diff1,sig_diff1$log2FoldChange>4.9)),10),
                 head(rownames(subset(sig_diff1,sig_diff1$log2FoldChange< -5.2)),10)),]
volcano_plot<- ggplot(data = DEG1, 
                      aes(x = log2FoldChange,
                          y = -log10(pvalue), 
                          color =change)) +
  scale_color_manual(values = c("blue", "darkgray","red")) +
  scale_x_continuous(breaks = c(-10,-6,-3,-1,0,1,3,6,10)) +
  scale_y_continuous(trans = "log1p",
                     breaks = c(0,1,5,10,20,50, 100,200)) +
  geom_point(size = 1.2, alpha = 0.4, na.rm=T) +
  theme_bw(base_size = 12, base_family = "Times") +
  geom_vline(xintercept = c(-1,1),
             lty = 4,
             col = "darkgray",
             lwd = 0.6)+
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
    size = 4,
    box.padding = unit(0.5, "lines"),
    min.segment.length = 0,
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )+
  labs(x = "log (Fold Change)",
       y = "-log10 (P.Val)")
volcano_plot
ggsave('volcano.png', volcano_plot,width = 8, height = 7)
ggsave('volcano.pdf', volcano_plot,width = 8, height = 7)
### 热图--------
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(gplots)
group_rt<-colData
group_rt<-group_rt[order(group_rt$group),]
rt<-dat1[,group_rt$sample]
group_rt<-group_rt$group%>%as.data.frame()
colnames(group_rt)<-'group'
rownames(group_rt)<-colnames(rt)
dat_rep<-dat_rep[order(dat_rep$change),]
heat<-rt[rownames(dat_rep),]
x<-log2(heat+0.01)
#x<-log2(heat+1)
#x<-t(scale(t(heat)))
ann_colors<-list(
  group = c(AIH="lightblue",PBC="darkorange"),
  Change=c(Up="#FF0000",Down="#436EEE")
)
#annotation_raw<-data.frame(row.names = rownames(heat),
#                           Change=factor(rep(c('Down','Up'),c(10,10))))
pheatmap(mat=x,
         annotation_col = group_rt,
         color=bluered(100),
         #       annotation_row = annotation_raw,
         scale = "row",
         annotation_colors = ann_colors,
         fontsize = 10,
         show_colnames = F,
         cluster_cols = F,
         cluster_rows = T,
         annotation_names_row = F)
## 02-1 AIH.HBV------

colData<-data.frame(sample=colnames(dat2),group=c(rep('AIH',3),rep('HBV',3)))
write.table(colData,file = 'group2.xls',sep = '\t',row.names = F,quote = F)
colData$group<-factor(colData$group,levels = c('AIH','HBV'))
dds<-DESeqDataSetFromMatrix(countData = dat2,colData=colData,design = ~group)
dds = dds[rownames(counts(dds)) > 1,]
dds<-estimateSizeFactors(dds)
##提取标准化后的数据 
#normalized_counts <- counts(dds,normalized=T) 
dds<-DESeq(dds)
#write.table(normalized_counts,file = 'normalized.counts.xls',sep = '\t',row.names = T,quote = F)
## 提取差异结果
## deseq2默认的是后一个比前一个
res =results(dds, contrast = c("group","HBV","AIH"))
res =res[order(res$padj),]
head(res)
summary(res)
table(res$pvalue<0.05)
DEG2 <- subset(res, pvalue < 0.05 & abs(log2FoldChange) >0 )
DEG2<-as.data.frame(res)
DEG2<-na.omit(DEG2)
dim(DEG2)
head(DEG2)
## 添加change列
logFC_cutoff<-1
DEG2$change=as.factor(
  ifelse(DEG2$pvalue<0.05&abs(DEG2$log2FoldChange)>logFC_cutoff,
         ifelse(DEG2$log2FoldChange>logFC_cutoff,'UP','DOWN'),'NOT'))
table(DEG2$change)
## DOWN   NOT    UP 
##  388 19865   397 

sig_diff2 <- subset(DEG2,
                   DEG2$pvalue < 0.05 & abs(DEG2$log2FoldChange) >= logFC_cutoff)
## 785
DEG2_write <- cbind(GeneSymbol=rownames(DEG2), DEG2)
write.table(DEG2_write, file = "DEG_all(AIH vs.HBV).xls",
            quote = F,
            sep = "\t",
            row.names = F)
sig_diff_write <- cbind(GeneSymbol=rownames(sig_diff2), sig_diff2)
write.table(sig_diff_write, file = "DEG_sig(AIH vs.HBV).xls",
            quote = F,
            sep = "\t",
            row.names = F)


dat_rep<-DEG2[rownames(DEG2)%in%
                c(head(rownames(subset(sig_diff2,sig_diff2$log2FoldChange>8)),10),
                  head(rownames(subset(sig_diff2,sig_diff2$log2FoldChange< -5.45)),10)),]
volcano_plot<- ggplot(data = DEG2, 
                      aes(x = log2FoldChange,
                          y = -log10(pvalue), 
                          color =change)) +
  scale_color_manual(values = c("blue", "darkgray","red")) +
  scale_x_continuous(breaks = c(-10,-6,-3,-1,0,1,3,6,10)) +
  scale_y_continuous(trans = "log1p",
                     breaks = c(0,1,5,10,20,50, 100,200)) +
  geom_point(size = 1.2, alpha = 0.4, na.rm=T) +
  theme_bw(base_size = 12, base_family = "Times") +
  geom_vline(xintercept = c(-1,1),
             lty = 4,
             col = "darkgray",
             lwd = 0.6)+
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
    size = 4,
    box.padding = unit(0.5, "lines"),
    min.segment.length = 0,
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )+
  labs(x = "log (Fold Change)",
       y = "-log10 (P.Val)")
volcano_plot
ggsave('volcano2.png', volcano_plot,width = 8, height = 7)
ggsave('volcano2.pdf', volcano_plot,width = 8, height = 7)
### 热图--------
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(gplots)
group_rt<-colData
group_rt<-group_rt[order(group_rt$group),]
rt<-dat2[,group_rt$sample]
group_rt<-group_rt$group%>%as.data.frame()
colnames(group_rt)<-'group'
rownames(group_rt)<-colnames(rt)
dat_rep<-dat_rep[order(dat_rep$change),]
heat<-rt[rownames(dat_rep),]
x<-log2(heat+0.01)
#x<-log2(heat+1)
#x<-t(scale(t(heat)))
ann_colors<-list(
  group = c(AIH="lightblue",HBV="darkorange"),
  Change=c(Up="#FF0000",Down="#436EEE")
)
#annotation_raw<-data.frame(row.names = rownames(heat),
#                           Change=factor(rep(c('Down','Up'),c(10,10))))
pheatmap(mat=x,
         annotation_col = group_rt,
         color=bluered(100),
         #       annotation_row = annotation_raw,
         scale = "row",
         annotation_colors = ann_colors,
         fontsize = 10,
         show_colnames = F,
         cluster_cols = F,
         cluster_rows = T,
         annotation_names_row = F)
