## 02 差异基因鉴定-----
rm(list = ls())
setwd("/data/nas1/luchunlin/project/LZZK-512/")
if (! dir.exists("./01_DEGs")){
  dir.create("./01_DEGs")
}
setwd("./01_DEGs")
library(lance)
library(tidyverse)
dat<-read.delim2('/data/nas1/luchunlin/project/LZZK-512//00_rawdata/dat.tcga.xls',row.names = 1)%>%lc.tableToNum()
colnames<-data.frame(sample=colnames(dat))
colnames$sample<-gsub('.','-',colnames$sample,fixed = T)
colnames(dat)<-colnames$sample
dat<-round(dat,digits = 0)

library(DESeq2)
colData<-data.frame(sample=colnames(dat),
                    group=c(rep('Normal',19),rep('Tumor',411)))
write.table(colData,file = 'group.xls',sep = '\t',row.names = F,quote = F)
colData$group<-factor(colData$group,levels = c('Normal','Tumor'))
dds<-DESeqDataSetFromMatrix(countData = dat,colData=colData,design = ~group)
dds = dds[rowSums(counts(dds)) > 10,]
dds<-estimateSizeFactors(dds)
##提取标准化后的数据 
#normalized_counts <- counts(dds,normalized=T) 
#write.table(normalized_counts,file = 'normalized.counts.xls',sep = '\t',row.names = T,quote = F)
dds<-DESeq(dds)
## 提取差异结果
res =results(dds, contrast = c("group","Tumor","Normal"))
res =res[order(res$padj),]
head(res)
summary(res)
table(res$padj<0.05)
DEG <- subset(res, padj < 0.05 & abs(log2FoldChange) >1 )
DEG<-as.data.frame(res)
DEG<-na.omit(DEG)
dim(DEG)
head(DEG)
## 添加change列
logFC_cutoff<-1
DEG$change=as.factor(
  ifelse(DEG$padj<0.05&abs(DEG$log2FoldChange)>logFC_cutoff,
         ifelse(DEG$log2FoldChange>logFC_cutoff,'UP','DOWN'),'NOT'))
table(DEG$change)
## DOWN   NOT    UP 
## 3213 27609  5651 
sig_diff <- subset(DEG,
                   DEG$padj < 0.05 & abs(DEG$log2FoldChange) >= logFC_cutoff)
## 8864
DEG_write <- cbind(GeneSymbol=rownames(DEG), DEG)
write.table(DEG_write, file = "DEG_all.xls",
            quote = F,
            sep = "\t",
            row.names = F)
sig_diff_write <- cbind(GeneSymbol=rownames(sig_diff), sig_diff)
write.table(sig_diff_write, file = "DEG_sig.xls",
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
dat_rep<-DEG[rownames(DEG)%in%
               c(head(rownames(subset(sig_diff,sig_diff$log2FoldChange>7.63)),10),
                 head(rownames(subset(sig_diff,sig_diff$log2FoldChange< -5.35)),10)),]
volcano_plot<- ggplot(data = DEG, 
                      aes(x = log2FoldChange,
                          y = -log10(padj), 
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
       y = "-log10 (adj.P.Val)")
volcano_plot
ggsave('volcano.png', volcano_plot,width = 9, height = 6)
ggsave('volcano.pdf', volcano_plot,width = 9, height = 6)
### 热图--------
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(gplots)


group_rt<-colData
group_rt<-group_rt[order(group_rt$group),]
rt<-dat[,group_rt$sample]
group_rt<-group_rt$group%>%as.data.frame()
colnames(group_rt)<-'group'
rownames(group_rt)<-colnames(rt)
dat_rep<-dat_rep[order(dat_rep$change),]
exprSet=log(edgeR::cpm(dat)+1)
heat <- dat[(rownames(dat_rep)),]
x<-log2(heat+0.001)
# x <- t(scale(t(heat)))
# x[x>2] <- 2
# x[x<-2] <- 2
ann_colors<-list(
  group = c(Normal="lightblue",Tumor="darkorange"),
  Change=c(Up="#FF0000",Down="#436EEE")
)
#annotation_raw<-data.frame(row.names = rownames(heat),
#                            Change=factor(rep(c('Down','Up'),c(10,10))))
pheatmap(mat=x,
         annotation_col = group_rt,
         color=bluered(100),
        # annotation_row = annotation_raw,
         scale = "row",
         annotation_colors = ann_colors,
         fontsize = 10,
         show_colnames = F,
         cluster_cols = F,
         cluster_rows = F,
         annotation_names_row = F,
         gaps_col = 19,
        cellwidth = 0.6
       #  gaps_row = 10
       )
