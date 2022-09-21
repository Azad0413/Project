## 02 差异基因鉴定-----
rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-258/")
if (! dir.exists("./02_DEGs")){
  dir.create("./02_DEGs")
}
setwd("./02_DEGs")
library(lance)
library(tidyverse)
dat<-read.delim2("/data/nas1/luchunlin/project/BJTC-258/00_rawdata/dat.tcga.xls", row.names = 1)%>% lc.tableToNum
colname<-data.frame(sample=colnames(dat))
colname$sample<-gsub('.','-',colname$sample,fixed = T)
colnames(dat)<-colname$sample
dat<-round(dat,digits = 0)
library(DESeq2)
cluster<-read.table('/data/nas1/luchunlin/project/BJTC-258/01_consensus/cluster.xls')
table(cluster$cluster)
cluster<-cluster[order(cluster$cluster),]
## cluster2 vs 1-------
dat1<-dat[,cluster$sample]
dat1<-dat1[,c(1:288)]
colData<-data.frame(sample=colnames(dat1),
                    group=c(rep('cluster1',146),rep('cluster2',142)))
colData$group<-factor(colData$group,levels = c('cluster1','cluster2'))
dds<-DESeqDataSetFromMatrix(countData = dat1,colData=colData,design = ~group)
dds = dds[rowSums(counts(dds)) > 1,]
dds<-estimateSizeFactors(dds)
##提取标准化后的数据 
#normalized_counts <- counts(dds,normalized=T) 
#write.table(normalized_counts,file = 'normalized.counts.xls',sep = '\t',row.names = T,quote = F)
dds<-DESeq(dds)
## 提取差异结果
res =results(dds, contrast = c("group","cluster2","cluster1"))
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
##   227 18435   117 
sig_diff <- subset(DEG,
                   DEG$padj < 0.05 & abs(DEG$log2FoldChange) >= logFC_cutoff)
## 344
DEG_write <- cbind(GeneSymbol=rownames(DEG), DEG)
write.table(DEG_write, file = "DEG_all(2vs1).xls",
            quote = F,
            sep = "\t",
            row.names = F)
sig_diff_write <- cbind(GeneSymbol=rownames(sig_diff), sig_diff)
write.table(sig_diff_write, file = "DEG_sig(2vs1).xls",
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
               rownames(rbind(head(sig_diff[order(sig_diff$log2FoldChange,decreasing = T),],10),
                              head(sig_diff[order(sig_diff$log2FoldChange,decreasing = F),],10))),]
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
ggsave('01.volcano(2v1).png', volcano_plot,width = 8, height = 6)
ggsave('01.volcano(2v1).pdf', volcano_plot,width = 8, height = 6)
## cluster3 vs 1-------
dat2<-dat[,cluster$sample]
dat2<-dat2[which(colnames(dat2)%in%cluster$sample[which(cluster$cluster=='cluster1'|cluster$cluster=='cluster3')])]
table(cluster$cluster)
colData2<-data.frame(sample=colnames(dat2),
                    group=c(rep('cluster1',146),rep('cluster3',90)))
colData2$group<-factor(colData2$group,levels = c('cluster1','cluster3'))
dds2<-DESeqDataSetFromMatrix(countData = dat2,colData=colData2,design = ~group)
dds2 = dds2[rowSums(counts(dds2)) > 1,]
dds2<-estimateSizeFactors(dds2)
##提取标准化后的数据 
#normalized_counts <- counts(dds,normalized=T) 
#write.table(normalized_counts,file = 'normalized.counts.xls',sep = '\t',row.names = T,quote = F)
dds2<-DESeq(dds2)
## 提取差异结果
res2 =results(dds2, contrast = c("group","cluster3","cluster1"))
res2 =res2[order(res2$padj),]
DEG2 <- subset(res2, padj < 0.05 & abs(log2FoldChange) >1 )
DEG2<-as.data.frame(res2)
DEG2<-na.omit(DEG2)
dim(DEG2)
head(DEG2)
## 添加change列
logFC_cutoff<-1
DEG2$change=as.factor(
  ifelse(DEG2$padj<0.05&abs(DEG2$log2FoldChange)>logFC_cutoff,
         ifelse(DEG2$log2FoldChange>logFC_cutoff,'UP','DOWN'),'NOT'))
table(DEG2$change)
## DOWN   NOT    UP 
##   455 18568   239 
sig_diff2 <- subset(DEG2,
                   DEG2$padj < 0.05 & abs(DEG2$log2FoldChange) >= logFC_cutoff)
## 694
DEG2_write <- cbind(GeneSymbol=rownames(DEG2), DEG2)
write.table(DEG2_write, file = "DEG_all(3vs1).xls",
            quote = F,
            sep = "\t",
            row.names = F)
sig_diff2_write <- cbind(GeneSymbol=rownames(sig_diff2), sig_diff2)
write.table(sig_diff2_write, file = "DEG_sig(3vs1).xls",
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
dat_rep<-DEG2[rownames(DEG2)%in%
               rownames(rbind(head(sig_diff2[order(sig_diff2$log2FoldChange,decreasing = T),],10),
                              head(sig_diff2[order(sig_diff2$log2FoldChange,decreasing = F),],10))),]
volcano_plot<- ggplot(data = DEG2, 
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
ggsave('02.volcano(3v1).png', volcano_plot,width = 8, height = 6)
ggsave('02.volcano(3v1).pdf', volcano_plot,width = 8, height = 6)

## cluster3 vs 2-------
dat3<-dat[,cluster$sample]
dat3<-dat3[which(colnames(dat3)%in%cluster$sample[which(cluster$cluster=='cluster2'|cluster$cluster=='cluster3')])]
table(cluster$cluster)
colData3<-data.frame(sample=colnames(dat3),
                     group=c(rep('cluster2',142),rep('cluster3',90)))
colData3$group<-factor(colData3$group,levels = c('cluster2','cluster3'))
dds3<-DESeqDataSetFromMatrix(countData = dat3,colData=colData3,design = ~group)
dds3 = dds3[rowSums(counts(dds3)) > 1,]
dds3<-estimateSizeFactors(dds3)
##提取标准化后的数据 
#normalized_counts <- counts(dds,normalized=T) 
#write.table(normalized_counts,file = 'normalized.counts.xls',sep = '\t',row.names = T,quote = F)
dds3<-DESeq(dds3)
## 提取差异结果
res3 =results(dds3, contrast = c("group","cluster3","cluster2"))
res3 =res3[order(res3$padj),]
DEG3 <- subset(res3, padj < 0.05 & abs(log2FoldChange) >1 )
DEG3<-as.data.frame(res3)
DEG3<-na.omit(DEG3)
dim(DEG3)
head(DEG3)
## 添加change列
logFC_cutoff<-1
DEG3$change=as.factor(
  ifelse(DEG3$padj<0.05&abs(DEG3$log2FoldChange)>logFC_cutoff,
         ifelse(DEG3$log2FoldChange>logFC_cutoff,'UP','DOWN'),'NOT'))
table(DEG3$change)
## DOWN   NOT    UP 
##  244 18362   160 
sig_diff3 <- subset(DEG3,
                    DEG3$padj < 0.05 & abs(DEG3$log2FoldChange) >= logFC_cutoff)
## 404
DEG3_write <- cbind(GeneSymbol=rownames(DEG3), DEG3)
write.table(DEG3_write, file = "DEG_all(3vs2).xls",
            quote = F,
            sep = "\t",
            row.names = F)
sig_diff3_write <- cbind(GeneSymbol=rownames(sig_diff3), sig_diff3)
write.table(sig_diff3_write, file = "DEG_sig(3vs2).xls",
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
dat_rep<-DEG3[rownames(DEG3)%in%
                rownames(rbind(head(sig_diff3[order(sig_diff3$log2FoldChange,decreasing = T),],10),
                               head(sig_diff3[order(sig_diff3$log2FoldChange,decreasing = F),],10))),]
volcano_plot<- ggplot(data = DEG3, 
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
ggsave('03.volcano(3v2).png', volcano_plot,width = 8, height = 6)
ggsave('03.volcano(3v2).pdf', volcano_plot,width = 8, height = 6)

### 合并------
sig.all<-data.frame(symbol=c(cbind(symbol=cbind(rownames(sig_diff),cbind(rownames(sig_diff2),rownames(sig_diff3))))))
length(unique(sig.all$symbol))
sig.all<-sig.all[!duplicated(sig.all$symbol),]%>%as.data.frame()
#1041
colnames(sig.all)<-'symbol'
write.table(sig.all,file = 'sig.all.xls',sep = '\t',row.names = F,quote = F)
