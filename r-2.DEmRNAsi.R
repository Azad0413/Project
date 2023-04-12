rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/SJZZK-428-10/")
if (! dir.exists("./02_DEmRNAsi")){
  dir.create("./02_DEmRNAsi")
}
setwd("./02_DEmRNAsi")
dat<-read.delim2('../00_rawdata/dat.tcga.xls',row.names = 1)%>%lc.tableToNum()
colnames(dat)<-gsub('.','-',colnames(dat),fixed = T)
dat<-round(dat,digits = 0)
library(DESeq2)
colData<-read.delim2('../01_mRNAsi/group(mRNAsi).xls')
colData <- colData[order(colData$group),]
dat <- dat[,colData$sample]
table(colData$group)
colData$group<-factor(colData$group,levels = c('Low','High'))
dds<-DESeqDataSetFromMatrix(countData = dat,colData=colData,design = ~group)
dds = dds[rownames(counts(dds)) > 1,]
dds<-estimateSizeFactors(dds)
##提取标准化后的数据 
#normalized_counts <- counts(dds,normalized=T) 
#write.table(normalized_counts,file = 'normalized.counts.xls',sep = '\t',row.names = T,quote = F)
dds<-DESeq(dds)
## 提取差异结果
res =results(dds, contrast = c("group","High","Low"))
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
##   560 18749   269   
sig_diff <- subset(DEG,
                   DEG$padj < 0.05 & abs(DEG$log2FoldChange) >= logFC_cutoff)
## 829
DEG_write <- cbind(GeneSymbol=rownames(DEG), DEG)
write.table(DEG_write, file = "DEG_all(mRNAsi).xls",
            quote = F,
            sep = "\t",
            row.names = F)
sig_diff_write <- cbind(GeneSymbol=rownames(sig_diff), sig_diff)
write.table(sig_diff_write, file = "DEG_sig(mRNAsi).xls",
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
ggsave('01.volcano(mRNAsi).png', volcano_plot,width = 8, height = 6)
ggsave('01.volcano(mRNAsi).pdf', volcano_plot,width = 8, height = 6)
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
heat<-rt[rownames(dat_rep),]
x<-log2(heat+0.1)
#x<-log2(heat+1)
#x<-t(scale(t(heat)))
ann_colors<-list(
  group = c(Low="lightblue",High="darkorange"),
  Change=c(Up="#FF0000",Down="#436EEE")
)
#annotation_raw<-data.frame(row.names = rownames(heat),
#                           Change=factor(rep(c('Down','Up'),c(10,10))))
pdf('02.heatmap(mRNAsi).pdf',w=7,h=6)
pheatmap(mat=x,
         annotation_col = group_rt,
         color=bluered(100),
         #       annotation_row = annotation_raw,
         scale = "row",
         annotation_colors = ann_colors,
         fontsize = 12,
         show_colnames = F,
         cluster_cols = F,
         cluster_rows = T,
         annotation_names_row = F)
dev.off()
png('02.heatmap(mRNAsi).png',w=600,h=500)
pheatmap(mat=x,
         annotation_col = group_rt,
         color=bluered(100),
         #       annotation_row = annotation_raw,
         scale = "row",
         annotation_colors = ann_colors,
         fontsize = 12,
         show_colnames = F,
         cluster_cols = F,
         cluster_rows = T,
         annotation_names_row = F)
dev.off()
