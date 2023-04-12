rm(list = ls())
setwd("/data/nas1/luchunlin/project/NN-0118-2/")
if (! dir.exists("./01_DEGs(TCGA)")){
  dir.create("./01_DEGs(TCGA)")
}
setwd("./01_DEGs(TCGA)")
library(tidyverse)
library(lance)
library(lance)
library(tidyverse)
library(magrittr)
library(stringr)
library(limma)
library(readxl)
##DESEQ2--------
dat<-read.delim2('../00_rawdata/dat.tcga.xls',row.names = 1)%>%lc.tableToNum()
colnames(dat)<-gsub('.','-',colnames(dat),fixed = T)
dat<-round(dat,digits = 0)
library(DESeq2)
colData<-data.frame(sample=colnames(dat),
                    group=c(rep('Normal',108),rep('Tumor',1006)))
colData$group<-factor(colData$group,levels = c('Normal','Tumor'))
write.table(colData,file = 'group.xls',sep = '\t',quote = F,row.names = F)
dds<-DESeqDataSetFromMatrix(countData = dat,colData=colData,design = ~group)
dds = dds[rownames(counts(dds)) > 1,]
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
DEG <- subset(res, padj < 0.05 & abs(log2FoldChange) >2 )
DEG<-as.data.frame(res)
# DEG<-na.omit(DEG)
dim(DEG)
head(DEG)
## 添加change列
logFC_cutoff<-2
DEG$change=as.factor(
  ifelse(DEG$padj<0.05&abs(DEG$log2FoldChange)>logFC_cutoff,
         ifelse(DEG$log2FoldChange>logFC_cutoff,'UP','DOWN'),'NOT'))
table(DEG$change)
## DOWN   NOT    UP
##   719 16933  2045 
sig_diff <- subset(DEG,
                   DEG$padj < 0.05 & abs(DEG$log2FoldChange) >= logFC_cutoff)
## 2764
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

hubgene <- c('TPX2','KIF14','KIF20A','MKI67','ANLN','TOP2A','SFTPB','SFTPC','SFTPD')

hub.diff <- DEG[hubgene,]

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
  scale_color_manual(values = c("#20B2AA", "darkgray","red")) +
  # scale_x_continuous(breaks = c(-10,-6,-3,-1,0,1,3,6,10)) +
  # scale_y_continuous(trans = "log1p",
  #                    breaks = c(0,1,5,10,20,50, 100,200)) +
  geom_point(size = 1.2, alpha = 0.4, na.rm=T) +
  theme_bw(base_size = 12, base_family = "Times") +
  geom_vline(xintercept = c(-2,2),
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
ggsave('01.volcano.png', volcano_plot,width = 8, height = 6)
ggsave('01.volcano.pdf', volcano_plot,width = 8, height = 6)
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
x<-log2(heat+1)
#x<-log2(heat+1)
#x<-t(scale(t(heat)))
ann_colors<-list(
  group = c(Normal="lightblue",Tumor="darkorange"),
  Change=c(Up="#FF0000",Down="#436EEE")
)
#annotation_raw<-data.frame(row.names = rownames(heat),
#                           Change=factor(rep(c('Down','Up'),c(10,10))))
pdf('02.heatmap.pdf',w=8,h=6)
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
png('02.heatmap.png',w=700,h=600)
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


# 
# df = read.delim2("../00_rawdata/dat.tcga.xls", row.names = 1) %>% lc.tableToNum
# colnames(df)<-gsub('.','-',colnames(df),fixed = T)
# pcg <- read.delim2('/data/nas1/luchunlin/pipeline/PCG/PCG.xls(v22)')
# df <- df[pcg$gene_name,]
# df <- na.omit(df)
# df <- edgeR::cpm(df)
# df <- log2(df+1)
# df.group = data.frame(sample=colnames(df),
#                       group=c(rep('Normal',108),rep('Tumor',1006)))
# df = df[,df.group$sample]
# table(df.group$group)
# df.group$group = factor(df.group$group, levels = c("Normal", "Tumor"))
# design.mat = cbind(Normal = ifelse(df.group$group == "Normal", 1, 0), 
#                    Tumor = ifelse(df.group$group == "Normal", 0, 1))
# contrast.mat = makeContrasts(contrasts="Tumor-Normal", levels=design.mat)
# fit = lmFit(df, design.mat)
# fit = contrasts.fit(fit, contrast.mat)
# fit = eBayes(fit)
# fit = topTable(fit, coef = 1, number = Inf, adjust.method = "fdr")
# #fit = fit[c(1,4,5)]
# DEG=na.omit(fit)
# logFC_cutoff <- 2
# DEG$change = as.factor(
#   ifelse(DEG$adj.P.Val <0.05 & abs(DEG$logFC) > logFC_cutoff,
#          ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
# )
# sig_diff <- subset(DEG,
#                    DEG$adj.P.Val< 0.05 & abs(DEG$logFC) > logFC_cutoff)
# 
# dim(DEG)
# dim(sig_diff)
# ## 900  6
# summary(sig_diff$change)
# # DOWN  NOT   UP 
# # 574    0  326 
# hubgene <- c('TPX2','KIF14','KIF20A','MKI67','ANLN','TOP2A','SFTPB','SFTPC','SFTPD')
# 
# hub.diff <- DEG[hubgene,]
# 
# write.table(DEG,file = "DEG_all.xls",quote = F,sep = "\t",row.names = T)
# write.table(sig_diff,file = "DEG_sig.xls",quote = F,sep = "\t",row.names = T)
# ## 火山图------
# library(ggplot2)
# library(ggthemes)
# library(RColorBrewer)
# library(Ipaper)
# library(scales)
# library(ggrepel)
# dat_rep<-DEG[rownames(DEG)%in%
#                rownames(rbind(head(sig_diff[order(sig_diff$logFC,decreasing = T),],10),
#                               head(sig_diff[order(sig_diff$logFC,decreasing = F),],10))),]
# 
# volcano_plot<- ggplot(data = DEG, 
#                       aes(x = logFC,
#                           y = -log10(adj.P.Val), 
#                           color =change)) +
#   scale_color_manual(values = c("#20B2AA", "darkgray","red")) +
#   # scale_x_continuous(breaks = c(-10,-6,-3,-1,0,1,3,6,10)) +
#   # scale_y_continuous(trans = "log1p",
#   #                    breaks = c(0,1,5,10,20,50, 100,200)) +
#   geom_point(size = 1.2, alpha = 0.4, na.rm=T) +
#   theme_bw(base_size = 12, base_family = "Times") +
#   geom_vline(xintercept = c(-2,2),
#              lty = 4,
#              col = "darkgray",
#              lwd = 0.6)+
#   geom_hline(yintercept = -log10(0.05),
#              lty = 4,
#              col = "darkgray",
#              lwd = 0.6)+
#   theme(legend.position = "right",
#         panel.grid = element_blank(),
#         legend.title = element_blank(),
#         legend.text = element_text(face="bold",
#                                    color="black",
#                                    family = "Times",
#                                    size=13),
#         plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_text(face = "bold",
#                                    color = "black",
#                                    size = 15),
#         axis.text.y = element_text(face = "bold",
#                                    color = "black",
#                                    size = 15),
#         axis.title.x = element_text(face = "bold",
#                                     color = "black",
#                                     size = 15),
#         axis.title.y = element_text(face = "bold",
#                                     color = "black",
#                                     size = 15)) +
#   geom_label_repel(
#     data = dat_rep,
#     aes(label = rownames(dat_rep)),
#     max.overlaps = 20,
#     size = 3,
#     box.padding = unit(0.5, "lines"),
#     min.segment.length = 0,
#     point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )+
#   labs(x = "log (Fold Change)",
#        y = "-log10 (adj.P.Val)")
# volcano_plot
# ggsave('01.volcano.png', volcano_plot,width = 8, height = 7)
# ggsave('01.volcano.pdf', volcano_plot,width = 8, height = 7)
# ## 热图-------
# library(pheatmap)
# library(ggplot2)
# library(RColorBrewer)
# library(gplots)
# group_rt<-df.group%>%as.data.frame()
# group_rt<-group_rt[order(group_rt$group),]
# rt<-df[,group_rt$sample]
# group_rt<-data.frame(group_rt$group)
# colnames(group_rt)<-'group'
# rownames(group_rt)<-colnames(rt)
# heat<-rt[rownames(dat_rep),]
# x <- heat
# #x<-log2(heat+1)
# #x<-t(scale(t(heat)))
# ann_colors<-list(
#   group = c(Normal="#20B2AA",Tumor="#F08080"))
# pdf('02.heatmap.pdf',w=8,h=6)
# pheatmap(mat=x,
#          annotation_col = group_rt,
#          color=bluered(100),
#          scale = "row",
#          annotation_colors = ann_colors,
#          fontsize = 10,
#          show_colnames = F,
#          cluster_cols = F,
#          cluster_rows = T,
#          show_rownames = T,
#          annotation_names_row = F)
# dev.off()
# png('02.heatmap.png',w=650,h=500)
# pheatmap(mat=x,
#          annotation_col = group_rt,
#          color=bluered(100),
#          scale = "row",
#          annotation_colors = ann_colors,
#          fontsize = 10,
#          show_colnames = F,
#          cluster_cols = F,
#          cluster_rows = T,
#          show_rownames = T,
#          annotation_names_row = F)
# dev.off()
