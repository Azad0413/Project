rm(list = ls())
# 02 差异分析----------
setwd("/data/nas1/luchunlin/project/BJTC-334")
if (! dir.exists("./01_DEGs")){
  dir.create("./01_DEGs")
}
setwd("./01_DEGs")
## 读取数据
library(magrittr)
library(stringr)
library(lance)
library(limma)
df = read.delim2("/data/nas1/luchunlin/project/BJTC-334/00_rawdata/dat.final.xls", row.names = 1) %>% lc.tableToNum
df<-log2(df+1)
df.group = read.delim2("/data/nas1/luchunlin/project/BJTC-334/00_rawdata/group.xls")
df = df[df.group$sample]
df.group$group = factor(df.group$group, levels = c("Control", "CAD"))
design.mat = cbind(Control = ifelse(df.group$group == "Control", 1, 0), 
                   CAD = ifelse(df.group$group == "Control", 0, 1))

contrast.mat = makeContrasts(contrasts="CAD-Control", levels=design.mat)

fit = lmFit(df, design.mat)
fit = contrasts.fit(fit, contrast.mat)
fit = eBayes(fit)
fit = topTable(fit, coef = 1, number = Inf, lfc=log2(1),adjust.method = "fdr")
#fit = fit[c(1,4,5)]
DEG=na.omit(fit)
logFC_cutoff <- 0
DEG$change = as.factor(
  ifelse(DEG$adj.P.Val <0.05 & abs(DEG$logFC) > logFC_cutoff,
         ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
sig_diff <- subset(DEG,
                   DEG$adj.P.Val < 0.05 & abs(DEG$logFC) > logFC_cutoff)

dim(DEG)
dim(sig_diff)
## 9485
summary(sig_diff$change)
# DOWN  NOT   UP 
# 326    0 9159 
write.table(DEG,file = "DEG_all.xls",quote = F,sep = "\t",row.names = T)
write.table(sig_diff,file = "DEG_sig.xls",quote = F,sep = "\t",row.names = T)

## 火山图------
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(Ipaper)
library(scales)
library(ggrepel)
dat_rep<-DEG[rownames(DEG)%in%
               c(head(rownames(subset(sig_diff,sig_diff$logFC>0.53)),10),
                 head(rownames(subset(sig_diff,sig_diff$logFC< -1.445)),10)),]
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
     # geom_vline(xintercept = c(-1,1),
     #            lty = 4,
     #            col = "darkgray",
     #            lwd = 0.6)+
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
ggsave('volcano.png', volcano_plot,width = 8, height = 7)
ggsave('volcano.pdf', volcano_plot,width = 8, height = 7)
## 热图-------
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(gplots)
group_rt<-df.group$group%>%as.data.frame()
rt<-df
colnames(group_rt)<-'group'
rownames(group_rt)<-df.group$sample
heat<-rt[rownames(dat_rep),]
#x<-log2(heat+1)
x<-t(scale(t(heat)))
ann_colors<-list(
  group = c(Control="#4682B4",CAD="#F4A460"))

pheatmap(mat=x,
         annotation_col = group_rt,
       #  color=colorRampPalette(color.key)(50),
         scale = "row",
         annotation_colors = ann_colors,
         fontsize = 10,
         show_colnames = F,
         cluster_cols = F,
         cluster_rows = T,
         show_rownames = T,
         annotation_names_row = F)
