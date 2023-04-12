rm(list = ls())
# 02 差异分析---------
setwd("/data/nas1/luchunlin/project/LZZK-519-10/")
if (! dir.exists("./03_DEGs(GSE25724)")){
  dir.create("./03_DEGs(GSE25724)")
}
setwd("./03_DEGs(GSE25724)")
## 读取数据
library(magrittr)
library(stringr)
library(lance)
library(limma)
df = read.delim2("../00_rawdata/dat(GSE25724).xls", row.names = 1) %>% lc.tableToNum
#df <- log2(df)
df.group = read.delim2("../00_rawdata/group(GSE25724).xls")
df = df[df.group$sample]
table(df.group$group)
df.group$group = factor(df.group$group, levels = c("control", "T2D"))
design.mat = cbind(control = ifelse(df.group$group == "control", 1, 0), 
                   T2D = ifelse(df.group$group == "control", 0, 1))

contrast.mat = makeContrasts(contrasts="T2D-control", levels=design.mat)
fit = lmFit(df, design.mat)
fit = contrasts.fit(fit, contrast.mat)
fit = eBayes(fit)
fit = topTable(fit, coef = 1, number = Inf, adjust.method = "fdr")
#fit = fit[c(1,4,5)]
DEG=na.omit(fit)
logFC_cutoff <- 0.6
DEG$change = as.factor(
  ifelse(DEG$adj.P.Val <0.05 & abs(DEG$logFC) > logFC_cutoff,
         ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
sig_diff <- subset(DEG,
                   DEG$adj.P.Val < 0.05 & abs(DEG$logFC) > logFC_cutoff)

dim(DEG)
dim(sig_diff)
## 2327    6
summary(sig_diff$change)
# DOWN  NOT   UP 
# 1555    0  772  
write.table(DEG,file = "DEG_all(GSE25724).xls",quote = F,sep = "\t",row.names = T)
write.table(sig_diff,file = "DEG_sig(GSE25724).xls",quote = F,sep = "\t",row.names = T)

library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(Ipaper)
library(scales)
library(ggrepel)
dat_rep<-DEG[rownames(DEG)%in%
               rownames(rbind(head(sig_diff[order(sig_diff$logFC,decreasing = T),],10),
                              head(sig_diff[order(sig_diff$logFC,decreasing = F),],10))),]
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
  geom_vline(xintercept = c(-0.6,0.6),
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
ggsave('01.volcano(GSE25724).png', volcano_plot,width = 8, height = 6)
ggsave('01.volcano(GSE25724).pdf', volcano_plot,width = 8, height = 6)
### 热图--------
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
x<-log2(heat+0.1)
#x<-t(scale(t(heat)))
ann_colors<-list(
  group = c(control="#00CED1",T2D="#F08080"))
pdf('02.heatmap.pdf',w=8,h=6)
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
png('02.heatmap.png',w=650,h=500)
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


