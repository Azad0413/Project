rm(list = ls())
# 01 获取数据集--------------
setwd("/data/nas1/luchunlin/project/BJTC-356/")
if (! dir.exists("./02_DEGs")){
  dir.create("./02_DEGs")
}
setwd("./02_DEGs")
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(Ipaper)
library(scales)
library(ggrepel)
library(lance)
library(tidyverse)
df <- read.delim2('../00_rawdata/dat.fpkm.xls',row.names = 1)%>%lc.tableToNum()
df <- log2(df+1)
df.group <- read.delim2('../01_GSVA/group.xls')
DEG <- read.delim2('../00_rawdata/core_table.xls')
colnames(DEG)
DEG <- DEG[,c(2,22:26)]%>%column_to_rownames(var = 'gene_id_alias')%>%lc.tableToNum()
DEG <- na.omit(DEG)
logFC_cutoff <- 0.5
DEG$change = as.factor(
  ifelse(DEG$diffexp_degseq_qvalue_Persistent.vs.Paroxysmal<0.05 & abs(DEG$diffexp_log2fc_Persistent.vs.Paroxysmal) > logFC_cutoff,
         ifelse(DEG$diffexp_log2fc_Persistent.vs.Paroxysmal > logFC_cutoff ,'UP','DOWN'),'NOT')
)
sig_diff <- subset(DEG,
                   DEG$diffexp_degseq_qvalue_Persistent.vs.Paroxysmal< 0.05 & abs(DEG$diffexp_log2fc_Persistent.vs.Paroxysmal) > logFC_cutoff)

colnames(sig_diff)

### logFC t p.value adj.Pvalue B
sig_diff <- sig_diff[,c(3,5,4,6)]
colnames(sig_diff) <- c('logFC','p.value','adj.Pvalue','change')

table(sig_diff$change)
# DOWN  NOT   UP 
# 1054    0  952
dim(sig_diff)
# 2006
write.table(sig_diff,file = 'DEG_sig.xls',sep = '\t',row.names = T,quote = F)
write.table(DEG,file = 'DEG_all.xls',sep = '\t',row.names = T,quote = F)
colnames(sig_diff)
dat_rep<-DEG[rownames(DEG)%in%
               rownames(rbind(head(sig_diff[order(sig_diff$logFC,decreasing = T),],10),
                              head(sig_diff[order(sig_diff$logFC,decreasing = F),],10))),]
volcano_plot<- ggplot(data = DEG, 
                      aes(x = diffexp_log2fc_Persistent.vs.Paroxysmal,
                          y = -log10(diffexp_degseq_qvalue_Persistent.vs.Paroxysmal), 
                          color =change)) +
  scale_color_manual(values = c("#20B2AA", "darkgray","red")) +
  scale_x_continuous(breaks = c(-10,-6,-3,-1,0,1,3,6,10)) +
  scale_y_continuous(trans = "log1p",
                     breaks = c(0,1,5,10,20,50, 100,200)) +
  geom_point(size = 1.2, alpha = 0.4, na.rm=T) +
  theme_bw(base_size = 12, base_family = "Times") +
  geom_vline(xintercept = c(-0.5,0.5),
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
  # geom_label_repel(
  #   data = dat_rep,
  #   aes(label = rownames(dat_rep)),
  #   max.overlaps = 20,
  #   size = 4,
  #   box.padding = unit(0.5, "lines"),
  #   min.segment.length = 0,
  #   point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )+
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
dat_rep2<-DEG[rownames(DEG)%in%
               rownames(rbind(head(sig_diff[order(sig_diff$logFC,decreasing = T),],952),
                              head(sig_diff[order(sig_diff$logFC,decreasing = F),],1054))),]
group_rt<-df.group%>%as.data.frame()
group_rt<-group_rt[order(group_rt$group),]
rt<-df[,group_rt$sample]
group_rt<-data.frame(group_rt$group)
colnames(group_rt)<-'group'
rownames(group_rt)<-colnames(rt)
heat<-rt[rownames(dat_rep2),]

x<-heat
x<-t(scale(t(heat)))
ann_colors<-list(
  group = c(Paroxysmal="#00CED1",Persistent="#F08080"))
pdf('02.heatmap.pdf',w=6,h=5)
pheatmap(mat=x,
         annotation_col = group_rt,
         color=bluered(100),
         scale = "row",
         annotation_colors = ann_colors,
         fontsize = 10,
         show_colnames = F,
         cluster_cols = F,
         cluster_rows = T,
         show_rownames = F,
         annotation_names_row = F)
dev.off()
png('02.heatmap.png',w=400,h=300)
pheatmap(mat=x,
         annotation_col = group_rt,
         color=bluered(100),
         scale = "row",
         annotation_colors = ann_colors,
         fontsize = 10,
         show_colnames = F,
         cluster_cols = F,
         cluster_rows = T,
         show_rownames = F,
         annotation_names_row = F)
dev.off()


