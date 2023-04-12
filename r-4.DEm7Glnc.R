rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-386-10/")
if (! dir.exists("./04_DEm7Glnc")){
  dir.create("./04_DEm7Glnc")
}
setwd("./04_DEm7Glnc")
diff <- read.delim2('../02_DEG/02.DEG_sig.xls',row.names = 1)
##1773
modgene <- read.delim2('../01_WGCNA/modGene.xls')
##575
m7Glnc <- read.delim2('../03_m7Glnc/04.correlation.final.xls')
intersect <- data.frame(symbol=intersect(rownames(diff),modgene$modgene))
##392
intersect <- data.frame(symbol=intersect(intersect$symbol,m7Glnc$lncRNA))
## 42
write.table(intersect,file = '01.intersect.xls',sep = '\t',row.names = F,quote = F)


library(ggvenn)

mydata<-list('DElncRNAs'=rownames(diff),'ModlncRNAs'=modgene$modgene,'m7GlncRNAs'=m7Glnc$lncRNA)
pdf(file = '01.intersect.pdf',w=6,h=6)
ggvenn(mydata,c('DElncRNAs','ModlncRNAs','m7GlncRNAs'),
       fill_color = c("#ffb2b2","#b2e7cb","#F0E442"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb","#F0E442"),
       text_color = 'black')
dev.off()
png(file = '01.intersect.png',w=400,h=400)
ggvenn(mydata,c('DElncRNAs','ModlncRNAs','m7GlncRNAs'),
       fill_color = c("#ffb2b2","#b2e7cb","#F0E442"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb","#F0E442"),
       text_color = 'black')
dev.off()

## 火山图------
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(Ipaper)
library(scales)
library(ggrepel)
DEG <- read.delim2('../02_DEG/01.DEG_all.xls',row.names = 1)%>%lc.tableToNum()
dat_rep<-DEG[intersect$symbol,]
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
    max.overlaps = 100,
    size = 2,
    box.padding = unit(0.5, "lines"),
    min.segment.length = 0,
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )+
  labs(x = "log (Fold Change)",
       y = "-log10 (adj.P.Val)")
volcano_plot
ggsave('02.volcano.png', volcano_plot,width = 8, height = 6)
ggsave('02.volcano.pdf', volcano_plot,width = 8, height = 6)
### 热图--------
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(gplots)
group_rt<-read.delim2('../01_WGCNA/group.xls')
group_rt<-group_rt[order(group_rt$group),]
dat<-read.delim2('../00_rawdata/dat.fpkm.xls',row.names = 1)%>%lc.tableToNum()
colnames(dat)<-gsub('.','-',colnames(dat),fixed = T)

rt<-dat[,group_rt$sample]
group_rt<-group_rt$group%>%as.data.frame()
colnames(group_rt)<-'group'
rownames(group_rt)<-colnames(rt)
dat_rep<-dat_rep[order(dat_rep$change),]
heat<-rt[rownames(dat_rep),]
x<-heat
x<-log2(heat+1)
#x<-t(scale(t(heat)))
ann_colors<-list(
  group = c(Normal="lightblue",Tumor="darkorange"),
  Change=c(Up="#FF0000",Down="#436EEE")
)
#annotation_raw<-data.frame(row.names = rownames(heat),
#                           Change=factor(rep(c('Down','Up'),c(10,10))))
pdf('03.heatmap.pdf',w=8,h=6)
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
dev.off()
png('03.heatmap.png',w=700,h=600)
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
dev.off()

