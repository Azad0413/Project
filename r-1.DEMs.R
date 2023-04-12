rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/GY0324-12/")
if (! dir.exists("./01_DEMs")){
  dir.create("./01_DEMs")
}
setwd("./01_DEMs")
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(Ipaper)
library(scales)
library(ggrepel)
library(lance)
library(tidyverse)
## A1 vs C5------
DEM <- read.csv('t_test(A1vs.C5).csv',row.names = 1)

colnames(DEM)
logFC_cutoff <- 1
DEM$change = as.factor(
  ifelse(DEM$p.value<0.05 & abs(DEM$log2.FC.) > logFC_cutoff &DEM$VIP>1,
         ifelse(DEM$log2.FC. > logFC_cutoff ,'UP','DOWN'),'NOT')
)
sig_diff <- subset(DEM,
                   DEM$p.value< 0.05 & abs(DEM$log2.FC.)> logFC_cutoff  & DEM$VIP>1)

colnames(sig_diff)

colnames(sig_diff) <- c("t.stat","p.value","-log10(p).","FDR","Fold.Change","log2.FC.","VIP","change")

table(sig_diff$change)
# DOWN  NOT   UP 
#  25    0   17
dim(sig_diff)
# 42
write.table(sig_diff,file = 'DEM_sig(A1vs.C5).xls',sep = '\t',row.names = T,quote = F)
# write.table(DEM,file = 'DEM_all(A1vs.C5).xls',sep = '\t',row.names = T,quote = F)
colnames(sig_diff)
dat_rep<-DEM[rownames(DEM)%in%
               rownames(rbind(head(sig_diff[order(sig_diff$log2.FC.,decreasing = T),],10),
                              head(sig_diff[order(sig_diff$log2.FC.,decreasing = F),],10))),]
volcano_plot<- ggplot(data = DEM, 
                      aes(x = log2.FC.,
                          y = -log10(p.value), 
                          color =change)) +
  scale_color_manual(values = c("#20B2AA", "darkgray","red")) +
  # scale_x_continuous(breaks = c(-10,-6,-3,-1,0,1,3,6,10)) +
  # scale_y_continuous(trans = "log1p",
  #                    breaks = c(0,1,5,10,20,50, 100,200)) +
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
        plot.title = element_text(hjust = 0.5,face = 'bold'),
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
                                    size = 15),) +
  # geom_label_repel(
  #   data = dat_rep,
  #   aes(label = rownames(dat_rep)),
  #   max.overlaps = 20,
  #   size = 4,
  #   box.padding = unit(0.5, "lines"),
  #   min.segment.length = 0,
  #   point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )+
  labs(x = "log (Fold Change)",
       y = "-log10 (p.value)",title = 'A1vs.C5')
volcano_plot
ggsave('08.volcano(A1vs.C5).png', volcano_plot,width = 6, height = 4)
ggsave('08.volcano(A1vs.C5).pdf', volcano_plot,width = 6, height = 4)


## A2 vs C5------
DEM <- read.csv('t_test(A2vs.C5).csv',row.names = 1)

colnames(DEM)
logFC_cutoff <- 1
DEM$change = as.factor(
  ifelse(DEM$p.value<0.05 & abs(DEM$log2.FC.) > logFC_cutoff &DEM$VIP>1,
         ifelse(DEM$log2.FC. > logFC_cutoff ,'UP','DOWN'),'NOT')
)
sig_diff <- subset(DEM,
                   DEM$p.value< 0.05 & abs(DEM$log2.FC.)> logFC_cutoff  & DEM$VIP>1)

colnames(sig_diff)

colnames(sig_diff) <- c("t.stat","p.value","-log10(p).","FDR","Fold.Change","log2.FC.","VIP","change")

table(sig_diff$change)
# DOWN  NOT   UP 
#   41    0   17 
dim(sig_diff)
# 58
write.table(sig_diff,file = 'DEM_sig(A2vs.C5).xls',sep = '\t',row.names = T,quote = F)
# write.table(DEM,file = 'DEM_all(A2vs.C5).xls',sep = '\t',row.names = T,quote = F)
colnames(sig_diff)
dat_rep<-DEM[rownames(DEM)%in%
               rownames(rbind(head(sig_diff[order(sig_diff$log2.FC.,decreasing = T),],10),
                              head(sig_diff[order(sig_diff$log2.FC.,decreasing = F),],10))),]
volcano_plot<- ggplot(data = DEM, 
                      aes(x = log2.FC.,
                          y = -log10(p.value), 
                          color =change)) +
  scale_color_manual(values = c("#20B2AA", "darkgray","red")) +
  # scale_x_continuous(breaks = c(-10,-6,-3,-1,0,1,3,6,10)) +
  # scale_y_continuous(trans = "log1p",
  #                    breaks = c(0,1,5,10,20,50, 100,200)) +
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
        plot.title = element_text(hjust = 0.5,face = 'bold'),
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
       y = "-log10 (p.value)",title = 'A2vs.C5')
volcano_plot
ggsave('09.volcano(A2vs.C5).png', volcano_plot,width = 6, height = 4)
ggsave('09.volcano(A2vs.C5).pdf', volcano_plot,width = 6, height = 4)


## A3 vs C5------
DEM <- read.csv('t_test(A3vs.C5).csv',row.names = 1)

colnames(DEM)
logFC_cutoff <- 1
DEM$change = as.factor(
  ifelse(DEM$p.value<0.05 & abs(DEM$log2.FC.) > logFC_cutoff &DEM$VIP>1,
         ifelse(DEM$log2.FC. > logFC_cutoff ,'UP','DOWN'),'NOT')
)
sig_diff <- subset(DEM,
                   DEM$p.value< 0.05 & abs(DEM$log2.FC.)> logFC_cutoff  & DEM$VIP>1)

colnames(sig_diff)

colnames(sig_diff) <- c("t.stat","p.value","-log10(p).","FDR","Fold.Change","log2.FC.","VIP","change")

table(sig_diff$change)
# DOWN  NOT   UP 
#   7    0   14 
dim(sig_diff)
# 21
write.table(sig_diff,file = 'DEM_sig(A3vs.C5).xls',sep = '\t',row.names = T,quote = F)
# write.table(DEM,file = 'DEM_all(A3vs.C5).xls',sep = '\t',row.names = T,quote = F)
colnames(sig_diff)
dat_rep<-DEM[rownames(DEM)%in%
               rownames(rbind(head(sig_diff[order(sig_diff$log2.FC.,decreasing = T),],10),
                              head(sig_diff[order(sig_diff$log2.FC.,decreasing = F),],10))),]
volcano_plot<- ggplot(data = DEM, 
                      aes(x = log2.FC.,
                          y = -log10(p.value), 
                          color =change)) +
  scale_color_manual(values = c("#20B2AA", "darkgray","red")) +
  # scale_x_continuous(breaks = c(-10,-6,-3,-1,0,1,3,6,10)) +
  # scale_y_continuous(trans = "log1p",
  #                    breaks = c(0,1,5,10,20,50, 100,200)) +
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
        plot.title = element_text(hjust = 0.5,face = 'bold'),
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
       y = "-log10 (p.value)",title = 'A3vs.C5')
volcano_plot
ggsave('10.volcano(A3vs.C5).png', volcano_plot,width = 6, height = 4)
ggsave('10.volcano(A3vs.C5).pdf', volcano_plot,width = 6, height = 4)


## B4 vs C5------
DEM <- read.csv('t_test(B4vs.C5).csv',row.names = 1)

colnames(DEM)
logFC_cutoff <- 1
DEM$change = as.factor(
  ifelse(DEM$p.value<0.05 & abs(DEM$log2.FC.) > logFC_cutoff &DEM$VIP>1,
         ifelse(DEM$log2.FC. > logFC_cutoff ,'UP','DOWN'),'NOT')
)
sig_diff <- subset(DEM,
                   DEM$p.value< 0.05 & abs(DEM$log2.FC.)> logFC_cutoff  & DEM$VIP>1)

colnames(sig_diff)

colnames(sig_diff) <- c("t.stat","p.value","-log10(p).","FDR","Fold.Change","log2.FC.","VIP","change")

table(sig_diff$change)
# DOWN  NOT   UP 
#  40    0    4 
dim(sig_diff)
# 44
write.table(sig_diff,file = 'DEM_sig(B4vs.C5).xls',sep = '\t',row.names = T,quote = F)
# write.table(DEM,file = 'DEM_all(B4vs.C5).xls',sep = '\t',row.names = T,quote = F)
colnames(sig_diff)
dat_rep<-DEM[rownames(DEM)%in%
               rownames(rbind(head(sig_diff[order(sig_diff$log2.FC.,decreasing = T),],10),
                              head(sig_diff[order(sig_diff$log2.FC.,decreasing = F),],10))),]
volcano_plot<- ggplot(data = DEM, 
                      aes(x = log2.FC.,
                          y = -log10(p.value), 
                          color =change)) +
  scale_color_manual(values = c("#20B2AA", "darkgray","red")) +
  # scale_x_continuous(breaks = c(-10,-6,-3,-1,0,1,3,6,10)) +
  # scale_y_continuous(trans = "log1p",
  #                    breaks = c(0,1,5,10,20,50, 100,200)) +
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
        plot.title = element_text(hjust = 0.5,face = 'bold'),
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
       y = "-log10 (p.value)",title = 'B4vs.C5')
volcano_plot
ggsave('11.volcano(B4vs.C5).png', volcano_plot,width = 6, height = 4)
ggsave('11.volcano(B4vs.C5).pdf', volcano_plot,width = 6, height = 4)
