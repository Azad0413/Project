rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/GY0324-12/")
if (! dir.exists("./08_DEGs2")){
  dir.create("./08_DEGs2")
}
setwd("./08_DEGs2")

# A1和A2-------
library(magrittr)
library(stringr)
library(lance)
library(limma)
df1 = read.delim2("../00_rawdata/datA1.xls", row.names = 1) %>% lc.tableToNum
df1 <- df1[,c(1:3)]
df2 = read.delim2('../00_rawdata/datA2.xls',row.names = 1)%>%lc.tableToNum()
df2 <- df2[,c(1:3)]
df <- cbind(df1,df2)

df.group = data.frame(sample=colnames(df),group=c(rep('A1',3),rep('A2',3)))
table(df.group$group)
df = df[df.group$sample]
df.group$group = factor(df.group$group, levels = c("A2", "A1"))
design.mat = cbind(A2 = ifelse(df.group$group == "A2", 1, 0), 
                   A1 = ifelse(df.group$group == "A2", 0, 1))

contrast.mat = makeContrasts(contrasts="A1-A2", levels=design.mat)
fit = lmFit(df, design.mat)
fit = contrasts.fit(fit, contrast.mat)
fit = eBayes(fit)
fit = topTable(fit, coef = 1, number = Inf, adjust.method = "fdr")
#fit = fit[c(1,4,5)]
DEG=na.omit(fit)
logFC_cutoff <- 1
DEG$change = as.factor(
  ifelse(DEG$adj.P.Val <0.05 & abs(DEG$logFC) > logFC_cutoff,
         ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
sig_diff <- subset(DEG,
                   DEG$adj.P.Val < 0.05 & abs(DEG$logFC) > logFC_cutoff)

dim(DEG)
dim(sig_diff)
## 3991    6
summary(sig_diff$change)
# DOWN  NOT   UP 
# 2012    0 1979
write.table(DEG,file = "DEG_all(A1vsA2).xls",quote = F,sep = "\t",row.names = T)
write.table(sig_diff,file = "DEG_sig(A1vsA2).xls",quote = F,sep = "\t",row.names = T)
## 火山图
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
  scale_color_manual(values = c("#20B2AA", "darkgray","red")) +
  scale_x_continuous(breaks = c(-10,-6,-3,-1,0,1,3,6,10)) +
  # scale_y_continuous(trans = "log1p",
  #                   breaks = c(0,1,5,10,20,50, 100,200)) +
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
  geom_label_repel(
    data = dat_rep,
    aes(label = rownames(dat_rep)),
    max.overlaps = 20,
    size = 3,
    box.padding = unit(0.5, "lines"),
    min.segment.length = 0,
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )+
  labs(x = "log (Fold Change)",
       y = "-log10 (adj.P.Val)",title = 'A1vs.A2')
volcano_plot
ggsave('01.volcano(A1vsA2).png', volcano_plot,width = 8, height = 7)
ggsave('01.volcano(A1vsA2).pdf', volcano_plot,width = 8, height = 7)

#A2和A3-------
df3 = read.delim2('../00_rawdata/datA3.xls',row.names = 1)%>%lc.tableToNum()
df3 <- df3[,c(1:3)]
df <- cbind(df2,df3)
df.group = data.frame(sample=colnames(df),group=c(rep('A2',3),rep('A3',3)))
table(df.group$group)
df = df[df.group$sample]
df.group$group = factor(df.group$group, levels = c("A3", "A2"))
design.mat = cbind(A3 = ifelse(df.group$group == "A3", 1, 0), 
                   A2 = ifelse(df.group$group == "A3", 0, 1))

contrast.mat = makeContrasts(contrasts="A2-A3", levels=design.mat)
fit = lmFit(df, design.mat)
fit = contrasts.fit(fit, contrast.mat)
fit = eBayes(fit)
fit = topTable(fit, coef = 1, number = Inf, adjust.method = "fdr")
#fit = fit[c(1,4,5)]
DEG=na.omit(fit)
logFC_cutoff <- 1
DEG$change = as.factor(
  ifelse(DEG$adj.P.Val <0.05 & abs(DEG$logFC) > logFC_cutoff,
         ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
sig_diff <- subset(DEG,
                   DEG$adj.P.Val < 0.05 & abs(DEG$logFC) > logFC_cutoff)

dim(DEG)
dim(sig_diff)
##  2066   6
summary(sig_diff$change)
# DOWN  NOT   UP 
#  1155    0  911
write.table(DEG,file = "DEG_all(A2vs.A3).xls",quote = F,sep = "\t",row.names = T)
write.table(sig_diff,file = "DEG_sig(A2vs.A3).xls",quote = F,sep = "\t",row.names = T)
## 火山图
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
  scale_color_manual(values = c("#20B2AA", "darkgray","red")) +
  scale_x_continuous(breaks = c(-10,-6,-3,-1,0,1,3,6,10)) +
  # scale_y_continuous(trans = "log1p",
  #                   breaks = c(0,1,5,10,20,50, 100,200)) +
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
  geom_label_repel(
    data = dat_rep,
    aes(label = rownames(dat_rep)),
    max.overlaps = 20,
    size = 3,
    box.padding = unit(0.5, "lines"),
    min.segment.length = 0,
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )+
  labs(x = "log (Fold Change)",
       y = "-log10 (adj.P.Val)",title = 'A2vs.A3')
volcano_plot
ggsave('02.volcano(A2vs.A3).png', volcano_plot,width = 8, height = 7)
ggsave('02.volcano(A2vs.A3).pdf', volcano_plot,width = 8, height = 7)
