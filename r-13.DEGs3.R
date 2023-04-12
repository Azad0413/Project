rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/GY0324-12/")
if (! dir.exists("./13_DEGs3")){
  dir.create("./13_DEGs3")
}
setwd("./13_DEGs3")

# A3和B4-------
library(magrittr)
library(stringr)
library(lance)
library(limma)
df3 = read.delim2("../00_rawdata/datA3.xls", row.names = 1) %>% lc.tableToNum
df3 <- df3[,c(1:3)]
df4 = read.delim2('../00_rawdata/datB4.xls',row.names = 1)%>%lc.tableToNum()
df4 <- df4[,c(1:3)]
df <- cbind(df3,df4)

df.group = data.frame(sample=colnames(df),group=c(rep('A3',3),rep('B4',3)))
table(df.group$group)
df = df[df.group$sample]
df.group$group = factor(df.group$group, levels = c("B4", "A3"))
design.mat = cbind(B4 = ifelse(df.group$group == "B4", 1, 0), 
                   A3 = ifelse(df.group$group == "B4", 0, 1))

contrast.mat = makeContrasts(contrasts="A3-B4", levels=design.mat)
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
## 2337    6
summary(sig_diff$change)
# DOWN  NOT   UP 
# 726    0 1611 
write.table(DEG,file = "DEG_all(A3vsB4).xls",quote = F,sep = "\t",row.names = T)
write.table(sig_diff,file = "DEG_sig(A3vsB4).xls",quote = F,sep = "\t",row.names = T)
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
       y = "-log10 (adj.P.Val)",title = 'A3vs.B4')
volcano_plot
ggsave('01.volcano(A3vsB4).png', volcano_plot,width = 8, height = 7)
ggsave('01.volcano(A3vsB4).pdf', volcano_plot,width = 8, height = 7)
