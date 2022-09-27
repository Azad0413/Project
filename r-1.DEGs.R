rm(list = ls())
# 02 差异分析---------
setwd("/data/nas1/luchunlin/project/JNZK-218-8/")
if (! dir.exists("./01_DEGs")){
  dir.create("./01_DEGs")
}
setwd("./01_DEGs")
## 读取数据
library(magrittr)
library(stringr)
library(lance)
library(limma)
df = read.delim2("../00_rawdata/dat(GSE27276).xls", row.names = 1) %>% lc.tableToNum
df.group = read.delim2("../00_rawdata/group(GSE27276).xls")
df = df[df.group$sample]
df.group$group = factor(df.group$group, levels = c("control", "POAG"))
design.mat = cbind(control = ifelse(df.group$group == "control", 1, 0), 
                   POAG = ifelse(df.group$group == "control", 0, 1))

contrast.mat = makeContrasts(contrasts="POAG-control", levels=design.mat)
fit = lmFit(df, design.mat)
fit = contrasts.fit(fit, contrast.mat)
fit = eBayes(fit)
fit = topTable(fit, coef = 1, number = Inf, adjust.method = "fdr")
#fit = fit[c(1,4,5)]
DEG=na.omit(fit)
logFC_cutoff <- 0.5
DEG$change = as.factor(
  ifelse(DEG$adj.P.Val <0.05 & abs(DEG$logFC) > logFC_cutoff,
         ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
sig_diff <- subset(DEG,
                   DEG$adj.P.Val < 0.05 & abs(DEG$logFC) > logFC_cutoff)

dim(DEG)
dim(sig_diff)
## 629 
summary(sig_diff$change)
# DOWN  NOT   UP 
#  321    0  308 
write.table(DEG,file = "DEG_all(GSE27276).xls",quote = F,sep = "\t",row.names = T)
write.table(sig_diff,file = "DEG_sig(GSE27276).xls",quote = F,sep = "\t",row.names = T)
## 火山图------
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
ggsave('01.volcano(GSE27276).png', volcano_plot,width = 8, height = 7)
ggsave('01.volcano(GSE27276).pdf', volcano_plot,width = 8, height = 7)



### AS-----
df2 = read.delim2("../00_rawdata/dat(GSE100927).xls", row.names = 1) %>% lc.tableToNum
df.group2 = read.delim2("../00_rawdata/group(GSE100927).xls")
df2 = df2[df.group2$sample]
df.group2$group = factor(df.group2$group, levels = c("control", "AS"))
design.mat2 = cbind(control = ifelse(df.group2$group == "control", 1, 0), 
                   AS = ifelse(df.group2$group == "control", 0, 1))

contrast.mat2 = makeContrasts(contrasts="AS-control", levels=design.mat2)
fit2 = lmFit(df2, design.mat2)
fit2 = contrasts.fit(fit2, contrast.mat2)
fit2 = eBayes(fit2)
fit2 = topTable(fit2, coef = 1, number = Inf, adjust.method = "fdr")
#fit = fit[c(1,4,5)]
DEG2=na.omit(fit2)
logFC_cutoff <- 0.5
DEG2$change = as.factor(
  ifelse(DEG2$adj.P.Val <0.05 & abs(DEG2$logFC) > logFC_cutoff,
         ifelse(DEG2$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
sig_diff2 <- subset(DEG2,
                   DEG2$adj.P.Val < 0.05 & abs(DEG2$logFC) > logFC_cutoff)

dim(DEG2)
dim(sig_diff2)
## 3541 
summary(sig_diff2$change)
# DOWN  NOT   UP 
#1573    0 1968 
write.table(DEG,file = "DEG_all(GSE100927).xls",quote = F,sep = "\t",row.names = T)
write.table(sig_diff,file = "DEG_sig(GSE100927).xls",quote = F,sep = "\t",row.names = T)
## 火山图------
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(Ipaper)
library(scales)
library(ggrepel)
dat_rep<-DEG2[rownames(DEG2)%in%
               rownames(rbind(head(sig_diff2[order(sig_diff2$logFC,decreasing = T),],10),
                              head(sig_diff2[order(sig_diff2$logFC,decreasing = F),],10))),]
volcano_plot<- ggplot(data = DEG2, 
                      aes(x = logFC,
                          y = -log10(adj.P.Val), 
                          color =change)) +
  scale_color_manual(values = c("blue", "darkgray","red")) +
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
ggsave('02.volcano(GSE100927).png', volcano_plot,width = 8, height = 7)
ggsave('02.volcano(GSE100927).pdf', volcano_plot,width = 8, height = 7)

### POAG-AS
intersection<-rownames(sig_diff)[rownames(sig_diff)%in%rownames(sig_diff2)]%>%as.data.frame()
write.table(intersection,file = 'interDEGs.xls',sep = '\t',row.names = F,quote = F)
library(ggvenn)
mydata<-list('DEGs(POAG)'=rownames(sig_diff),'DEGs(AS)'=rownames(sig_diff2))
pdf('03.DEGs_venn.pdf',w=5,h=5)
ggvenn(mydata,c('DEGs(POAG)','DEGs(AS)'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()

png('03.DEGs_venn.png',w=400,h=400)
ggvenn(mydata,c('DEGs(POAG)','DEGs(AS)'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()

