rm(list = ls())
#08 estimate-------------
setwd("/data/nas1/luchunlin/project/LZZK-512/")
if (! dir.exists("./09_estimate")){
  dir.create("./09_estimate")
}
setwd("./09_estimate")
library(tidyverse)

dat.tcga<-read.delim2("/data/nas1/luchunlin/project/LZZK-512/00_rawdata/dat.fpkm.xls", row.names = 1)%>% lc.tableToNum
colname<-data.frame(sample=colnames(dat.tcga))
colname$sample<-gsub('.','-',colname$sample,fixed = T)
colnames(dat.tcga)<-colname$sample

library(estimate)
expr_train <- log2(dat.tcga+1)
write.table(expr_train, 
            'expr_sample408_log2.txt', 
            col.names = T, 
            row.names = T, 
            quote = F, sep="\t")

# 生成expr_train.gct
filterCommonGenes(input.f = './expr_sample408_log2.txt', 
                  output.f = 'expr_train.gct', 
                  id = 'GeneSymbol')
# [1] "Merged dataset includes 10221 genes (191 mismatched)."
# 生成train_purity.gct
estimateScore('expr_train.gct', 'train_purity.gct', platform="affymetrix")
# [1] "1 gene set: StromalSignature  overlap= 139"
# [1] "2 gene set: ImmuneSignature  overlap= 141"
es_score <- read.table('train_purity.gct', skip = 2, header = T, check.names = F)
immu_score <- es_score[,3:length(es_score)]
rownames(immu_score) <- es_score$NAME

write.table(es_score,
            file = "es_score.xls",
            sep = "\t",
            quote = F,
            row.names = F)

#violin_dat <- data.frame(t(immu_score))
#violin_dat$sample <- rownames(violin_dat)
risk<-read.delim2('/data/nas1/luchunlin/project/LZZK-512/06_risk/risk.xls')
high.sample<-risk$id[which(risk$risk==0)]
low.sample<-risk$id[which(risk$risk==1)]
#violin_dat$group <- ifelse(violin_dat$sample %in% gsub("-",".",high.sample),
#                           "High risk", "Low risk")

library(tidyr)
library(ggplot2)
library(ggpubr)
library(Ipaper)
es_cor_dat <- data.frame(t(immu_score))
rownames(es_cor_dat)<-gsub('.','-',rownames(es_cor_dat),fixed = T)
es_cor_dat <- es_cor_dat[risk$id,]
es_cor_dat$riskScore <- risk$riskScore

es_cor_dat$group <- ifelse(rownames(es_cor_dat) %in% high.sample,
                           "High", "Low")
es_cor_dat$riskScore<-as.numeric(es_cor_dat$riskScore)
imm_cor <- ggplot(es_cor_dat, aes(x = riskScore,
                                  y = ImmuneScore)) +
  geom_point(aes(x = riskScore,
                 y = ImmuneScore,
                 color = group),
             es_cor_dat) +
  stat_cor(method = "pearson") +
  scale_color_manual(values = c("#357EBDFF", "grey"), name = "Group") +
  geom_smooth(method = "lm",
              se = T,
              color = "red",
              formula = 'y ~ x') +
  theme_bw()+
  theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
        panel.grid = element_blank()) + 
  theme(axis.text.x=element_text(family="Times",size=12,face="bold"), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Times",size=12,face="bold"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 15,face="bold"), #设置y轴标题的字体属性
        axis.title.x=element_text(family="Times",size = 15,face="bold"),
        plot.title = element_text(hjust = 0.5,family="Times",size = 15,face="bold"),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "top")
imm_cor
write_fig(imm_cor,
          file = "immune_cor1.pdf",
          width = 6,
          height = 6,
          devices = NULL,
          res = 600,
          show = F)
write_fig(imm_cor,
          file = "immune_cor1.png",
          width = 6,
          height = 6,
          devices = NULL,
          res = 600,
          show = F)
stromal_cor <- ggplot(es_cor_dat, aes(x = riskScore,
                                      y = StromalScore)) +
  geom_point(aes(x = riskScore,
                 y = StromalScore,
                 color = group),
             es_cor_dat) +
  stat_cor(method = "pearson") +
  scale_color_manual(values = c("#357EBDFF", "grey"), name = "Group") +
  geom_smooth(method = "lm",
              se = T,
              color = "red",
              formula = 'y ~ x') +
  theme_bw()+
  theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
        panel.grid = element_blank()) + 
  theme(axis.text.x=element_text(family="Times",size=12,face="bold"), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Times",size=12,face="bold"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 15,face="bold"), #设置y轴标题的字体属性
        axis.title.x=element_text(family="Times",size = 15,face="bold"),
        plot.title = element_text(hjust = 0.5,family="Times",size = 15,face="bold"),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "top")
stromal_cor
write_fig(stromal_cor,
          file = "stromal_cor1.pdf",
          width = 6,
          height = 6,
          devices = NULL,
          res = 600,
          show = F)
write_fig(stromal_cor,
          file = "stromal_cor1.png",
          width = 6,
          height = 6,
          devices = NULL,
          res = 600,
          show = F)
estimate_cor <- ggplot(es_cor_dat, aes(x = riskScore,
                                       y = ESTIMATEScore)) +
  geom_point(aes(x = riskScore,
                 y = ESTIMATEScore,
                 color = group),
             es_cor_dat) +
  stat_cor(method = "pearson") +
  scale_color_manual(values = c("#357EBDFF", "grey"), name = "Group") +
  geom_smooth(method = "lm",
              se = T,
              color = "red",
              formula = 'y ~ x') +
  theme_bw()+
  theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
        panel.grid = element_blank()) + 
  theme(axis.text.x=element_text(family="Times",size=12,face="bold"), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Times",size=12,face="bold"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 15,face="bold"), #设置y轴标题的字体属性
        axis.title.x=element_text(family="Times",size = 15,face="bold"),
        plot.title = element_text(hjust = 0.5,family="Times",size = 15,face="bold"),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "top")
estimate_cor
write_fig(estimate_cor,
          file = "estimate_cor1.pdf",
          width = 6,
          height = 6,
          devices = NULL,
          res = 600,
          show = F)
write_fig(estimate_cor,
          file = "estimate_cor1.png",
          width = 6,
          height = 6,
          devices = NULL,
          res = 600,
          show = F)
purity_cor <- ggplot(es_cor_dat, aes(x = riskScore,
                                     y = TumorPurity)) +
  geom_point(aes(x = riskScore,
                 y = TumorPurity,
                 color = group),
             es_cor_dat) +
  stat_cor(method = "pearson") +
  scale_color_manual(values = c("#357EBDFF", "grey"), name = "Group") +
  geom_smooth(method = "lm",
              se = T,
              color = "red") +
  theme_bw()+
  theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
        panel.grid = element_blank()) + 
  theme(axis.text.x=element_text(family="Times",size=12,face="bold"), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
        axis.text.y=element_text(family="Times",size=12,face="bold"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.y=element_text(family="Times",size = 15,face="bold"), #设置y轴标题的字体属性
        axis.title.x=element_text(family="Times",size = 15,face="bold"),
        plot.title = element_text(hjust = 0.5,family="Times",size = 15,face="bold"),
        legend.text = element_text(face = "bold"),
        legend.title = element_text(face = "bold"),
        legend.position = "top")
purity_cor
write_fig(purity_cor,
          file = "purity_cor1.pdf",
          width = 6,
          height = 6,
          devices = NULL,
          res = 600,
          show = F)
write_fig(purity_cor,
          file = "purity_cor1.png",
          width = 6,
          height = 6,
          devices = NULL,
          res = 600,
          show = F)

riskscore_dis <- ggplot(es_cor_dat, aes(x = group, 
                                        y = riskScore,
                                        fill = group)) +
  geom_boxplot(key_glyph = draw_key_label) +
  scale_fill_manual(values = c("#357EBDFF", "grey"), name = "Group") +
  theme_minimal() +
  labs(title = "", x = "", y = "") +
  theme(legend.title = element_text(face = "bold"),
        legend.text = element_text(face = "bold"),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        legend.position = "top",
        plot.margin = margin(t = 0,
                             l = 0,
                             b = 0, 
                             r = 0)) +
  coord_flip()
riskscore_dis

my_comparisons = list(c("High", "Low"))
imm_dis <- ggboxplot(es_cor_dat,
                     x = "group",
                     y = "ImmuneScore",
                     fill = "group",
                     palette =c( "grey","#357EBDFF")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  theme_minimal() +
  labs(title = "", x = "", y = "") +
  theme(legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        legend.position = "none",
        plot.margin = margin(t = 0,
                             l = 0,
                             b = 0, 
                             r = 0))
imm_dis
stromal_dis <- ggboxplot(es_cor_dat,
                         x = "group",
                         y = "StromalScore",
                         fill = "group",
                         palette =c("grey","#357EBDFF")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  theme_minimal() +
  labs(title = "", x = "", y = "") +
  theme(legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        legend.position = "none",
        plot.margin = margin(t = 0,
                             l = 0,
                             b = 0, 
                             r = 0))
stromal_dis
estimate_dis <- ggboxplot(es_cor_dat,
                          x = "group",
                          y = "ESTIMATEScore",
                          fill = "group",
                          palette =c("grey","#357EBDFF")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  theme_minimal() +
  labs(title = "", x = "", y = "") +
  theme(legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        legend.position = "none",
        plot.margin = margin(t = 0,
                             l = 0,
                             b = 0, 
                             r = 0))
estimate_dis
purity_dis <- ggboxplot(es_cor_dat,
                        x = "group",
                        y = "TumorPurity",
                        fill = "group",
                        palette =c("grey","#357EBDFF")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  theme_minimal() +
  labs(title = "", x = "", y = "") +
  theme(legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        legend.position = "none",
        plot.margin = margin(t = 0,
                             l = 0,
                             b = 0, 
                             r = 0))
purity_dis

library(patchwork)

{
  imm_cor2 <- imm_cor + theme(legend.position = "none",
                              plot.margin = margin(t = 0,
                                                   l = 0,
                                                   b = 0, 
                                                   r = 0))
  imm_cor2
  
  # layout <- "
  #   1##
  #   23#
  #   23#
  #   23#
  #   23#
  # "
  # imm_all <- riskscore_dis + imm_cor2  + imm_dis + 
  #   plot_layout(design = layout)
  # imm_all
  imm_all <- riskscore_dis + plot_spacer() + imm_cor2  + imm_dis +
    plot_layout(ncol = 2, heights = c(1,6), widths = c(8,1)) 
  # &
  #   theme(plot.margin = margin(t =-1, l = 1, r=-1, b =-1))
  imm_all
  write_fig(imm_all,
            file = "immune_all.pdf",
            width = 10,
            height = 8,
            devices = NULL,
            res = 600,
            show = F)  
  
  write_fig(imm_all,
            file = "immune_all.png",
            width = 10,
            height = 8,
            devices = NULL,
            res = 600,
            show = F)  
}

{
  stromal_cor2 <- stromal_cor + theme(legend.position = "none",
                                      plot.margin = margin(t = 0,
                                                           l = 0,
                                                           b = 0, 
                                                           r = 0))
  stromal_cor2
  stromal_all <- riskscore_dis + plot_spacer() + stromal_cor2 + stromal_dis +
    plot_layout(ncol =2, heights = c(1,6), widths = c(8,1))
  stromal_all
  write_fig(stromal_all,
            file = "stromal_all.pdf",
            width = 10,
            height = 8,
            devices = NULL,
            res = 600,
            show = F)  
  
  write_fig(stromal_all,
            file = "stromal_all.png",
            width = 10,
            height = 8,
            devices = NULL,
            res = 600,
            show = F)  
}

{
  estimate_cor2 <- estimate_cor + theme(legend.position = "none",
                                        plot.margin = margin(t = 0,
                                                             l = 0,
                                                             b = 0, 
                                                             r = 0))
  estimate_cor2
  estimate_all <- riskscore_dis + plot_spacer() + estimate_cor2 + estimate_dis +
    plot_layout(ncol =2, heights = c(1,6), widths = c(8,1))
  estimate_all
  write_fig(estimate_all,
            file = "estimate_all.pdf",
            width = 10,
            height = 8,
            devices = NULL,
            res = 600,
            show = F)  
  
  write_fig(estimate_all,
            file = "estimate_all.png",
            width = 10,
            height = 8,
            devices = NULL,
            res = 600,
            show = F)  
}

{
  purity_cor2 <- purity_cor + theme(legend.position = "none",
                                    plot.margin = margin(t = 0,
                                                         l = 0,
                                                         b = 0, 
                                                         r = 0))
  purity_cor2
  purity_all <- riskscore_dis + plot_spacer() + purity_cor2 + purity_dis +
    plot_layout(ncol = 2, heights = c(1,6), widths = c(8,1))
  purity_all
  write_fig(purity_all,
            file = "purity_all.pdf",
            width = 10,
            height = 8,
            devices = NULL,
            res = 600,
            show = F)  
  
  write_fig(purity_all,
            file = "purity_all.png",
            width = 10,
            height = 8,
            devices = NULL,
            res = 600,
            show = F)  
  
}

{
  es_cor_all <- imm_all / estimate_all | stromal_all / purity_all
  es_cor_all
  write_fig(es_cor_all,
            file = "es_allplot.pdf",
            width = 25,
            height = 20,
            devices = NULL,
            res = 600,
            show = F)  
  
  write_fig(es_cor_all,
            file = "es_allplot.png",
            width = 25,
            height = 20,
            devices = NULL,
            res = 600,
            show = F)  
}

