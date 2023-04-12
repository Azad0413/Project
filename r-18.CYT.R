rm(list = ls())
setwd("/data/nas1/luchunlin/project/SJZZK-428-10/")
if (! dir.exists("./18_CYT")){
  dir.create("./18_CYT")
}
setwd("./18_CYT")
dat.tcga<-read.delim2("/data/nas1/luchunlin/project/SJZZK-428-10/00_rawdata/dat.fpkm.xls", row.names = 1)%>% lc.tableToNum
colnames(dat.tcga)<-gsub('.','-',colnames(dat.tcga),fixed = T)
risk<-read.delim2('../09_risk/risk.xls')
dat.tcga<-dat.tcga[,risk$id]
dat.tcga <- log2(dat.tcga+1)
gene <- c('GZMA','PRF1')
gene.dat <- dat.tcga[gene,]
cyt.dat <- (gene.dat[1,]+gene.dat[2,])/2
rownames(cyt.dat) <- 'CYT'
cyt.dat <- t(rbind(gene.dat,cyt.dat))%>%as.data.frame()
cyt.dat$Patient <- rownames(cyt.dat)
cyt_plot_dat <- data.frame(Patient=risk$id,
                            riskScore=risk$riskScore,
                            risk_group=risk$risk)

cyt_plot_dat <- merge(cyt_plot_dat,cyt.dat,by='Patient')
cyt_plot_dat$risk_group <- factor(cyt_plot_dat$risk_group,
                                  levels = c(0,1),
                                  labels = c("High Risk", "Low Risk"))
dim(cyt_plot_dat)
library(ggplot2)
library(ggpubr)

cyt_plot_dat$riskScore<-as.numeric(cyt_plot_dat$riskScore)

## 1 GZMA------
library(ggstatsplot)
gene_cor1 <- ggscatterstats(data = cyt_plot_dat,
                           x = riskScore,
                           y = GZMA,
                           centrality.para = "mean",
                           margins = "both",
                           xfill = "#A73030FF",
                           yfill = "#0073C2FF",
                           type = "pearson",
                           ylab = "GZMA expression level",
                           marginal.type = "histogram",
                           title = "Relationship between GZMA expression level and riskScore"
)+theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(size = 12))
gene_cor1
# ggsave(filename = "01.tide_riskscore_cor.png", height = 6, width = 6,tide_cor)
# ggsave(filename = "01.tide_riskscore_cor.pdf", height = 6, width = 6,tide_cor)

gene_box1 <- ggboxplot(cyt_plot_dat,
                      x = "risk_group",
                      y = "GZMA",
                      fill = "risk_group",
                      palette =c("#A73030FF", "#0073C2FF")) +
  stat_compare_means(label.y = 6) +
  theme_bw() +
  labs(title = "", x = "", y = "GZMA expression level") +
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=12),
        axis.text.x=element_text(angle = 45,hjust = 1,colour="black",face="bold",size=12), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12),
        axis.title.y=element_text(size=12,face="bold"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
gene_box1
# ggsave(filename = "02.tide_riskscore_boxplot.png", height = 5, width = 4,tide_box)
# ggsave(filename = "02.tide_riskscore_boxplot.pdf", height = 5, width = 4,tide_box)
library(patchwork)
gene_plot1 <- gene_cor1+gene_box1+
  plot_layout(ncol = 2,heights = c(8,1), widths = c(3,1))
gene_plot1
ggsave(filename = '01.GZMA_plot.pdf',gene_plot1,w=9,h=7)
ggsave(filename = '01.GZMA_plot.png',gene_plot1,w=9,h=7)

## 2 PRF1------
library(ggstatsplot)
gene_cor2 <- ggscatterstats(data = cyt_plot_dat,
                            x = riskScore,
                            y = PRF1,
                            centrality.para = "mean",
                            margins = "both",
                            xfill = "#A73030FF",
                            yfill = "#0073C2FF",
                            type = "pearson",
                            ylab = "PRF1 expression level",
                            marginal.type = "histogram",
                            title = "Relationship between PRF1 expression level and riskScore"
)+theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(size = 12))
gene_cor2
# ggsave(filename = "01.tide_riskscore_cor.png", height = 6, width = 6,tide_cor)
# ggsave(filename = "01.tide_riskscore_cor.pdf", height = 6, width = 6,tide_cor)

gene_box2 <- ggboxplot(cyt_plot_dat,
                       x = "risk_group",
                       y = "PRF1",
                       fill = "risk_group",
                       palette =c("#A73030FF", "#0073C2FF")) +
  stat_compare_means(label.y = 4.5) +
  theme_bw() +
  labs(title = "", x = "", y = "PRF1 expression level") +
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=12),
        axis.text.x=element_text(angle = 45,hjust = 1,colour="black",face="bold",size=12), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12),
        axis.title.y=element_text(size=12,face="bold"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
gene_box2
# ggsave(filename = "02.tide_riskscore_boxplot.png", height = 5, width = 4,tide_box)
# ggsave(filename = "02.tide_riskscore_boxplot.pdf", height = 5, width = 4,tide_box)
library(patchwork)
gene_plot2 <- gene_cor2+gene_box2+
  plot_layout(ncol = 2,heights = c(8,1), widths = c(3,1))
gene_plot2
ggsave(filename = '02.PRF1_plot.pdf',gene_plot2,w=9,h=7)
ggsave(filename = '02.PRF1_plot.png',gene_plot2,w=9,h=7)

## 3 CYT-------

library(ggstatsplot)
gene_cor3 <- ggscatterstats(data = cyt_plot_dat,
                            x = riskScore,
                            y = CYT,
                            centrality.para = "mean",
                            margins = "both",
                            xfill = "#A73030FF",
                            yfill = "#0073C2FF",
                            type = "pearson",
                            ylab = "CYT score",
                            marginal.type = "histogram",
                            title = "Relationship between CYT score and riskScore"
)+theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(size = 12))
gene_cor3
# ggsave(filename = "01.tide_riskscore_cor.png", height = 6, width = 6,tide_cor)
# ggsave(filename = "01.tide_riskscore_cor.pdf", height = 6, width = 6,tide_cor)

gene_box3 <- ggboxplot(cyt_plot_dat,
                       x = "risk_group",
                       y = "CYT",
                       fill = "risk_group",
                       palette =c("#A73030FF", "#0073C2FF")) +
  stat_compare_means(label.y = 4.5) +
  theme_bw() +
  labs(title = "", x = "", y = "CYT score") +
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=12),
        axis.text.x=element_text(angle = 45,hjust = 1,colour="black",face="bold",size=12), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12),
        axis.title.y=element_text(size=12,face="bold"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
gene_box3
# ggsave(filename = "02.tide_riskscore_boxplot.png", height = 5, width = 4,tide_box)
# ggsave(filename = "02.tide_riskscore_boxplot.pdf", height = 5, width = 4,tide_box)
library(patchwork)
gene_plot3 <- gene_cor3+gene_box3+
  plot_layout(ncol = 2,heights = c(8,1), widths = c(3,1))
gene_plot3
ggsave(filename = '03.CYT_plot.pdf',gene_plot3,w=9,h=7)
ggsave(filename = '03.CYT_plot.png',gene_plot3,w=9,h=7)

##4 2个基因和差异免疫细胞相关性--------
library(Hmisc)
load('../15_CIBERSORT/cibersort.Rdata')
DE.cibersort <- read.delim2('../15_CIBERSORT/DE.cibersort.xls')
diff_score <- res.cibersort%>%column_to_rownames(var = 'cell_type')
diff_score <- diff_score[DE.cibersort$ImmuneCell,]
gene <- as.data.frame(gene)
nc<-t(rbind(diff_score,gene.dat))
m=rcorr(nc)$r[1:nrow(diff_score),(ncol(nc)-length(gene$gene)+1):ncol(nc)]
m<-t(m)
p=rcorr(nc)$P[1:nrow(diff_score),(ncol(nc)-length(gene$gene)+1):ncol(nc)]
p<-t(p)
library(dplyr)
library(dplyr)
tmp = matrix(case_when(p<0.01~"**",
                       p<0.05~"*",
                       T~""),nrow = nrow(p))
source("modified_pheatmap.R")
pdf('04.heatmap.pdf',w=5,h=4)
pheatmap(m,
         display_numbers =tmp,
         angle_col =45,
         color = colorRampPalette(c("#92b7d1", "white", "#d71e22"))(100),
         border_color = "white",
         treeheight_col = 0,
         treeheight_row = 0,
         fontsize_row = 12,
         fontsize_col = 12,
         fontsize = 10)

dev.off()

png('04.heatmap.png',w=400,h=300)
pheatmap(m,
         display_numbers =tmp,
         angle_col =45,
         color = colorRampPalette(c("#92b7d1", "white", "#d71e22"))(100),
         border_color = "white",
         treeheight_col = 0,
         treeheight_row = 0,
         fontsize_row = 12,
         fontsize_col = 12,
         fontsize = 10)

dev.off()
write.table(m,'corr.xls',
            sep = '\t',
            row.names = T)
write.table(p,'pvalue.xls',
            sep = '\t',
            row.names = T)

