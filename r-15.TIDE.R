rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-399-11/")
if (! dir.exists("./15_TIDE")){
  dir.create("./15_TIDE")
}
setwd("./15_TIDE")
dat.tcga<-read.delim2("../00_rawdata/dat.fpkm.xls", row.names = 1)%>% lc.tableToNum()
colnames(dat.tcga)<-gsub('.','-',colnames(dat.tcga),fixed = T)
risk<-read.delim2('../08_risk/risk.xls')
dat.tcga<-dat.tcga[,risk$id]
FPKM2TPM <- function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
# tide_dat[which(tide_dat<0)] <- 0
tide_dat <- apply(dat.tcga,2,FPKM2TPM)
tide_dat <- log2(tide_dat + 1)
rownmean <- apply(tide_dat,1,mean)
tide_dat2 <- sweep(tide_dat, 1, rownmean)
dim(tide_dat2)
write.table(tide_dat2,
            file ="tide_dat.txt",
            sep = "\t",
            quote = F,
            row.names = T)

tide_result <- read.csv("tide_results.csv",header = T)
colnames(tide_result)
#View(tide_result)
tide_result2 <- subset(tide_result, select = c("Patient", "TIDE"))
tide_plot_dat <- data.frame(Patient=risk$id,
                            riskScore=risk$riskScore,
                            risk_group=risk$risk)
tide_plot_dat <- merge(tide_plot_dat, tide_result2, by = "Patient")
tide_plot_dat$risk_group <- factor(tide_plot_dat$risk_group,
                                   level = c(0, 1),
                                   labels = c("High Risk", "Low Risk"))

dim(tide_plot_dat)
library(ggplot2)
library(ggpubr)
tide_plot_dat$riskScore<-as.numeric(tide_plot_dat$riskScore)
#install.packages('ggside')
#install.packages('ggstatsplot')
##1 TIDE-----
library(ggstatsplot)
tide_cor <- ggscatterstats(data = tide_plot_dat,
                           y = riskScore,
                           x = TIDE,
                           centrality.para = "mean",
                           margins = "both",
                           xfill = "#A73030FF",
                           yfill = "#0073C2FF",
                           type = "pearson",
                           xlab = "TIDE Prediction Score",
                           marginal.type = "histogram",
                           title = "Relationship between TIDE Prediction Score and riskScore"
)+theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(size = 12))
tide_cor
# ggsave(filename = "01.tide_riskscore_cor.png", height = 6, width = 6,tide_cor)
# ggsave(filename = "01.tide_riskscore_cor.pdf", height = 6, width = 6,tide_cor)

tide_box <- ggboxplot(tide_plot_dat,
                      x = "risk_group",
                      y = "TIDE",
                      fill = "risk_group",
                      palette =c("#A73030FF", "#0073C2FF")) +
  stat_compare_means(label.y = 3.2) +
  theme_bw() +
  labs(title = "", x = "", y = "TIDE Prediction Score") +
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=12),
        axis.text.x=element_text(angle = 45,hjust = 1,colour="black",face="bold",size=12), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12),
        axis.title.y=element_text(size=12,face="bold"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
tide_box
# ggsave(filename = "02.tide_riskscore_boxplot.png", height = 5, width = 4,tide_box)
# ggsave(filename = "02.tide_riskscore_boxplot.pdf", height = 5, width = 4,tide_box)
library(patchwork)
tide_plot <- tide_cor+tide_box+
  plot_layout(ncol = 2,heights = c(8,1), widths = c(3,1))
tide_plot

ggsave(filename = '01.TIDE_plot.pdf',tide_plot,w=9,h=7)
ggsave(filename = '01.TIDE_plot.png',tide_plot,w=9,h=7)

