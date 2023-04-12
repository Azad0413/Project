rm(list = ls())
# 免疫治疗反应预测---------
setwd("/data/nas1/luchunlin/project/BJTC-308")
if (! dir.exists("./10_immunotherapy")){
  dir.create("./10_immunotherapy")
}
setwd("./10_immunotherapy")
##  TIDE-------
##肿瘤样本
dat<-read.delim2("/data/nas1/luchunlin/project/BJTC-308/00_rawdata/dat2.xls", row.names = 1)%>% lc.tableToNum
risk<-read.delim2('/data/nas1/luchunlin/project/BJTC-308/05_risk/risk.xls')

tide_dat <- dat[,risk$id]
# fpkm转TPM
# FPKM2TPM <- function(fpkm){
#   exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
# }
# 
# tide_dat[which(tide_dat<0)] <- 0
# tide_dat <- apply(tide_dat,2,FPKM2TPM)
# tide_dat <- log2(tide_dat + 1)
rownmean <- apply(tide_dat,1,mean)
tide_dat2 <- sweep(tide_dat, 1, rownmean)
dim(tide_dat2)
write.table(tide_dat2,
            file ="tide_dat.txt",
            sep = "\t",
            quote = F,
            row.names = T)


tide_result <- read.csv("TIDE.score.csv",header = T)
View(tide_result)
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
class(tide_plot_dat$riskScore)
tide_plot_dat$riskScore<-as.numeric(tide_plot_dat$riskScore)
#install.packages('ggside')
#install.packages('ggstatsplot')
{
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
  ggsave(filename = "tide_riskscore_cor.png", height = 4, width = 6,tide_cor)
  ggsave(filename = "tide_riskscore_cor.pdf", height = 4, width = 6,tide_cor)
  
  tide_box <- ggboxplot(tide_plot_dat,
                        x = "risk_group",
                        y = "TIDE",
                        fill = "risk_group",
                        palette =c("#FF6A6A", "#20B2AA")) +
    stat_compare_means(label.y = 3.2) +
    theme_bw() +
    labs(title = "", x = "", y = "TIDE Prediction Score") +
    theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18),
          axis.text.x=element_text(colour="black",face="bold",size=15), 
          axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12),
          axis.title.y=element_text(size=18,face="bold"),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  tide_box
  ggsave(filename = "tide_riskscore_boxplot.png", height = 5, width = 4,tide_box)
  ggsave(filename = "tide_riskscore_boxplot.pdf", height = 5, width = 4,tide_box)
}
## IPS--------
## IOBR包

#if (!requireNamespace("IOBR", quietly = TRUE))
#  devtools::install_github("IOBR/IOBR")
library(IOBR)
ips_dat<-tide_dat
ips<-deconvo_tme(eset = ips_dat, method = "ips", plot= FALSE)
# ips<-IPS_calculation(eset = ips_dat,plot = F)
head(ips)
ips_result<-subset(ips,select=c('ID','IPS_IPS'))
colnames(ips_result)<-c('ID','IPS')
ips_plot_dat<-data.frame(ID=risk$id,
                         riskScore=risk$riskScore,
                         risk_group=risk$risk)
ips_plot_dat<-merge(ips_plot_dat,ips_result,by='ID')
ips_plot_dat$risk_group<-factor(ips_plot_dat$risk_group,
                                levels = c(0,1),
                                labels = c("High Risk", "Low Risk"))
# BiocManager::install('ggridges')
library('ggridges')
my_comparisons<-list(c('High Risk','Low Risk'))
ggplot(ips_plot_dat, aes(x = IPS, y = risk_group)) +
  geom_density_ridges_gradient(aes(fill = risk_group),
                               scale = 1, size = 0.3,
                               position = 'identity') +
  geom_signif(comparisons = my_comparisons,
              map_signif_level = F,
              y_position = c(6))+
  theme(legend.position = "none")+
  labs(y = "")+
  theme(axis.text.y=element_text(hjust=0.5,colour="black",size=12),
        axis.text.x=element_text(hjust=0.5,colour="black",size=12))+
  scale_fill_manual(values= c("#FF6A6A", "#20B2AA"))
