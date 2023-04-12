rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-441-3/")
if (! dir.exists("./05_survival")){
  dir.create("./05_survival")
}
setwd("./05_survival")
library(survival)
library(survminer)

survival<-read.delim2('/data/nas1/luchunlin/TCGA_survival/TCGA-LAML.survival.tsv')
dat<-read.csv('../00_rawdata/dat.fpkm.xls',sep = '\t',row.names = 1)
colnames(dat) <- gsub('.','-',colnames(dat),fixed = T)

hub_exp<-dat['UGCG',]
hub.dat <- log2(t(hub_exp)+1)%>%as.data.frame()
group <- data.frame(sample=rownames(hub.dat),group=ifelse(hub.dat$UGCG>median(hub.dat$UGCG),'High UGCG','Low UGCG'))
write.table(group,file = 'group(UGCG).xls',sep = '\t',row.names = F,quote = F)

##KM-----------
km.dat <- merge(survival,group,by='sample')
kmfit<-survfit(Surv(OS.time, OS) ~ group, data =  km.dat)
table(km.dat$group)
surv_pvalue(kmfit,method = 'Gehan-Breslow')
cluster_survival_median <- ggsurvplot(kmfit,
                                      pval = TRUE, 
                                      conf.int = F,
                                      legend.labs=c("High UGCG","Low UGCG" ),
                                      legend.title="group",
                                      title="KM(UGCG)",
                                      font.main = c(15,"bold"),
                                      risk.table = TRUE, 
                                      risk.table.col = "strata", 
                                      linetype = "strata", 
                                      surv.median.line = "hv", 
                                      ggtheme = theme_bw(), 
                                      palette = c("#A73030FF", "#0073C2FF"))
cluster_survival_median


##roc-----
roc.dat <- cbind(hub.dat,km.dat)
multi_ROC <- function(time_vector, roc_table){
  library(survivalROC)
  single_ROC <- function(single_time){
    for_ROC <- survivalROC(Stime=roc_table$OS.time,
                           status=roc_table$OS,
                           marker=roc_table$UGCG,
                           predict.time=single_time,
                           method = 'KM')
    data.frame('True_positive'=for_ROC$TP, 'False_positive'=for_ROC$FP,
               'Cut_values'=for_ROC$cut.values, 'Time'=rep(single_time, length(for_ROC$TP)),
               'AUC'=rep(for_ROC$AUC, length(for_ROC$TP)))
  }
  multi_ROC_list <- lapply(time_vector, single_ROC)
  do.call(rbind, multi_ROC_list)
}
for_multi_ROC <- multi_ROC(time_vector = c(365*seq(1,5,2)), 
                           roc_table = roc.dat)
table(for_multi_ROC$AUC)
for_multi_ROC$Time <- factor(for_multi_ROC$Time)
# 画ROC曲线 
library(scales)
library(geomROC)
library(plotROC)
auc_y1 <- round(for_multi_ROC[which(for_multi_ROC$Time==365),5][1],2)
auc_y3 <- round(for_multi_ROC[which(for_multi_ROC$Time==1095),5][1],2)
auc_y5 <- round(for_multi_ROC[which(for_multi_ROC$Time==1825),5][1],2)

ROC <- ggplot(for_multi_ROC, aes(x=False_positive,
                                 y=True_positive, 
                                 label=Cut_values, 
                                 color=Time)) + 
  scale_color_manual(breaks = c("365","1095", "1825"),
                     labels = c("1 years", "3 years","5 years"),
                     values = c("#4682B4", "#FF4040", "#20B2AA")) +
  geom_roc(labels = F, stat = 'identity') + 
  style_roc() + 
  geom_abline(slope = 1, intercept = 0, color = 'gray', linetype=2) +
  theme_bw() +
  labs(title = "Train ROC") +
  theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
        panel.grid = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold")) +
  # annotate('text', x=.75, y=.25, label=paste('AUC of 1 years =', round(auc_y1,2))) + 
  # annotate('text', x=.75, y=.15, label=paste('AUC of 3 years =', round(auc_y3,2))) + 
  # annotate('text', x=.75, y=.05, label=paste('AUC of 5 years =', round(auc_y5,2))) +
  annotate("text", x=0.75, y=c(0.25, 0.15, 0.05),
           label = c(paste('AUC of 1 years =', format(auc_y1,nsmall=2)),
                     paste('AUC of 3 years =', format(auc_y3,nsmall=2)),
                     paste('AUC of 5 years =', format(auc_y5,nsmall=2))))
ROC
ggsave('02.ROC.png', ROC,width = 5, height = 4)
ggsave('02.ROC.pdf', ROC,width = 5, height = 4)

