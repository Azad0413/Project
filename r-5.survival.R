rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-258/")
if (! dir.exists("./05_survival")){
  dir.create("./05_survival")
}
setwd("./05_survival")
library(tidyverse)
library(lance)
dat<-read.delim2("/data/nas1/luchunlin/project/BJTC-258/00_rawdata/dat.fpkm.xls", row.names = 1)%>% lc.tableToNum
colname<-data.frame(sample=colnames(dat))
colname$sample<-gsub('.','-',colname$sample,fixed = T)
colnames(dat)<-colname$sample
## 基于高低pcascore进行生存分析-------
library("survival")
library("survminer")
pcaScore<-read.delim2('/data/nas1/luchunlin/project/BJTC-258/04_PCA/pcaScore.xls')
survival<-read.delim2('TCGA-OV.survival.tsv')
survival<-survival[,-3]
dat <- merge(pcaScore,survival,by = "sample")
##计算最佳截断值
dat$scores<-as.numeric(dat$scores)
dat$OS.time <- as.numeric(dat$OS.time)
dat$OS <- as.numeric(dat$OS)
# res.cut<-surv_cutpoint(dat,time = 'OS.time',event = 'OS',variables = 'scores')
# summary(res.cut)
# cutpoint<-'-0.2371429'
dat$group<-ifelse(dat$scores>median(dat$scores),'High','Low')
# table(dat$group)
write.table(dat,file = 'clinical.xls',sep = '\t',row.names = F,quote = F)
# 去除掉"score"列
dat <- dat[,-2]
group<-dat[,-c(2,3)]
write.table(group,file = 'group.xls',sep = '\t',row.names = F,quote = F)
dat$group <- factor(dat$group)
# 进行KM生存分析
kmfit <- survfit(Surv(OS.time,OS) ~ group, data = dat)
# 绘制KM曲线
KM <- ggsurvplot(kmfit,
                 pval = TRUE, 
                 conf.int = F,
                 legend.labs=c("High","Low" ),
                 legend.title="group",
                 title="KM",
                 font.main = c(15,"bold"),
                 risk.table = TRUE, 
                 risk.table.col = "strata", 
                 linetype = "strata", 
                 surv.median.line = "hv", 
                 ggtheme = theme_bw(), 
                 palette = c("#A73030FF", "#0073C2FF"))
KM
pcaScore<-read.delim2('/data/nas1/luchunlin/project/BJTC-258/04_PCA/pcaScore.xls')
survival<-read.delim2('TCGA-OV.survival.tsv')
survival<-survival[,-3]
dat <- merge(pcaScore,survival,by = "sample")
dat$scores<-as.numeric(dat$scores)
dat$OS.time <- as.numeric(dat$OS.time)
dat$OS <- as.numeric(dat$OS)
dat$group<-ifelse(dat$scores>median(dat$scores),'High','Low')

library(survivalROC)
library(tidyverse)
rt = subset(dat, select = c(OS, OS.time, scores))
rt$OS.time <- rt$OS.time / 365
survivalROC_helper <- function(t) {
  
  survivalROC(Stime=rt$OS.time, status=rt$OS, marker = rt$scores, 
              
              predict.time =t, method="KM")
}
survivalROC_data <- data_frame(t = c(5,7,9)) %>%
  
  mutate(survivalROC = map(t, survivalROC_helper),
         
         ## Extract scalar AUC
         
         auc = map_dbl(survivalROC, magrittr::extract2, "AUC"),
         
         ## Put cut off dependent values in a data_frame
         
         df_survivalROC = map(survivalROC, function(obj) {
           
           as_data_frame(obj[c("cut.values","TP","FP")])
           
         })) %>%
  
  dplyr::select(-survivalROC) %>%
  
  unnest() %>%
  
  arrange(t, FP, TP)
survivalROC_data1 <- survivalROC_data %>% 
  
  mutate(auc =sprintf("%.2f",auc))%>% 
  
  unite(year, t,auc,sep = " year AUC: ")


year =factor(survivalROC_data1$year)
survivalROC_data1 %>%
  
  ggplot(mapping = aes(x = FP, y = TP)) +
  scale_color_manual(
    values = c("#BC392F", "#076DB2", "#DD822D")) +
  
  geom_path(aes(color= year))+
  
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
  
  facet_wrap( ~ year) +
  
  theme_bw() +
  labs(x = "False positive fraction", y = "True positive fraction", title = "Train ROC") +
  theme(axis.text.x = element_text(vjust = 0.5),
        
        legend.key = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 16),
        legend.position = "none",
        panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
        panel.grid = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        text = element_text(face = "bold"))

ggsave(filename = "02.train_ROC.pdf", height = 4, width = 8)
ggsave(filename = "02.train_ROC.png", height = 4, width = 8)
