rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-300-8/")
if (! dir.exists("./12_KM")){
  dir.create("./12_KM")
}
setwd("./12_KM")
dat<-read.delim2('../11_clinical/clinical_risk.xls')
library("survival")
library("survminer")
group <- read.delim2('../08_risk/risk.xls')%>%select(c('id','OS.time'))
colnames(group) <- c('sample','OS.time')
dat <- merge(group,dat,by='sample')
#dat$metastasis <- factor(dat$metastasis,levels = c('Metastatic','Non-metastatic'))
# cutpoint<-'0.5793966'
# dat$group<-ifelse(dat$group>cutpoint,'High','Low')
library(survival)
dat$OS.time <- as.numeric(dat$OS.time)
dat$OS <- as.numeric(dat$OS)
dat$group <- factor(dat$group,levels = c('High risk','Low risk'))


## Metastatic------
dat2 <- dat[which(dat$metastasis=='Metastatic'),]
# 进行KM生存分析
kmfit <- survfit(Surv(OS.time,OS) ~ group, data = dat2)

# 绘制KM曲线
KM <- ggsurvplot(kmfit,
                 pval = TRUE, 
                 conf.int = F,
                 legend.labs=c("High risk","Lowrisk" ),
                 legend.title="group",
                 title="Metastatic",
                 font.main = c(15,"bold"),
                 risk.table = TRUE, 
                 risk.table.col = "strata",  
                 linetype = "strata", 
                 surv.median.line = "hv", 
                 ggtheme = theme_bw(), 
                 palette = c("#A73030FF", "#0073C2FF"))
KM

## Non-metastatic------
dat3 <- dat[which(dat$metastasis=='Non-metastatic'),]
# 进行KM生存分析
kmfit <- survfit(Surv(OS.time,OS) ~ group, data = dat3)

# 绘制KM曲线
KM <- ggsurvplot(kmfit,
                 pval = TRUE, 
                 conf.int = F,
                 legend.labs=c("High risk","Lowrisk" ),
                 legend.title="group",
                 title="Non-Metastatic",
                 font.main = c(15,"bold"),
                 risk.table = TRUE, 
                 risk.table.col = "strata", 
                 linetype = "strata", 
                 surv.median.line = "hv", 
                 ggtheme = theme_bw(), 
                 palette = c("#A73030FF", "#0073C2FF"))
KM

## 亚型------
cluster <- read.delim2('../05_cluster/cluster.xls')
dat <- merge(cluster,dat,by='sample')
dat4 <- dat[which(dat$cluster=='cluster 1'),]
# 进行KM生存分析
kmfit <- survfit(Surv(OS.time,OS) ~ group, data = dat4)

# 绘制KM曲线
KM <- ggsurvplot(kmfit,
                 pval = TRUE, 
                 conf.int = F,
                 legend.labs=c("High risk","Lowrisk" ),
                 legend.title="group",
                 title="Cluster 1",
                 font.main = c(15,"bold"),
                 risk.table = TRUE, 
                 risk.table.col = "strata", 
                 linetype = "strata", 
                 surv.median.line = "hv", 
                 ggtheme = theme_bw(), 
                 palette = c("#A73030FF", "#0073C2FF"))
KM
dat5 <- dat[which(dat$cluster=='cluster 2'),]

# 进行KM生存分析
kmfit <- survfit(Surv(OS.time,OS) ~ group, data = dat5)

# 绘制KM曲线
KM <- ggsurvplot(kmfit,
                 pval = TRUE, 
                 conf.int = F,
                 legend.labs=c("High risk","Lowrisk" ),
                 legend.title="group",
                 title="Cluster 2",
                 font.main = c(15,"bold"),
                 risk.table = TRUE, 
                 risk.table.col = "strata", 
                 linetype = "strata", 
                 surv.median.line = "hv", 
                 ggtheme = theme_bw(), 
                 palette = c("#A73030FF", "#0073C2FF"))
KM
