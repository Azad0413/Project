## 04 KM生存分析-----
rm(list = ls())
setwd("/data/nas1/luchunlin/project/LLZK-505")
if (! dir.exists("./03_KM")){
  dir.create("./03_KM")
}
setwd("./03_KM")
library(survival)
library(survminer)
dat<-read.delim2('/data/nas1/luchunlin/project/LLZK-505/00_rawdata/dat.tcga.xls',row.names = 1)%>%lc.tableToNum()
colnames<-data.frame(sample=colnames(dat))
colnames$sample<-gsub('.','-',colnames$sample,fixed = T)
colnames(dat)<-colnames$sample
survival<-read_tsv(file = 'TCGA-OV.survival.tsv') 
survival<-survival[survival$sample%in%colnames(dat),]
## 04-1 GRIM-19(NDUFA13)-------
group1<-read.delim2('/data/nas1/luchunlin/project/LLZK-505/02_clinical/group(GRIM-19).xls')
group1$group<-as.vector(group1$group)
km.dat<-t(dat)%>%as.data.frame()
km.dat$group<-as.vector(group1$group,)
km.dat$sample<-rownames(km.dat)
km.dat<-merge(survival,km.dat,by='sample')
rownames(km.dat)<-km.dat$sample
km.dat<-km.dat[,-c(1,3)]
table(km.dat$group)
kmfit<-survfit(Surv(OS.time, OS) ~ group, data =  km.dat)
cluster_survival_median <- ggsurvplot(kmfit,
                                      pval = TRUE, 
                                      conf.int = F,
                                      legend.labs=c("High GRIM-19","Low GRIM-19" ),
                                      legend.title="group",
                                      title="Train KM",
                                      font.main = c(15,"bold"),
                                      risk.table = TRUE, 
                                      risk.table.col = "strata", 
                                      linetype = "strata", 
                                      surv.median.line = "hv", 
                                      ggtheme = theme_bw(), 
                                      palette = c("#A73030FF", "#0073C2FF"))
cluster_survival_median
## 04-2 NDUFS3------
group2<-read.delim2('/data/nas1/luchunlin/project/LLZK-505/02_clinical/group(NDUFS3).xls')
group2$group<-as.vector(group2$group)
km.dat<-t(dat)%>%as.data.frame()
km.dat$group<-as.vector(group2$group,)
km.dat$sample<-rownames(km.dat)
km.dat<-merge(survival,km.dat,by='sample')
rownames(km.dat)<-km.dat$sample
km.dat<-km.dat[,-c(1,3)]
table(km.dat$group)
kmfit<-survfit(Surv(OS.time, OS) ~ group, data =  km.dat)
cluster_survival_median <- ggsurvplot(kmfit,
                                      pval = TRUE, 
                                      conf.int = F,
                                      legend.labs=c("High NDUFS3","Low NDUFS3" ),
                                      legend.title="group",
                                      title="Train KM",
                                      font.main = c(15,"bold"),
                                      risk.table = TRUE, 
                                      risk.table.col = "strata", 
                                      linetype = "strata", 
                                      surv.median.line = "hv", 
                                      ggtheme = theme_bw(), 
                                      palette = c("#A73030FF", "#0073C2FF"))
cluster_survival_median
## 04-3 NDUFA4------
group3<-read.delim2('/data/nas1/luchunlin/project/LLZK-505/02_clinical/group(NDUFA4).xls')
group3$group<-as.vector(group3$group)
km.dat<-t(dat)%>%as.data.frame()
km.dat$group<-as.vector(group3$group,)
km.dat$sample<-rownames(km.dat)
km.dat<-merge(survival,km.dat,by='sample')
rownames(km.dat)<-km.dat$sample
km.dat<-km.dat[,-c(1,3)]
table(km.dat$group)
kmfit<-survfit(Surv(OS.time, OS) ~ group, data =  km.dat)
cluster_survival_median <- ggsurvplot(kmfit,
                                      pval = TRUE, 
                                      conf.int = F,
                                      legend.labs=c("High NDUFA4","Low NDUFA4" ),
                                      legend.title="group",
                                      title="Train KM",
                                      font.main = c(15,"bold"),
                                      risk.table = TRUE, 
                                      risk.table.col = "strata", 
                                      linetype = "strata", 
                                      surv.median.line = "hv", 
                                      ggtheme = theme_bw(), 
                                      palette = c("#A73030FF", "#0073C2FF"))
cluster_survival_median
## 04-4 LRPPRC------
group4<-read.delim2('/data/nas1/luchunlin/project/LLZK-505/02_clinical/group(LRPPRC).xls')
group4$group<-as.vector(group4$group)
km.dat<-t(dat)%>%as.data.frame()
km.dat$group<-as.vector(group4$group,)
km.dat$sample<-rownames(km.dat)
km.dat<-merge(survival,km.dat,by='sample')
rownames(km.dat)<-km.dat$sample
km.dat<-km.dat[,-c(1,3)]
table(km.dat$group)
kmfit<-survfit(Surv(OS.time, OS) ~ group, data =  km.dat)
cluster_survival_median <- ggsurvplot(kmfit,
                                      pval = TRUE, 
                                      conf.int = F,
                                      legend.labs=c("High LRPPRC","Low LRPPRC" ),
                                      legend.title="group",
                                      title="Train KM",
                                      font.main = c(15,"bold"),
                                      risk.table = TRUE, 
                                      risk.table.col = "strata", 
                                      linetype = "strata", 
                                      surv.median.line = "hv", 
                                      ggtheme = theme_bw(), 
                                      palette = c("#A73030FF", "#0073C2FF"))
cluster_survival_median
