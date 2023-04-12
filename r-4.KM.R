## 04 KM生存分析-----
rm(list = ls())
setwd("/data/nas1/luchunlin/project/LZZK-503")
if (! dir.exists("./04_KM")){
  dir.create("./04_KM")
}
setwd("./04_KM")
library(survival)
library(survminer)
dat<-read.delim2('/data/nas1/luchunlin/project/LLZK-505/00_rawdata/dat.tcga.xls',row.names = 1)%>%lc.tableToNum()
colnames<-data.frame(sample=colnames(dat))
colnames$sample<-gsub('.','-',colnames$sample,fixed = T)
colnames(dat)<-colnames$sample
survival<-read_tsv(file = 'TCGA-OV.survival.tsv') 
survival<-survival[survival$sample%in%colnames(dat),]
## 04-1 GRIM-19(NDUFA13)-------
group<-read.delim2('/data/nas1/luchunlin/project/LZZK-503/01_DEGs/group.xls')
group$group<-as.vector(group$group)
km.dat<-t(dat)%>%as.data.frame()
km.dat$group<-as.vector(group$group,)
km.dat$sample<-rownames(km.dat)
km.dat<-merge(survival,km.dat,by='sample')
rownames(km.dat)<-km.dat$sample
km.dat<-km.dat[,-c(1,3)]
table(km.dat$group)
kmfit<-survfit(Surv(OS.time, OS) ~ group, data =  km.dat)
cluster_survival_median <- ggsurvplot(kmfit,
                                      pval = TRUE, 
                                      conf.int = F,
                                      legend.labs=c("High ZEB2-AS1","Low ZEB2-AS1" ),
                                      legend.title="group",
                                      title="KM(ZEB2-AS1)",
                                      font.main = c(15,"bold"),
                                      risk.table = TRUE, 
                                      risk.table.col = "strata", 
                                      linetype = "strata", 
                                      surv.median.line = "hv", 
                                      ggtheme = theme_bw(), 
                                      palette = c("#A73030FF", "#0073C2FF"))
cluster_survival_median
