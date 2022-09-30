rm(list = ls())
setwd("/data/nas1/luchunlin/project/YQ444-8/")
if (! dir.exists("./04_KM(cluster)")){
  dir.create("./04_KM(cluster)")
}
setwd("./04_KM(cluster)")
library(survival)
library(survminer)
dat<-read.delim2("/data/nas1/luchunlin/project/YQ444-8/00_rawdata/dat.fpkm.xls", row.names = 1)%>% lc.tableToNum
colnames(dat)<-gsub('.','-',colnames(dat),fixed = T)                                          
dat<-dat[,-c(1:8)]
survival<-read_tsv(file = '/data/nas1/luchunlin/project/YQ444-8/02_uncox/survival.xls') 
survival$OS.time<-as.numeric(survival$OS.time)
group<-read.delim2('/data/nas1/luchunlin/project/YQ444-8/03_cluster/cluster.xls')

colnames(group)<-c('group','sample')
group$group<-as.vector(group$group)

km.dat<-t(dat)%>%as.data.frame()
km.dat<-km.dat[group$sample,]
#km.dat<-t(scale(t(km.dat)))%>%as.data.frame()
km.dat<-log2(km.dat+1)
km.dat$group<-as.vector(group$group,)
km.dat$sample<-rownames(km.dat)
km.dat<-merge(survival,km.dat,by='sample')
rownames(km.dat)<-km.dat$sample
head(km.dat[1:6])
km.dat<-km.dat[,-1]
table(km.dat$group)

kmfit<-survfit(Surv(OS.time, OS) ~ group, data =  km.dat)
cluster_survival_median <- ggsurvplot(kmfit,
                                      pval = TRUE, 
                                      conf.int = F,
                                      legend.labs=c("Cluster1","Cluster2" ),
                                      legend.title="group",
                                      title="KM curve",
                                      font.main = c(15,"bold"),
                                      risk.table = TRUE, 
                                      risk.table.col = "strata", 
                                      linetype = "strata", 
                                      surv.median.line = "hv", 
                                      ggtheme = theme_bw(), 
                                      palette = c("#A73030FF", "#0073C2FF"))
cluster_survival_median
