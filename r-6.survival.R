## 生存分析-----
rm(list = ls())
setwd("/data/nas1/luchunlin/project/LZZK-503")
if (! dir.exists("./06_survival")){
  dir.create("./06_survival")
}
setwd("./06_survival")
## 匹配生存数据
dat<-read.delim2('/data/nas1/luchunlin/project/LLZK-505/00_rawdata/dat.fpkm.xls',row.names = 1)%>%lc.tableToNum()
colnames<-data.frame(sample=colnames(dat))
colnames$sample<-gsub('.','-',colnames$sample,fixed = T)
colnames(dat)<-colnames$sample
survival<-read.delim2('/data/nas1/luchunlin/project/LLZK-505/03_KM/TCGA-OV.survival.tsv') 
survival<-survival[survival$sample%in%colnames(dat),]
hub<-read.delim2('/data/nas1/luchunlin/project/LZZK-503/02_correlation/hubgene.xls')
hub
dat<-t(dat)%>%as.data.frame()
## KM（MS4A7）-----
group1<-data.frame(sample=rownames(dat),group=ifelse(dat$MS4A7>median(dat$MS4A7),'High MS4A7','Low MS4A7'))
write.table(group1,file = 'group(MS4A7).xls',row.names = F,sep = '\t',quote = F)
group1$group<-as.vector(group1$group)

km.dat<-dat%>%as.data.frame()
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
                                      legend.labs=c("High MS4A7","Low MS4A7" ),
                                      legend.title="group",
                                      title="KM(MS4A7)",
                                      font.main = c(15,"bold"),
                                      risk.table = TRUE, 
                                      risk.table.col = "strata", 
                                      linetype = "strata", 
                                      surv.median.line = "hv", 
                                      ggtheme = theme_bw(), 
                                      palette = c("#A73030FF", "#0073C2FF"))
cluster_survival_median
## IFFO1------
group2<-data.frame(sample=rownames(dat),group=ifelse(dat$IFFO1>median(dat$IFFO1),'High IFFO1','Low IFFO1'))
group2$group<-as.vector(group2$group)
write.table(group2,file = 'group(IFFO1).xls',row.names = F,sep = '\t',quote = F)
km.dat<-dat%>%as.data.frame()
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
                                      legend.labs=c("High IFFO1","Low IFFO1" ),
                                      legend.title="group",
                                      title="KM(IFFO1)",
                                      font.main = c(15,"bold"),
                                      risk.table = TRUE, 
                                      risk.table.col = "strata", 
                                      linetype = "strata", 
                                      surv.median.line = "hv", 
                                      ggtheme = theme_bw(), 
                                      palette = c("#A73030FF", "#0073C2FF"))
cluster_survival_median
## ZEB2------
group3<-data.frame(sample=rownames(dat),group=ifelse(dat$ZEB2>median(dat$ZEB2),'High ZEB2','Low ZEB2'))
group3$group<-as.vector(group3$group)
write.table(group3,file = 'group(ZEB2).xls',row.names = F,sep = '\t',quote = F)
km.dat<-dat%>%as.data.frame()
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
                                      legend.labs=c("High ZEB2","Low ZEB2" ),
                                      legend.title="group",
                                      title="KM(ZEB2)",
                                      font.main = c(15,"bold"),
                                      risk.table = TRUE, 
                                      risk.table.col = "strata", 
                                      linetype = "strata", 
                                      surv.median.line = "hv", 
                                      ggtheme = theme_bw(), 
                                      palette = c("#A73030FF", "#0073C2FF"))
cluster_survival_median
## GAB3-----
group4<-data.frame(sample=rownames(dat),group=ifelse(dat$GAB3>median(dat$GAB3),'High GAB3','Low GAB3'))
group4$group<-as.vector(group4$group)
write.table(group4,file = 'group(GAB3).xls',row.names = F,sep = '\t',quote = F)
km.dat<-dat%>%as.data.frame()
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
                                      legend.labs=c("High GAB3","Low GAB3" ),
                                      legend.title="group",
                                      title="KM(GAB3)",
                                      font.main = c(15,"bold"),
                                      risk.table = TRUE, 
                                      risk.table.col = "strata", 
                                      linetype = "strata", 
                                      surv.median.line = "hv", 
                                      ggtheme = theme_bw(), 
                                      palette = c("#A73030FF", "#0073C2FF"))
cluster_survival_median
## GIMAP1-----
group5<-data.frame(sample=rownames(dat),group=ifelse(dat$GIMAP1>median(dat$GIMAP1),'High GIMAP1','Low GIMAP1'))
group5$group<-as.vector(group5$group)
write.table(group5,file = 'group(GIMAP1).xls',row.names = F,sep = '\t',quote = F)
km.dat<-dat%>%as.data.frame()
km.dat$group<-as.vector(group5$group,)
km.dat$sample<-rownames(km.dat)
km.dat<-merge(survival,km.dat,by='sample')
rownames(km.dat)<-km.dat$sample
km.dat<-km.dat[,-c(1,3)]
table(km.dat$group)
kmfit<-survfit(Surv(OS.time, OS) ~ group, data =  km.dat)
cluster_survival_median <- ggsurvplot(kmfit,
                                      pval = TRUE, 
                                      conf.int = F,
                                      legend.labs=c("High GIMAP1","Low GIMAP1" ),
                                      legend.title="group",
                                      title="KM(GIMAP1)",
                                      font.main = c(15,"bold"),
                                      risk.table = TRUE, 
                                      risk.table.col = "strata", 
                                      linetype = "strata", 
                                      surv.median.line = "hv", 
                                      ggtheme = theme_bw(), 
                                      palette = c("#A73030FF", "#0073C2FF"))
cluster_survival_median
## CD37-------
group6<-data.frame(sample=rownames(dat),group=ifelse(dat$CD37>median(dat$CD37),'High CD37','Low CD37'))
group6$group<-as.vector(group6$group)
write.table(group6,file = 'group(CD37).xls',row.names = F,sep = '\t',quote = F)
km.dat<-dat%>%as.data.frame()
km.dat$group<-as.vector(group6$group,)
km.dat$sample<-rownames(km.dat)
km.dat<-merge(survival,km.dat,by='sample')
rownames(km.dat)<-km.dat$sample
km.dat<-km.dat[,-c(1,3)]
table(km.dat$group)
kmfit<-survfit(Surv(OS.time, OS) ~ group, data =  km.dat)
cluster_survival_median <- ggsurvplot(kmfit,
                                      pval = TRUE, 
                                      conf.int = F,
                                      legend.labs=c("High CD37","Low CD37" ),
                                      legend.title="group",
                                      title="KM(CD37)",
                                      font.main = c(15,"bold"),
                                      risk.table = TRUE, 
                                      risk.table.col = "strata", 
                                      linetype = "strata", 
                                      surv.median.line = "hv", 
                                      ggtheme = theme_bw(), 
                                      palette = c("#A73030FF", "#0073C2FF"))
cluster_survival_median
