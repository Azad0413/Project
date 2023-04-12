rm(list = ls())
setwd("/data/nas1/luchunlin/project/NN-0118-2/")
if (! dir.exists("./10_KM")){
  dir.create("./10_KM")
}
setwd("./10_KM")

library(survival)
library(survminer)
dat<-read.delim2('../00_rawdata/dat.tpm.xls',row.names = 1)%>%lc.tableToNum()
colnames(dat)<-gsub('.','-',colnames(dat),fixed = T)
group <- read.delim2('../01_DEGs(TCGA)/group.xls')
tumor.sample <- group[which(group$group=='Tumor'),]
dat <- dat[,tumor.sample$sample]
survival<-read_tsv(file = '../00_rawdata/survival.xls') 
survival<-survival[survival$sample%in%colnames(dat),]

## TOP2A-------
group <- dat['CCNB1',]%>%t%>%as.data.frame()
group$group <- ifelse(group$CCNB1>median(group$CCNB1),'High CCNB1','Low CCNB1')
group <- rownames_to_column(group,var = 'sample')
write.table(group,file = 'group(CCNB1).xls',sep = '\t',row.names = F,quote = F)

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
                                      legend.labs=c("High CCNB1","Low CCNB1" ),
                                      legend.title="group",
                                      title="CCNB1 KM",
                                      font.main = c(15,"bold"),
                                      risk.table = TRUE, 
                                      risk.table.col = "strata", 
                                      linetype = "strata", 
                                      surv.median.line = "hv", 
                                      ggtheme = theme_bw(), 
                                      palette = c("#A73030FF", "#0073C2FF"))
cluster_survival_median
