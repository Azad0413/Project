rm(list = ls())
setwd("/data/nas1/luchunlin/project/JNZK-214-8/")
if (! dir.exists("./cellkm")){
  dir.create("./cellkm")
}
setwd("./cellkm")
library(tidyverse)
dat<-read.csv('../00_rawdata/dat.fpkm.xls',sep = '\t')
colnames(dat)<-gsub('.','-',colnames(dat),fixed = T)
group <- read.delim2('../01_CIBERSORT/group.xls')
tumor.sample <- group$sample[which(group$group=='Tumor')]
dat <- dat[,tumor.sample]
library(immunedeconv)
set_cibersort_binary("CIBERSORT.R")
set_cibersort_mat("LM22.txt")
res.cibersort<-deconvolute(as.matrix(dat),method = 'cibersort')
write.table(res.cibersort,'cibersort.txt',sep = '\t',col.names = T,row.names = F,quote = F)

tiics <- res.cibersort%>%column_to_rownames(var = 'cell_type')
# colnames(tiics) <- gsub('.','-',colnames(tiics),fixed = T)
de.cell <- read.csv('../02_CIBERSORT(GEO)/DEcell.xls',sep = '\t')


tumor.tiics <- tiics

tumor.tiics <- tumor.tiics[de.cell$ImmuneCell,]
## 3 KM曲线------
survival <- read.csv('../07_univariate_cox/survival.xls',sep = '\t')%>%column_to_rownames(var = 'sample')

dat.km <- t(tumor.tiics[,colnames(tumor.tiics)%in%rownames(survival)])%>%as.data.frame()
survival <- survival[rownames(dat.km),]
dat.km <- cbind(dat.km,survival)
dat.km$OS <- as.numeric(dat.km$OS)
dat.km$OS.time <- as.numeric(dat.km$OS.time)

library(survival)
library(survminer)
for (i in c(1:8)) {
  library(survival)
  res.cut <- surv_cutpoint(dat.km,time = 'OS.time',event = 'OS',variables = colnames(dat.km[i]))
  cutpoint <- res.cut$cutpoint[1]
  dat.km$group <- factor(ifelse(dat.km[,i]>cutpoint$cutpoint,'High','Low'),levels = c('High','Low'))
  kmfit <- survfit(Surv(OS.time,OS) ~ group, data = dat.km)
  KM <- ggsurvplot(kmfit,
                   pval = TRUE, 
                   conf.int = F,
                   legend.labs=c("High","Low" ),
                   legend.title="group",
                   title=colnames(dat.km[i]),
                   font.main = c(15,"bold"),
                   risk.table = TRUE, 
                   risk.table.col = "strata", 
                   linetype = "strata", 
                   surv.median.line = "hv", 
                   ggtheme = theme_bw(), 
                   palette = c("#A73030FF", "#0073C2FF"))
  print(KM)
  fn1 = paste0(i,'.',colnames(dat.km)[i],'.pdf')
  fn2 = paste0(i,'.',colnames(dat.km)[i],'.png')
  ggsave(fn1,KM$plot,width = 4,h=3,units = 'in',limitsize = 300)
  ggsave(fn2,KM$plot,width = 4,h=3,units = 'in',limitsize = 300)
  i <- i+1
}




library(survival)
library(survminer)
for (i in c(1:8)) {
  library(survival)
  dat.km$group <- factor(ifelse(dat.km[,i]>median(dat.km[,i]),'High','Low'),levels = c('High','Low'))
  kmfit <- survfit(Surv(OS.time,OS) ~ group, data = dat.km)
  KM <- ggsurvplot(kmfit,
                   pval = TRUE, 
                   conf.int = F,
                   legend.labs=c("High","Low" ),
                   legend.title="group",
                   title=colnames(dat.km[i]),
                   font.main = c(15,"bold"),
                   risk.table = TRUE, 
                   risk.table.col = "strata", 
                   linetype = "strata", 
                   surv.median.line = "hv", 
                   ggtheme = theme_bw(), 
                   palette = c("#A73030FF", "#0073C2FF"))
  print(KM)
  fn1 = paste0(i,'.',colnames(dat.km)[i],'.pdf')
  fn2 = paste0(i,'.',colnames(dat.km)[i],'.png')
  ggsave(fn1,KM$plot,width = 4,h=3,units = 'in',limitsize = 300)
  ggsave(fn2,KM$plot,width = 4,h=3,units = 'in',limitsize = 300)
  i <- i+1
}
