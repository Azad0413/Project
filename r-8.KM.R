rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-386-10/")
if (! dir.exists("./08_KM")){
  dir.create("./08_KM")
}
setwd("./08_KM")
dat<-read.delim2('../00_rawdata/dat.fpkm.xls',row.names = 1)%>%lc.tableToNum()
colnames(dat)<-gsub('.','-',colnames(dat),fixed = T)
## 保留tumor
dat <- dat[,-c(1:72)]
hubgene <- read.delim2('../07_hubgene/hubgene.xls')
## KM curve
survival <- read.delim2('/data/nas1/luchunlin/TCGA_survival/TCGA-KIRC.survival.tsv')%>%select(-'X_PATIENT')%>%
  column_to_rownames(var = 'sample')
survival <- survival[colnames(dat),]
write.table(survival,file = 'survival.xls',sep = '\t',row.names = T,quote = F)
dat.km <- dat[hubgene$symbol,]
dat.km <- t(dat.km[,rownames(survival)])%>%as.data.frame()

dat.km <- cbind(dat.km,survival)
dat.km$OS <- as.numeric(dat.km$OS)
dat.km$OS.time <- as.numeric(dat.km$OS.time)
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

##PTCSC3 RP11-321G12.1