rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-386-10/")
if (! dir.exists("./15_roc(validation)")){
  dir.create("./15_roc(validation)")
}
setwd("./15_roc(validation)")
dat<-read.delim2('../00_rawdata/dat.va.xls',row.names = 1)%>%lc.tableToNum()
dat <- na.omit(dat)
group<-read.delim2('../00_rawdata/group.va.xls')
table(group$group)
group$group = factor(group$group, levels = c("Tumor", "Normal"))
hub_exp<-dat
library(ROCR)
library(ggplot2)
hub_exp2<-t(hub_exp)

## 绘制ROC曲线
library(pROC)
for (i in c(1:7)) {
  roc<-roc(group$group,hub_exp2[,i],levels=c("Tumor", "Normal"))
  png(paste0(i, ".", colnames(hub_exp2)[i],".png"),width = 300,height = 300)
  plot(roc,
       print.auc=T,
       print.auc.x=0.4,print.auc.y=0.5,
       #auc.polygon=T,
       #auc.polygon.con="#fff7f7",
       grid=c(0.5,0.2),
       grid.col=c("black","black"),
       #print.thres=T,
       main=colnames(hub_exp2)[i],
       col="#FF2E63",
       legacy.axes=T)
  dev.off()
  pdf(paste0(i, ".", colnames(hub_exp2)[i],".pdf"),width = 4,height = 4)
  plot(roc,
       print.auc=T,
       print.auc.x=0.4,print.auc.y=0.5,
       #auc.polygon=T,
       #auc.polygon.con="#fff7f7",
       grid=c(0.5,0.2),
       grid.col=c("black","black"),
       #print.thres=T,
       main=colnames(hub_exp2)[i],
       col="#FF2E63",
       legacy.axes=T)
  dev.off()
  i<-i+1
}

