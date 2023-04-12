rm(list = ls())
setwd("/data/nas1/luchunlin/project/CD-0601-2/")
if (! dir.exists("./08_ROC")){
  dir.create("./08_ROC")
}
setwd("./08_ROC")

##单基因ROC 训练集-------
library(lance)
library(tidyverse)
hubgene <- read.delim2('../07_features/features.xls')
dat<-read.delim2('../00_rawdata/dat(GSE65682).xls',row.names = 1)%>%lc.tableToNum()
group<-read.delim2("../00_rawdata/group(GSE65682).xls")
table(group$group)
group$group = factor(group$group, levels = c("control", "Sepsis"))
hub_exp<-dat[hubgene$symbol,]
# hub_exp <- log2(hub_exp+1)
library(ROCR)
library(ggplot2)

hub_exp2<-t(hub_exp)
## 绘制ROC曲线
library(pROC)
for (i in c(1:7)) {
  roc<-roc(group$group,hub_exp2[,i],levels=c("Sepsis", "control"))
  png(paste0(i, ".", colnames(hub_exp2)[i],".png"),width = 300,height = 300)
  plot(roc,
       print.auc=T,
       print.auc.x=0.4,print.auc.y=0.5,
       print.auc.pattern='AUC=%.3f',
       #auc.polygon=T,
       #auc.polygon.con="#fff7f7",
       grid=c(0.5,0.2),
       grid.col=c("black","black"),
       #print.thres=T,
       main=paste0(colnames(hub_exp2)[i]),
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
       main=paste0(colnames(hub_exp2)[i]),
       col="#FF2E63",
       legacy.axes=T)
  dev.off()
  i<-i+1
}
