rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/XA-0214-1/")
if (! dir.exists("./05_ROC")){
  dir.create("./05_ROC")
}
setwd("./05_ROC")

library(lance)
library(tidyverse)
hubgene<-read.delim2('../04_model/hubgene.xls')
dat<-read.delim2('../00_rawdata/dat(GSE97537).xls',row.names = 1)%>%lc.tableToNum()
group<-read.delim2("../00_rawdata/group(GSE97537).xls")
table(group$group)
group$group = factor(group$group, levels = c("CIRI", "control"))
hub_exp<-dat[hubgene$symbol,]
library(ROCR)
library(ggplot2)

hub_exp2<-t(hub_exp)
## 绘制ROC曲线
library(pROC)
for (i in c(1:5)) {
  roc<-roc(group$group,hub_exp2[,i],levels=c("CIRI", "control"))
 # png(paste0(i, ".", colnames(hub_exp2)[i],".png"),width = 300,height = 300)
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
  # dev.off()
  # pdf(paste0(i, ".", colnames(hub_exp2)[i],".pdf"),width = 4,height = 4)
  # plot(roc,
  #      print.auc=T,
  #      print.auc.x=0.4,print.auc.y=0.5,
  #      #auc.polygon=T,
  #      #auc.polygon.con="#fff7f7",
  #      grid=c(0.5,0.2),
  #      grid.col=c("black","black"),
  #      #print.thres=T,
  #      main=colnames(hub_exp2)[i],
  #      col="#FF2E63",
  #      legacy.axes=T)
  # dev.off()
  i<-i+1
}

