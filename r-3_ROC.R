## ROC-----
rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-327/")
if (! dir.exists("./03_ROC")){
  dir.create("./03_ROC")
}
setwd("./03_ROC")
library(lance)
library(tidyverse)
# AIH.PBC------
##将基因的矩阵提取出来
hubgene<-read.delim2('/data/nas1/luchunlin/project/BJTC-327/02_DEIRG/DEIRG.xls')
dat1<-read.delim2('/data/nas1/luchunlin/project/BJTC-327/00_rawdata/AIH.PBC.fpkm.xls',row.names = 1)%>%lc.tableToNum()
dat1<-na.omit(dat1)
group1<-read.delim2("/data/nas1/luchunlin/project/BJTC-327/01_DEGs/group1.xls")
table(group1$group)
group1$group = factor(group1$group, levels = c("AIH", "PBC"))
hub_exp1<-dat1[hubgene$.,]
library(ROCR)
library(ggplot2)
hub_exp2<-t(hub_exp1)
#hub_exp2<-cbind(group1,hub_exp2)
## 绘制ROC曲线
library(pROC)
for (i in c(1:18)) {
  roc<-roc(group1$group,hub_exp2[,i],levels=c("AIH", "PBC"))
  png(paste0(i, ".", colnames(hub_exp2)[i],"(AIH vs.PBC)",".png"),width = 300,height = 300)
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
    pdf(paste0(i, ".", colnames(hub_exp2)[i], "(AIH vs.PBC)",".pdf"),width = 4,height = 4)
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
##  AIH.HBV-------
dat2<-read.delim2('/data/nas1/luchunlin/project/BJTC-327/00_rawdata/AIH.HBV.fpkm.xls',row.names = 1)%>%lc.tableToNum()
dat2<-na.omit(dat2)
group2<-read.delim2("/data/nas1/luchunlin/project/BJTC-327/01_DEGs/group2.xls")
table(group2$group)
group2$group = factor(group2$group, levels = c("AIH", "HBV"))
hub_exp3<-dat2[hubgene$.,]
library(ROCR)
library(ggplot2)
hub_exp4<-t(hub_exp3)
#hub_exp2<-cbind(group2,hub_exp2)
## 绘制ROC曲线
library(pROC)
for (i in c(1:18)) {
  roc<-roc(group2$group,hub_exp4[,i],levels=c("AIH", "HBV"))
  png(paste0(i, ".", colnames(hub_exp4)[i],"(AIH vs.HBV)",".png"),width = 300,height = 300)
  plot(roc,
       print.auc=T,
       print.auc.x=0.4,print.auc.y=0.5,
       #auc.polygon=T,
       #auc.polygon.con="#fff7f7",
       grid=c(0.5,0.2),
       grid.col=c("black","black"),
       #print.thres=T,
       main=colnames(hub_exp4)[i],
       col="#FF2E63",
       legacy.axes=T)
  dev.off()
  pdf(paste0(i, ".", colnames(hub_exp4)[i], "(AIH vs.HBV)",".pdf"),width = 4,height = 4)
  plot(roc,
       print.auc=T,
       print.auc.x=0.4,print.auc.y=0.5,
       #auc.polygon=T,
       #auc.polygon.con="#fff7f7",
       grid=c(0.5,0.2),
       grid.col=c("black","black"),
       #print.thres=T,
       main=colnames(hub_exp4)[i],
       col="#FF2E63",
       legacy.axes=T)
  dev.off()
  i<-i+1
}
