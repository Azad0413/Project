rm(list = ls())
setwd("/data/nas1/luchunlin/project/GY0315/")
if (! dir.exists("./06_roc")){
  dir.create("./06_roc")
}
setwd("./06_roc")
library(lance)
library(tidyverse)
hubgene<-read.delim2('/data/nas1/luchunlin/project/GY0315/03_DEmodERS/intersect.xls')
dat<-read.delim2('/data/nas1/luchunlin/project/GY0315/00_rawdata/dat.xls',row.names = 1)%>%lc.tableToNum()
group<-read.delim2("/data/nas1/luchunlin/project/GY0315/00_rawdata/group.xls")
table(group$group)
group$group = factor(group$group, levels = c("High BMD", "Low BMD"))
hub_exp<-dat[hubgene$.,]
library(ROCR)
library(ggplot2)
hub_exp2<-t(hub_exp)
## 绘制ROC曲线
library(pROC)
for (i in c(1:10)) {
  roc<-roc(group$group,hub_exp2[,i],levels=c("High BMD", "Low BMD"))
  png(paste0(i, ".", colnames(hub_exp2)[i],".png"),width = 300,height = 300)
  plot(roc,
       print.auc=T,
       print.auc.x=0.4,print.auc.y=0.5,
       print.auc.pattern='%.2f',
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
       print.auc.pattern='%.2f',
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
# 1    RAC1
# 2    RPN2
# 3   FOXO3
# 4     GLA
# 5   ABCD4
# 6   LPIN1
# 7  ERGIC2
# 8  SEC22B
# 10  MYO9A


