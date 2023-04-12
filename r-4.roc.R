rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-441-3/")
if (! dir.exists("./04_roc")){
  dir.create("./04_roc")
}
setwd("./04_roc")

library(tidyverse)
dat<-read.csv('../00_rawdata/dat(GSE9476).xls',sep = '\t',row.names = 1)
group<-read.delim2('../00_rawdata/group(GSE9476).xls')
control.sample<-group$sample[which(group$group=='control')]
hub_exp<-dat['UGCG',]
# hub_exp <- log2(hub_exp+1)
library(ROCR)
library(ggplot2)

hub_exp2<-t(hub_exp)%>%as.data.frame()
## 绘制ROC曲线
library(pROC)
roc<-roc(group$group,hub_exp2$UGCG,levels=c("AML", "control"))
png('01.roc.png',width = 300,height = 300)
plot(roc,
       print.auc=T,
       print.auc.x=0.4,print.auc.y=0.5,
       print.auc.pattern='AUC=%.3f',
       #auc.polygon=T,
       #auc.polygon.con="#fff7f7",
       grid=c(0.5,0.2),
       grid.col=c("black","black"),
       #print.thres=T,
       main='UGCG ROC curve',
       col="#FF2E63",
       legacy.axes=T)
dev.off()
pdf('01.roc.pdf',width = 4,height = 4)
plot(roc,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     print.auc.pattern='AUC=%.3f',
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main='UGCG ROC curve',
     col="#FF2E63",
     legacy.axes=T)
dev.off()

