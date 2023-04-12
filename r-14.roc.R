rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-386-10/")
if (! dir.exists("./14_roc")){
  dir.create("./14_roc")
}
setwd("./14_roc")
dat<-read.delim2('../00_rawdata/dat.fpkm.xls',row.names = 1)%>%lc.tableToNum()
colnames(dat)<-gsub('.','-',colnames(dat),fixed = T)
hubgene <- read.delim2('../07_hubgene/hubgene.xls')
group<-read.delim2("../01_WGCNA/group.xls")
table(group$group)
group$group = factor(group$group, levels = c("Tumor", "Normal"))
hub_exp<-dat[hubgene$symbol,]
library(ROCR)
library(ggplot2)
hub_exp2<-t(hub_exp)

## 绘制ROC曲线
library(pROC)
for (i in c(1:8)) {
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
