rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-302")
if (! dir.exists("./08_roc")){
  dir.create("./08_roc")
}
setwd("./08_roc")
hubgene<-read.delim2('/data/nas1/luchunlin/project/BJTC-302/05_PPI/hubgene.xls')
dat<-read.delim2('/data/nas1/luchunlin/project/BJTC-302/00_rawdata/dat.xls',row.names = 1)%>%lc.tableToNum()
group<-read.delim2("/data/nas1/luchunlin/project/BJTC-302/00_rawdata/group.xls")
table(group$group)
group$group = factor(group$group, levels = c("Burn", "control"))
hub_exp<-dat[hubgene$hubgene,]
library(ROCR)
library(ggplot2)
hub_exp2<-t(hub_exp)
#hub_exp2<-cbind(group,hub_exp2)
## 绘制ROC曲线
library(pROC)
for (i in c(1:8)) {
  roc<-roc(group$group,hub_exp2[,i],levels=c("Burn", "control"))
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

