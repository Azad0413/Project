## ROC----------
rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-320")
if (! dir.exists("./07_roc")){
  dir.create("./07_roc")
}
setwd("./07_roc")
## 07-1 ROC曲线验证--------
## 将7个基因的矩阵提取出来
hubgene<-read.delim2('/data/nas1/luchunlin/project/BJTC-320/06_expression/hubgene.xls')
dat = read.delim2("/data/nas1/luchunlin/project/BJTC-320/00_rawdata/dat.final.xls", row.names = 1) %>% lc.tableToNum
group = read.delim2("/data/nas1/luchunlin/project/BJTC-320/00_rawdata/group.xls")
group$group = factor(group$group, levels = c("control", "HT"))
hub_exp<-dat[hubgene$x,]
library(pROC)
library(ggplot2)
hub_exp2<-t(hub_exp)
hub_exp2<-cbind(group,hub_exp2)
## 绘制ROC曲线

roc_FOS<-roc(hub_exp2$group,hub_exp2$FOS,
              levels=c("control","HT"))
roc_TNFAIP3<-roc(hub_exp2$group,hub_exp2$TNFAIP3,
               levels=c("control","HT"))
roc_PTK2B<-roc(hub_exp2$group,hub_exp2$PTK2B,
               levels=c("control","HT"))
roc_STAT1<-roc(hub_exp2$group,hub_exp2$STAT1,
                levels=c("control","HT"))
roc_MMP9<-roc(hub_exp2$group,hub_exp2$MMP9,
               levels=c("control","HT"))

plot(roc_FOS,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="FOS ROC curve",
     col="#FF2E63",
     legacy.axes=T)
plot(roc_TNFAIP3,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="TNFAIP3 ROC curve",
     col="#FF2E63",
     legacy.axes=T)
plot(roc_PTK2B,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="PTK2B ROC curve",
     col="#FF2E63",
     legacy.axes=T)
plot(roc_STAT1,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="STAT1 ROC curve",
     col="#FF2E63",
     legacy.axes=T)
plot(roc_MMP9,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="MMP9 ROC curve",
     col="#FF2E63",
     legacy.axes=T)
