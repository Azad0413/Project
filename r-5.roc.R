## ROC----------
rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-334")
if (! dir.exists("./06_roc")){
  dir.create("./06_roc")
}
setwd("./06_roc")
## 06-1 ROC曲线验证--------
## 将基因的矩阵提取出来
hubgene<-read.delim2('/data/nas1/luchunlin/project/BJTC-334/05_expression/hubgene.xls')
dat = read.delim2("/data/nas1/luchunlin/project/BJTC-334/00_rawdata/dat.final.xls", row.names = 1) %>% lc.tableToNum
dat<-log2(dat+1)
group = read.delim2("/data/nas1/luchunlin/project/BJTC-334/00_rawdata/group.xls")
group$group = factor(group$group, levels = c("Control", "CAD"))
hub_exp<-dat[hubgene$hubgene,]
library(pROC)
library(ggplot2)
hub_exp2<-t(hub_exp)
hub_exp2<-cbind(group,hub_exp2)
## 绘制ROC曲线

roc_IL6<-roc(hub_exp2$group,hub_exp2$IL6,
                 levels=c("Control","CAD"))
roc_TP53<-roc(hub_exp2$group,hub_exp2$TP53,
               levels=c("Control","CAD"))
roc_FOS<-roc(hub_exp2$group,hub_exp2$FOS,
              levels=c("Control","CAD"))
roc_IL10<-roc(hub_exp2$group,hub_exp2$IL10,
             levels=c("Control","CAD"))
roc_JUN<-roc(hub_exp2$group,hub_exp2$JUN,
             levels=c("Control","CAD"))
roc_SRC<-roc(hub_exp2$group,hub_exp2$SRC,
             levels=c("Control","CAD"))
roc_MAPK3<-roc(hub_exp2$group,hub_exp2$MAPK3,
                levels=c("Control","CAD"))
roc_MMP9<-roc(hub_exp2$group,hub_exp2$MMP9,
              levels=c("Control","CAD"))
roc_CASP3<-roc(hub_exp2$group,hub_exp2$CASP3,
              levels=c("Control","CAD"))
roc_EGFR<-roc(hub_exp2$group,hub_exp2$EGFR,
               levels=c("Control","CAD"))
plot(roc_IL6,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="IL6 ROC curve",
     col="#FF2E63",
     legacy.axes=T)
plot(roc_TP53,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="TP53 ROC curve",
     col="#FF2E63",
     legacy.axes=T)
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
plot(roc_IL10,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="IL10 ROC curve",
     col="#FF2E63",
     legacy.axes=T)

plot(roc_JUN,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="JUN ROC curve",
     col="#FF2E63",
     legacy.axes=T)
plot(roc_SRC,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="SRC ROC curve",
     col="#FF2E63",
     legacy.axes=T)
plot(roc_MAPK3,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="MAPK3 ROC curve",
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

plot(roc_CASP3,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="CASP3 ROC curve",
     col="#FF2E63",
     legacy.axes=T)

plot(roc_EGFR,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="EGFR ROC curve",
     col="#FF2E63",
     legacy.axes=T)
