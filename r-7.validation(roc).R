rm(list = ls())
setwd("/data/nas1/luchunlin/project/GY0315/")
if (! dir.exists("./07_validation(roc)")){
  dir.create("./07_validation(roc)")
}
setwd("./07_validation(roc)")
library(lance)
library(tidyverse)
hubgene<-read.delim2('/data/nas1/luchunlin/project/GY0315/03_DEmodERS/intersect.xls')
dat<-read.delim2('/data/nas1/luchunlin/project/GY0315/05_validation/dat(GSE62402).xls',row.names = 1)%>%lc.tableToNum()
group<-read.delim2("/data/nas1/luchunlin/project/GY0315/05_validation/group(GSE62402).xls")
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
  #png(paste0(i, ".", colnames(hub_exp2)[i],".png"),width = 300,height = 300)
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
  # dev.off()
  # pdf(paste0(i, ".", colnames(hub_exp2)[i],".pdf"),width = 4,height = 4)
  # plot(roc,
  #      print.auc=T,
  #      print.auc.x=0.4,print.auc.y=0.5,
  #      print.auc.pattern='%.2f',
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
## RPN2  FOXO3  ERGIC2  MYO9A  
hubgene<-data.frame(symbol=c('RPN2','FOXO3','ERGIC2','MYO9A'))
write.table(hubgene,file = 'hub.final.xls',sep = '\t',row.names = F,quote = F)
##CTNNB1 SEC22B LPIN1  GLA

roc<-roc(group$group,hub_exp2[,9],levels=c("High BMD", "Low BMD"),direction='<')
png(paste0(9, ".", colnames(hub_exp2)[9],".png"),width = 300,height = 300)
plot(roc,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     print.auc.pattern='%.2f',
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main=colnames(hub_exp2)[9],
     col="#FF2E63",
     legacy.axes=T)
dev.off()
pdf(paste0(9, ".", colnames(hub_exp2)[9],".pdf"),width = 4,height = 4)
plot(roc,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     print.auc.pattern='%.2f',
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main=colnames(hub_exp2)[9],
     col="#FF2E63",
     legacy.axes=T)
dev.off()



roc<-roc(group$group,hub_exp2[,8],levels=c("High BMD", "Low BMD"),direction='<')
png(paste0(8, ".", colnames(hub_exp2)[8],".png"),width = 300,height = 300)
plot(roc,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     print.auc.pattern='%.2f',
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main=colnames(hub_exp2)[8],
     col="#FF2E63",
     legacy.axes=T)
dev.off()
pdf(paste0(8, ".", colnames(hub_exp2)[8],".pdf"),width = 4,height = 4)
plot(roc,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     print.auc.pattern='%.2f',
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main=colnames(hub_exp2)[8],
     col="#FF2E63",
     legacy.axes=T)
dev.off()


roc<-roc(group$group,hub_exp2[,4],levels=c("High BMD", "Low BMD"),direction='>')
png(paste0(4, ".", colnames(hub_exp2)[4],".png"),width = 300,height = 300)
plot(roc,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     print.auc.pattern='%.2f',
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main=colnames(hub_exp2)[4],
     col="#FF2E63",
     legacy.axes=T)
dev.off()
pdf(paste0(4, ".", colnames(hub_exp2)[4],".pdf"),width = 4,height = 4)
plot(roc,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     print.auc.pattern='%.2f',
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main=colnames(hub_exp2)[4],
     col="#FF2E63",
     legacy.axes=T)
dev.off()


roc<-roc(group$group,hub_exp2[,6],levels=c("High BMD", "Low BMD"),direction='>')
png(paste0(6, ".", colnames(hub_exp2)[6],".png"),width = 300,height = 300)
plot(roc,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     print.auc.pattern='%.2f',
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main=colnames(hub_exp2)[6],
     col="#FF2E63",
     legacy.axes=T)
dev.off()
pdf(paste0(6, ".", colnames(hub_exp2)[6],".pdf"),width = 4,height = 4)
plot(roc,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     print.auc.pattern='%.2f',
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main=colnames(hub_exp2)[6],
     col="#FF2E63",
     legacy.axes=T)
dev.off()
