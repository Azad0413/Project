rm(list = ls())
setwd("/data/nas1/luchunlin/project/JNZK-218-8/")
if (! dir.exists("./11_validation(GSE138125)")){
  dir.create("./11_validation(GSE138125)")
}
setwd("./11_validation(GSE138125)")
library(lance)
library(tidyverse)
hubgene<-read.delim2('/data/nas1/luchunlin/project/JNZK-218-8/06_PPI/hubgene.xls')
dat<-read.delim2('/data/nas1/luchunlin/project/JNZK-218-8/00_rawdata/dat(GSE138125).xls',row.names = 1)%>%lc.tableToNum()
group<-read.delim2("/data/nas1/luchunlin/project/JNZK-218-8/11_validation(GSE138125)/group(GSE138125).xls")
table(group$group)
group$group = factor(group$group, levels = c("POAG", "control"))
hub_exp<-dat[hubgene$symbol,]
library(ROCR)
library(ggplot2)
hub_exp2<-t(hub_exp)
## 绘制ROC曲线
library(pROC)
for (i in c(1:17)) {
  roc<-roc(group$group,hub_exp2[,i],levels=c("POAG", "control"))
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
##PDLIM1
library(pROC)
roc<-roc(group$group,hub_exp2[,10],levels=c("POAG", "control"),direction='>')
png("10.PDLIM1.png",width = 300,height = 300)
plot(roc,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main=colnames(hub_exp2)[10],
     col="#FF2E63",
     legacy.axes=T)
dev.off()
pdf("10.PDLIM1.pdf",width = 4,height = 4)
plot(roc,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main=colnames(hub_exp2)[10],
     col="#FF2E63",
     legacy.axes=T)
dev.off()

##AS  MEIS1  C3  HBG1  KRT19  JAM3  DYNC1I1  MYH10
##POAG  TNNT3  DYNC1I1  JAM3  KRT19  HBD  TNNC1  HBG1 HBB  GRP HBA2  MEIS1
##ALL  MEIS1  HBG1   KRT19  JAM3  DYNC1I1
hub.final<-data.frame(symbol=c('MEIS1','HBG1','KRT19','JAM3','DYNC1I1'))
write.table(hub.final,file = 'hub.final.xls',sep = '\t',row.names = F,quote = F)
