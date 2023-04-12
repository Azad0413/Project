rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-370-8/")
if (! dir.exists("./07_logistic")){
  dir.create("./07_logistic")
}
setwd("./07_logistic")

hub_gene <- read.delim2('../06_XGBoost/hubgene.xls')
train.dat<-read.delim2("../03_DEGs/dat_final.xls", row.names = 1)  %>% lc.tableToNum
group = read.delim2("../03_DEGs/group.xls")
dat<-train.dat[hub_gene$symbol,group$sample]%>%t%>%as.data.frame()
dat<-merge(dat,group,by.x='row.names',by.y='sample')%>%tibble::column_to_rownames(var = 'Row.names')
dat$group<-factor(dat$group,levels = c('RB','control'))
hub_gene
mod<-glm(group~CP+RRM2+NT5DC2,family="binomial",data= dat)
predicidata <- data.frame(dat[,c('CP','RRM2','NT5DC2')])
###AUC###
library(pROC)
df.pred<-predict(mod, newdata=predicidata,type = 'link')
df.pred = data.frame(sample = rownames(dat), score = df.pred)
df.pred$group = factor(dat$group,levels = c('control','RB'))
roc <- roc(df.pred$group, df.pred$score)
pdf(file = '01.ROC.pdf',w=4,h=4)
plot(roc,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="ROC curve",
     col="#FF2E63",
     legacy.axes=T)
dev.off()
png(file = '01.ROC.png',w=300,h=300)
plot(roc,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="ROC curve",
     col="#FF2E63",
     legacy.axes=T)
dev.off()
# ### 基因的ROC曲线
# CP_roc <- roc(dat$group, dat$CP)
# pdf(file = '02.CP_ROC.pdf',w=4,h=4)
# plot(roc,
#      print.auc=T,
#      print.auc.x=0.4,print.auc.y=0.5,
#      #auc.polygon=T,
#      #auc.polygon.con="#fff7f7",
#      grid=c(0.5,0.2),
#      grid.col=c("black","black"),
#      #print.thres=T,
#      main="CP ROC curve",
#      col="#FF2E63",
#      legacy.axes=T)
# 
# dev.off()
# png(file = '02.CP_ROC.png',w=300,h=300)
# plot(roc,
#      print.auc=T,
#      print.auc.x=0.4,print.auc.y=0.5,
#      #auc.polygon=T,
#      #auc.polygon.con="#fff7f7",
#      grid=c(0.5,0.2),
#      grid.col=c("black","black"),
#      #print.thres=T,
#      main="CP ROC curve",
#      col="#FF2E63",
#      legacy.axes=T)
# dev.off()
# RRM2_roc <- roc(dat$group, dat$RRM2)
# pdf(file = '03.RRM2_ROC.pdf',w=4,h=4)
# plot(roc,
#      print.auc=T,
#      print.auc.x=0.4,print.auc.y=0.5,
#      # add = T,
#      #auc.polygon=T,
#      #auc.polygon.con="#fff7f7",
#      grid=c(0.5,0.2),
#      grid.col=c("black","black"),
#      #print.thres=T,
#      main="RRM2 ROC curve",
#      col="#FF2E63",
#      legacy.axes=T)
# dev.off()
# png(file = '03.RRM2_ROC.png',w=300,h=300)
# plot(roc,
#      print.auc=T,
#      print.auc.x=0.4,print.auc.y=0.5,
#      # add = T,
#      #auc.polygon=T,
#      #auc.polygon.con="#fff7f7",
#      grid=c(0.5,0.2),
#      grid.col=c("black","black"),
#      #print.thres=T,
#      main="RRM2 ROC curve",
#      col="#FF2E63",
#      legacy.axes=T)
# dev.off()
# 
# NT5DC2_roc <- roc(dat$group, dat$NT5DC2)
# pdf(file = '04.NT5DC2_ROC.pdf',w=4,h=4)
# plot(roc,
#      print.auc=T,
#      print.auc.x=0.4,print.auc.y=0.5,
#      # add = T,
#      #auc.polygon=T,
#      #auc.polygon.con="#fff7f7",
#      grid=c(0.5,0.2),
#      grid.col=c("black","black"),
#      #print.thres=T,
#      main="NT5DC2 ROC curve",
#      col="#FF2E63",
#      legacy.axes=T)
# dev.off()
# png(file = '04.NT5DC2_ROC.png',w=300,h=300)
# plot(roc,
#      print.auc=T,
#      print.auc.x=0.4,print.auc.y=0.5,
#      # add = T,
#      #auc.polygon=T,
#      #auc.polygon.con="#fff7f7",
#      grid=c(0.5,0.2),
#      grid.col=c("black","black"),
#      #print.thres=T,
#      main="NT5DC2 ROC curve",
#      col="#FF2E63",
#      legacy.axes=T)
# dev.off()
