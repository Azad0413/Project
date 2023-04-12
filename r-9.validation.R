rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-370-8/")
if (! dir.exists("./09_validation")){
  dir.create("./09_validation")
}
setwd("./09_validation")
##GSE24673

hub_gene <- read.delim2('../06_XGBoost/hubgene.xls')
train.dat<-read.delim2("/data/nas1/luchunlin/project/BJTC-385-10/00_rawdata/dat(GSE111168).xls", row.names = 1)  %>% lc.tableToNum
group = read.delim2("/data/nas1/luchunlin/project/BJTC-385-10/00_rawdata/group(GSE111168).xls")
train.dat <- edgeR::cpm(train.dat)
train.dat <- log2(train.dat+1)
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
