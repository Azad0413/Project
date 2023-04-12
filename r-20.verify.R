rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-321")
if (! dir.exists("./20_verify")){
  dir.create("./20_verify")
}
setwd("./20_verify")

### 验证集-------
dat<-read.delim2('/data/nas1/luchunlin/project/BJTC-321/08_va/dat_va.xls',row.names = 1)%>%lc.tableToNum()
group<-read.delim2('/data/nas1/luchunlin/project/BJTC-321/08_va/group.va.xls')
table(group$group)
dat.va<-dat[,group$sample]
hub_gene<-read.delim2('/data/nas1/luchunlin/project/BJTC-321/08_va/hub_final.xls')
test.dat<-dat.va[hub_gene$hubgene,]%>%t%>%as.data.frame()
colnames(test.dat)
test.dat$group<-factor(group$group)
## RF-----
res.rf <- readRDS("../19_model/rf.models.rds")
ntree = which.min(res.rf$err.rate[,1])
res.rf = randomForest::randomForest(group ~ ., data = test.dat, ntree = ntree)
df.rf = data.frame(sample = rownames(res.rf$votes), vote = res.rf$votes[,2])
df.rf$group = group$group[match(df.rf$sample, group$sample)]
library(pROC)
rf_roc <- roc(df.rf$group, df.rf$vote)
ggroc(rf_roc,color = "red",
      linetype = 1,
      size = 1,
      alpha = 1,
      legacy.axes = T)+
  geom_abline(intercept = 0,
              slope = 1,
              color = "grey",
              size = 1,
              linetype = 1)+
  labs(x = "False Postive Rate(1 - Specificity)",
       y = "True Positive Rate(Sensivity or Recall)",
       title = "RF ROC curve (GSE4302)")+
  annotate("text",x = 0.70,y = 0.30,
           label = paste("AUC =",signif(auc(rf_roc),2)),
           size = 5,family = "Times")+
  theme_bw()+
  theme(panel.background = element_rect(fill = "transparent"),
        panel.grid = element_blank(),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(size = 0.5,colour = "black"),
        axis.title = element_text(colour = "black",size = 15),
        axis.text = element_text(colour = "black",size = 10),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        text = element_text(size = 8,color = "black",family = "Times"))
ggsave(filename = "01.RF.ROC.png", width = 5, height = 5)
ggsave(filename = "051.RF.ROC.pdf", width = 5, height = 5)

## LASSO-----
models <- readRDS("../19_model/lasso.models.rds")
df.pred <- predict(models$res.lasso, newx = as.matrix(test.dat[-ncol(test.dat)]), type = "link") %>% as.data.frame %>% 
  tibble::rownames_to_column(var = "sample")
colnames(df.pred)[2] <- "score"
df.pred$group = group$group[match(df.pred$sample, group$sample)]
library(pROC)
lasso_roc <- roc(df.pred$group, df.pred$score)
ggroc(lasso_roc,color = "red",
      linetype = 1,
      size = 1,
      alpha = 1,
      legacy.axes = T)+
  geom_abline(intercept = 0,
              slope = 1,
              color = "grey",
              size = 1,
              linetype = 1)+
  labs(x = "False Postive Rate(1 - Specificity)",
       y = "True Positive Rate(Sensivity or Recall)",
       title = "Lasso ROC curve (GSE4302)")+
  annotate("text",x = 0.70,y = 0.30,
           label = paste("AUC =", signif(auc(lasso_roc),2)),
           size = 5,family = "Times")+
  theme_bw()+
  theme(panel.background = element_rect(fill = "transparent"),
        panel.grid = element_blank(),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(size = 0.5,colour = "black"),
        axis.title = element_text(colour = "black",size = 15),
        axis.text = element_text(colour = "black",size = 10),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        text = element_text(size = 8,color = "black",family = "Times"))
ggsave(filename = "02.lasso.ROC.verify.png", width = 5, height = 5)
ggsave(filename = "02.lasso.ROC.verify.pdf", width = 5, height = 5)
## svm------
res.svm <- readRDS("../19_model/svm.models.rds")
df.svm=predict(res.svm,newdata = as.matrix(test.dat[-ncol(test.dat)]),type='link')
df.svm<-data.frame(sample=rownames(test.dat),group=test.dat$group,pred=df.svm$res.svm)
df.svm$pred<-ifelse(df.svm$pred=='asthma',1,0)
df.svm$group<-ifelse(df.svm$group=='asthma',1,0)
library(pROC)
svm_roc <- roc(df.svm$pred, df.svm$group,direction='>')
ggroc(svm_roc,color = "red",
      linetype = 1,
      size = 1,
      alpha = 1,
      legacy.axes = T)+
  geom_abline(intercept = 0,
              slope = 1,
              color = "grey",
              size = 1,
              linetype = 1)+
  labs(x = "False Postive Rate(1 - Specificity)",
       y = "True Positive Rate(Sensivity or Recall)",
       title = "SVM ROC curve (GSE4302)")+
  annotate("text",x = 0.70,y = 0.30,
           label = paste("AUC =",signif(auc(svm_roc),2)),
           size = 5,family = "Times")+
  theme_bw()+
  theme(panel.background = element_rect(fill = "transparent"),
        panel.grid = element_blank(),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(size = 0.5,colour = "black"),
        axis.title = element_text(colour = "black",size = 15),
        axis.text = element_text(colour = "black",size = 10),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        text = element_text(size = 8,color = "black",family = "Times"))
ggsave(filename = "03.SVM.ROC.png", width = 5, height = 5)
ggsave(filename = "03.SVM.ROC.pdf", width = 5, height = 5)
