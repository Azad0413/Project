rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-321")
if (! dir.exists("./19_model")){
  dir.create("./19_model")
}
setwd("./19_model")
library(magrittr)
library(ggplot2)
hub_gene <- read.delim2('/data/nas1/luchunlin/project/BJTC-321/08_validation/hub_final.xls')
train.dat<-read.delim2("/data/nas1/luchunlin/project/BJTC-321/00_rawdata/dat.final.xls", row.names = 1)  %>% lc.tableToNum
group = read.delim2("/data/nas1/luchunlin/project/BJTC-321/00_rawdata/group.xls")
dat<-train.dat[hub_gene$hubgene,group$sample]%>%t%>%as.data.frame()
dat<-merge(dat,group,by.x='row.names',by.y='sample')%>%tibble::column_to_rownames(var = 'Row.names')
dat$group<-factor(dat$group,levels = c('asthma','control'))
# 随机森林
library(randomForest)
set.seed(6)
res.rf <- randomForest(group ~ ., data = dat,ntrees=1000,mtry = 2,importance = T)
plot(res.rf, main = NULL)
png("01.RF.ntree.png", width = 4, height = 4, bg = "white", units = "in", res = 300, family = "Times")
plot(res.rf, main = NULL)
dev.off()
pdf("01.RF.ntree.pdf", width = 4, height = 4, family = "Times")
plot(res.rf, main = NULL)
dev.off()
saveRDS(res.rf, "rf.models.rds")
##混淆矩阵----
ntree = which.min(res.rf$err.rate[,1])
res.rf = randomForest::randomForest(group ~ ., data = dat, ntree = ntree)
df.rf = data.frame(sample = rownames(res.rf$votes), vote = res.rf$votes[,2])
df.rf$group = group$group[match(df.rf$sample, group$sample)]
data <- res.rf$confusion[2:1,2:1] %>% t %>% as.data.frame() %>% tibble::rownames_to_column(var = "true_label")
colnames(data)[2:3] <- c("pred_patient", "pred_normal")
data <- tidyr::gather(data, pred_label, value, -1)
data$true_label<-factor(data$true_label,levels = c('control','asthma'))
#data$pred_label<-factor(data$pred_label,levels = c('asthma','control'))
p1 <- ggplot(data, aes(x = pred_label, y = true_label, fill = value)) +
  geom_tile() +
  geom_text(aes(label = value), family = "Times", size = 6) +
  scale_x_discrete(breaks = c("pred_normal", "pred_patient"),
                   labels = c("control", "asthma")) +
  scale_fill_gradient(low = "white", high = "royalblue") +
  coord_fixed() +
  labs(x="Predicted Classification", y="True Classification") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 13, face = "bold", family = "Times"),
        axis.text.y = element_text(size = 13, face = "bold", family = "Times"),
        axis.title = element_text(size = 15, face = "bold", family = "Times"),
        legend.position = "none")
p1
ggsave(filename = "02.RF.Confusion.png", height = 6, width = 6, p1)
ggsave(filename = "02.RF.Confusion.pdf", height = 6, width = 6, p1)

#根据随机森林中的不同数来绘制误差率

#变量重要性评分
importance(res.rf,type = 1)
#重要性绘图
varImpPlot(res.rf,main = 'importance')
#使用测试集对构建好的随机森林进行测试
irispred<-predict(res.rf,newdata = dat)
table(irispred,dat$group)
#数据点的边距为正确归类的比例减去被归到其他类的最大比例  #观测值被判断正确的概率图
pdf('03.RF.probality.pdf',w=6,h=6)
plot(margin(res.rf,dat$group),main='Probality graph of the observed value being judged correct',
     ylab='Probality',
     xlab='Predict the number of samples')
dev.off()
png('03.RF.probality.png',w=600,h=600)
plot(margin(res.rf,dat$group),main='Probality graph of the observed value being judged correct',
     ylab='Probality',
     xlab='Predict the number of samples')
dev.off()

## PR-----
library(PRROC)
pr.en = pr.curve(df.rf$vote, weights.class0 = df.rf$group == "control", curve = T)
png("04.RF.PR.png", width = 4, height = 4, res = 300, units = "in", bg = "white")
plot(pr.en, auc.main = T, legend = F, color = 'red', asp = 1)
dev.off()
pdf("04.RF.PR.pdf", width = 4, height = 4)
plot(pr.en, auc.main = T, legend = F, color = 'red', asp = 1)
dev.off()

## roc------
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
       title = "RF ROC curve")+
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
ggsave(filename = "05.RF.ROC.png", width = 5, height = 5)
ggsave(filename = "05.RF.ROC.pdf", width = 5, height = 5)

## DCA-----
library(dcurves)
df.rf$prob = (df.rf$vote - min(df.rf$vote))/(max(df.rf$vote) - min(df.rf$vote))
dca(group ~ prob, data = df.rf,
    label = list(prob = "RF")) %>% plot(smooth = T) +
  scale_color_manual(values = c("#7FC97F","#BEAED4","#FDC086")) +
  theme(aspect.ratio = 1, legend.position = "top",
        text = element_text(size = 12, family = "Times"),
        panel.grid = element_blank())
ggsave("06.RF.DCA.png", width = 5, height = 4, dpi = 300)
ggsave("06.RF.DCA.pdf", width = 5, height = 4, dpi = 300)

# Lasso-logistic---------
set.seed(123)
library(glmnet)
res.lasso <- cv.glmnet(as.matrix(dat[-ncol(dat)]), dat$group, family = "binomial", 
                       type.measure = "auc")

plot(res.lasso)
plot(res.lasso$glmnet.fit, xvar = 'lambda')
ggsave("01.lasso.CV.png", plot(res.lasso), width = 6, height = 5, dpi = 300, units = "in", bg = "white")
ggsave("01.lasso.CV.pdf", plot(res.lasso), width = 6, height = 5, dpi = 300, units = "in", bg = "white")
ggsave("02.lasso.Coef.png", plot(res.lasso$glmnet.fit, xvar = 'lambda'), width = 6, height = 5, dpi = 300, units = "in", bg = "white")
ggsave("02.lasso.Coef.pdf", plot(res.lasso$glmnet.fit, xvar = 'lambda'), width = 6, height = 5, dpi = 300, units = "in", bg = "white")
l.coef<-coef(res.lasso$glmnet.fit,s=res.lasso$lambda.min,exact= F)
l.coef
res.lasso$lambda.min
#0.03467889
colnames(dat)
#dat<-dat[,-c(3,4,6)]
set.seed(123)
res.lasso = glmnet(as.matrix(dat[-ncol(dat)]), dat$group, family = "binomial",  lambda = res.lasso$lambda.min)
saveRDS(list(res.lasso = res.lasso), "lasso.models.rds")

df.pred = predict.glmnet(res.lasso, newx = as.matrix(dat[-ncol(dat)]), type = "link")[,1]
df.pred = data.frame(sample = rownames(dat), score = df.pred)
df.pred$group = group$group[match(df.pred$sample, group$sample)]

cs.en = confusion.glmnet(res.lasso, newx = as.matrix(dat[-ncol(dat)]), newy = dat$group) %>% as.matrix()
cs.en = cs.en[2:1,2:1]
cs.en
cs.en <- cs.en %>% t %>% as.data.frame()
# 真阳性放在第一位
cs.en$True <- factor(cs.en$True, levels = c("control","asthma"))
cs.en$Predicted<-factor(cs.en$Predicted,levels = c('asthma','control'))
p1 <- ggplot(cs.en, aes(x = Predicted, y = True, fill = Freq)) +
  geom_tile(color = "grey") +
  geom_text(aes(label = Freq), family = "Times", size = 6) +
  scale_fill_gradient(low = "white", high = "royalblue") +
  coord_fixed() +
  labs(x="Predicted Classification", y="True Classification") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 13, face = "bold", family = "Times"),
        axis.text.y = element_text(size = 13, face = "bold", family = "Times"),
        axis.title = element_text(size = 15, face = "bold", family = "Times"),
        legend.position = "none",
        panel.grid = element_blank())
p1
ggsave(filename = "03.lasso.Confusion.png", height = 6, width = 6, p1)
ggsave(filename = "03.lasso.Confusion.pdf", height = 6, width = 6, p1)

library(PRROC)
plot(pr.en, auc.main = T, legend = F, color = 'red', asp = 1)
pr.en = pr.curve(df.pred$score, weights.class0 = df.pred$group == "control", curve = T)
png("05.lasso.PR.png", width = 4, height = 4, res = 300, units = "in", bg = "white")
plot(pr.en, auc.main = T, legend = F, color = 'red', asp = 1)
dev.off()
pdf("05.lasso.PR.pdf", width = 4, height = 4)
plot(pr.en, auc.main = T, legend = F, color = 'red', asp = 1)
dev.off()

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
       title = "lasso ROC curve")+
  annotate("text",x = 0.70,y = 0.30,
           label = paste("AUC =",signif(auc(lasso_roc),2)),
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
ggsave(filename = "04.lasso.ROC.png", width = 5, height = 5)
ggsave(filename = "04.lasso.ROC.pdf", width = 5, height = 5)

# https://www.icode9.com/content-4-1399138.html
library(dcurves)
df.pred$prob = (df.pred$score - min(df.pred$score))/(max(df.pred$score) - min(df.pred$score))
dca(group ~ prob, data = df.pred,
    label = list(prob = "Lasso Logistic")) %>% plot(smooth = T) +
  scale_color_manual(values = c("#7FC97F","#BEAED4","#FDC086")) +
  theme(aspect.ratio = 1, legend.position = "top",
        text = element_text(size = 12, family = "Times"),
        panel.grid = element_blank())
ggsave("06.lasso.DCA.png", width = 5, height = 4, dpi = 300)
ggsave("06.lasso.DCA.pdf", width = 5, height = 4, dpi = 300)

## SVM-------
library(e1071)
res.svm<-svm(group~.,data = dat,kernel="radial",cost=10,gamma=1/ncol(dat))
summary(res.svm)
saveRDS(list(res.svm = res.svm), "svm.models.rds")
library(caret)
df.svm=predict(res.svm,newdata = as.matrix(dat[-ncol(dat)]),type='link')
df.svm<-data.frame(sample=rownames(dat),group=dat$group,pred=df.svm)
df.svm$group<-ifelse(df.svm$group=='asthma',1,0)
df.svm$pred<-ifelse(df.svm$pred=='asthma',1,0)
Actual<-factor(df.svm$group,levels = c(1,0),labels = c('asthma','control'))
Predict<-factor(df.svm$pred,levels = c(1,0),labels = c('asthma','control'))
xtab<-table(Predict,Actual)
xtab
xtab<-xtab%>%t%>%as.data.frame()
xtab
# 真阳性放在第一位
xtab$Actual <- factor(xtab$Actual, levels = c("control","asthma"))
xtab$Predict<-factor(xtab$Predict,levels = c('asthma','control'))
p <- ggplot(xtab, aes(x = Predict, y = Actual, fill = Freq)) +
  geom_tile(color = "grey") +
  geom_text(aes(label = Freq), family = "Times", size = 6) +
  scale_fill_gradient(low = "white", high = "royalblue") +
  coord_fixed() +
  labs(x="Predicted Classification", y="True Classification") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 13, face = "bold", family = "Times"),
        axis.text.y = element_text(size = 13, face = "bold", family = "Times"),
        axis.title = element_text(size = 15, face = "bold", family = "Times"),
        legend.position = "none",
        panel.grid = element_blank())
p
ggsave(filename = "01.SVM.Confusion.png", height = 6, width = 6, p)
ggsave(filename = "01.SVM.Confusion.pdf", height = 6, width = 6, p)

df.svm$value<-res.svm$decision.values
library(PRROC)
pr.en = pr.curve(df.svm$pred, weights.class0 = df.svm$group == "1", curve = T)
png("02.SVM.PR.png", width = 4, height = 4, res = 300, units = "in", bg = "white")
plot(pr.en, auc.main = T, legend = F, color = 'red', asp = 1)
dev.off()
pdf("02.SVM.PR.pdf", width = 4, height = 4)
plot(pr.en, auc.main = T, legend = F, color = 'red', asp = 1)
dev.off()

library(pROC)
svm_roc <- roc(df.svm$group, df.svm$pred)
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
       title = "SVM ROC curve")+
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

library(rmda)
dat$group<-ifelse(dat$group=='asthma',1,0)
colnames(dat)
complex<-decision_curve(group ~CST1+CST4+POSTN+CST2+CCL26+SERPINB2,data =dat,family = binomial(link ='logit'),
                        thresholds = seq(0,1, by = 0.01),
                        confidence.intervals= 0.95,
                        study.design = 'case-control',
                        population.prevalence= 0.3
)
pdf("04.SVM.DCA.pdf", width = 5, height = 5)
plot_decision_curve(complex,
                    curve.names=c('SVM'),
                    cost.benefit.axis =FALSE,col= c("#FDC086","#7FC97F","#BEAED4"),
                    confidence.intervals=FALSE,
                    standardize = FALSE)
dev.off()
png("04.SVM.DCA.png", width = 500, height = 500)
plot_decision_curve(complex,
                    curve.names=c('SVM'),
                    cost.benefit.axis =FALSE,col= c("#FDC086","#7FC97F","#BEAED4"),
                    confidence.intervals=FALSE,
                    standardize = FALSE)
dev.off()
library(rms)
lrm <-lrm(group ~ CST1+CST4+POSTN+CST2+CCL26+SERPINB2, data=dat, x=TRUE, y=TRUE,maxit=1000)
cal1 <- calibrate(lrm, cmethod='KM', method='boot', B=30)#建模组中绘制校准曲线
pdf("05.SVM.calibrate.pdf",width=8,height=6, family = "Times")
par(mar = c(6,5,2,2))
plot(cal1, lwd=2, lty=1, 
     cex.lab=1.5, cex.axis=1.2, cex.main=1.5, cex.sub=1.2, 
     xlim=c(0, 1), ylim= c(0, 1), 
     xlab="Nomogram-Predicted Probability of death risk", 
     ylab="Actual death (proportion)", 
     col=c("#00468BFF", "#ED0000FF", "#42B540FF"),
     legend=FALSE)
lines(cal1[, c(1:3)], type ="l", lwd=2, pch=16, col=c("#00468BFF"))
abline(0, 1, lty=3, lwd=2) 
legend(x=.6, y=.4, legend=c("Apparent", "Bias-corrected", "Ideal"), 
       lty=c(1, 1, 2), lwd = 2, col=c("#00468BFF", "black", "black"), bty="n")
#text(x = 0.2, y = 0.8, paste0("Hosmer-Lemeshow "))
#text(x = 0.34, y = 0.8, as.expression(bquote(italic('p')==.(pval))))
dev.off()
png("05.SVM.calibrate.png",width=8,height=6, units = "in", res = 600, family = "Times")
par(mar = c(6,5,2,2))
plot(cal1, lwd=2, lty=1, 
     cex.lab=1.5, cex.axis=1.2, cex.main=1.5, cex.sub=1.2, 
     xlim=c(0, 1), ylim= c(0, 1), 
     xlab="Nomogram-Predicted Probability of death risk", 
     ylab="Actual death (proportion)", 
     col=c("#00468BFF", "#ED0000FF", "#42B540FF"),
     legend=FALSE)
lines(cal1[, c(1:3)], type ="l", lwd=2, pch=16, col=c("#00468BFF"))
abline(0, 1, lty=3, lwd=2) 
legend(x=.6, y=.4, legend=c("Apparent", "Bias-corrected", "Ideal"), 
       lty=c(1, 1, 2), lwd = 2, col=c("#00468BFF", "black", "black"), bty="n")
#text(x = 0.2, y = 0.8, paste0("Hosmer-Lemeshow "))
#text(x = 0.34, y = 0.8, as.expression(bquote(italic('p')==.(pval))))
dev.off()
