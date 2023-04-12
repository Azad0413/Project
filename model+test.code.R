rm(list = ls())
library(glmnet)
library(foreign)
library(rms)
library(pROC)
## lasso-logsitic模型------
input = read.table("input.txt",sep="\t",header=T,check.names=F) 
input$sample <- factor(input$sample)
y<-data.matrix(input[,1])
x<-as.matrix(input[,c(2:18)])
f1 =glmnet(x, y, family="binomial", nlambda=100, alpha=1)
pdf(file="01.lambda.pdf",width = 8,height = 6)
plot(f1,xvar="lambda", label=TRUE)
dev.off()
png(file="01.lambda.png",width = 800,height = 600)
plot(f1,xvar="lambda", label=TRUE)
dev.off()
predict(f1,newx=x[2:5,], type = "response")
y <- as.numeric(unlist(y))
cvfit=cv.glmnet(x,y)
pdf(file="02.cvfit.pdf",width = 8,height = 6)
plot(cvfit)
dev.off()
png(file="02.cvfit.png",width = 800,height = 600)
plot(cvfit)
dev.off()
mod<-glm(sample~AREG+ZFP36+ATF3+DUSP1,family="binomial",data= input)

predicidata <- data.frame(input[,c('AREG','ATF3','DUSP1','ZFP36')])

###AUC###
library(pROC)
df.pred<-predict(mod, newdata=predicidata,type = 'link')
df.pred = data.frame(sample = rownames(input), score = df.pred)
df.pred$group = factor(ifelse(input$sample=='0','control','OSA'),levels = c('control','OSA'))
roc <- roc(df.pred$group, df.pred$score)
plot(roc,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     print.thres=T,
     main="ROC curve",
     col="#FF2E63",
     legacy.axes=T)
### 1.903
ggroc(roc,color = "red",
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
       title = "ROC curve")+
  annotate("text",x = 0.70,y = 0.30,
           label = paste("AUC =",signif(auc(roc),2)),
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
ggsave(filename = "03.ROC.png", width = 5, height = 5)
ggsave(filename = "03.ROC.pdf", width = 5, height = 5)
##PR-------
library(PRROC)
pr.en = pr.curve(df.pred$score, weights.class0 = df.pred$group == "OSA", curve = T)
png("04.PR.png", width = 4, height = 4, res = 300, units = "in", bg = "white")
plot(pr.en, auc.main = T, legend = F, color = 'red', asp = 1)
dev.off()
pdf("04.PR.pdf", width = 4, height = 4)
plot(pr.en, auc.main = T, legend = F, color = 'red', asp = 1)
dev.off()

## 混淆矩阵-------
## 值-1.554
confusion<-df.pred
confusion$pred<-ifelse(confusion$score>1.903,'OSA','control')
table(confusion$group)
table(confusion$pred)
cs.en = data.frame(True=c('control','OSA','control','OSA'),
                   Predicted=c('OSA','OSA','control','control'),
                   Freq=c('0','25','8','9'))
cs.en
cs.en$Freq<-as.numeric(cs.en$Freq)
# 真阳性放在第一位
cs.en$True <- factor(cs.en$True, levels = c("control", "OSA"))
cs.en$Predicted<-factor(cs.en$Predicted,levels = c('OSA','control'))
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
ggsave(filename = "04.Confusion.png", height = 6, width = 6, p1)
ggsave(filename = "04.Confusion.pdf", height = 6, width = 6, p1)

### 验证集 ---------
library(readxl)
library(tidyverse)
test<-read_xlsx('BJTC-158_qPCR.xlsx')
test<-column_to_rownames(test,var = 'sample')
colnames(test)
library(pROC)
test.pred<-predict(mod, newdata=test[-ncol(test)],type = 'link')
test.pred = data.frame(sample = rownames(test), score = test.pred)
test.pred$group = factor(test$group,levels = c('Normal','OSA'))
write.table(test.pred,file = 'test.pred.xls',sep = '\t',row.names = F,quote = F)
roc <- roc(test.pred$group, test.pred$score)
plot(roc,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     print.thres=T,
     main="ROC curve",
     col="#FF2E63",
     legacy.axes=T)
ggroc(roc,color = "red",
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
       title = "ROC curve")+
  annotate("text",x = 0.70,y = 0.30,
           label = paste("AUC =",signif(auc(roc),2)),
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
ggsave(filename = "05.Test_ROC.png", width = 5, height = 5)
ggsave(filename = "05.Test_ROC.pdf", width = 5, height = 5)

## 16.010
test.confusion<-test.pred
test.confusion$pred<-ifelse(test.confusion$score>16.010,'control','OSA')
cs.en = data.frame(True=c('control','OSA','control','OSA'),
                   Predicted=c('OSA','OSA','control','control'),
                   Freq=c('3','13','8','7'))
cs.en
cs.en$Freq<-as.numeric(cs.en$Freq)
# 真阳性放在第一位
cs.en$True <- factor(cs.en$True, levels = c("control", "OSA"))
cs.en$Predicted<-factor(cs.en$Predicted,levels = c('OSA','control'))
p <- ggplot(cs.en, aes(x = Predicted, y = True, fill = Freq)) +
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
ggsave(filename = "06.test.confusion.png", height = 6, width = 6, p)
ggsave(filename = "06.test.confusion.pdf", height = 6, width = 6, p)

write.table(test.confusion,file = 'test_confusion.xls',sep = '\t',quote = F,row.names = F)
