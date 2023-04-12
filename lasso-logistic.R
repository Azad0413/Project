rm(list = ls())
library(glmnet)
library(foreign)
library(rms)
library(pROC)
library(magrittr)
library(ggplot2)
library(lance)
library(tidyverse)
# df<-read.delim2('/data/nas1/luchunlin/project/BJTC-158(返修）/diff.xls')
# gene<-df[df$id%in%colnames(input),]
input = read.table("input.txt",sep="\t",header=T,check.names=F) 
input<-input[,-1]
train.dat<-read.delim2("/data/nas1/luchunlin/project/BJTC-158(返修）/原始数据.xls", row.names = 1)  %>% lc.tableToNum
dat<-train.dat[colnames(input),]%>%t%>%as.data.frame()
group = data.frame(sample=rownames(dat),group=c(rep('control',8),rep('OSA',34)))
dat<-merge(dat,group,by.x='row.names',by.y='sample')%>%tibble::column_to_rownames(var = 'Row.names')
dat$group<-factor(dat$group,levels = c('OSA','control'))

class(dat)
# Lasso-logistic
set.seed(123)
library(glmnet)
res.lasso <- cv.glmnet(as.matrix(dat[-ncol(dat)]), dat$group, family = "binomial",nfolds = 5)

plot(res.lasso)
plot(res.lasso$glmnet.fit, xvar = 'lambda')
ggsave("01.lasso.CV.png", plot(res.lasso), width = 6, height = 6, dpi = 300, units = "in", bg = "white")
ggsave("01.lasso.CV.pdf", plot(res.lasso), width = 6, height = 6, dpi = 300, units = "in", bg = "white")
ggsave("02.lasso.Coef.png", plot(res.lasso$glmnet.fit, xvar = 'lambda'), width = 6, height = 6, dpi = 300, units = "in", bg = "white")
ggsave("02.lasso.Coef.pdf", plot(res.lasso$glmnet.fit, xvar = 'lambda'), width = 6, height = 6, dpi = 300, units = "in", bg = "white")
l.coef<-coef(res.lasso$glmnet.fit,s=res.lasso$lambda.min,exact= F)
l.coef
res.lasso$lambda.min
## 0.01101454
colnames(dat)
#dat<-dat[,-c(3,4,6)]
set.seed(15)
res.lasso = glmnet(as.matrix(dat[-ncol(dat)]), dat$group, family = "binomial",  lambda = res.lasso$lambda.min)
saveRDS(list(res.lasso = res.lasso), "models.rds")

df.pred = predict.glmnet(res.lasso, newx = as.matrix(dat[-ncol(dat)]), type = "link")[,1]
df.pred = data.frame(sample = rownames(dat), score = df.pred)
df.pred$group = group$group[match(df.pred$sample, group$sample)]

library(pROC)
lasso_roc <- roc(df.pred$group, df.pred$score)
plot(lasso_roc,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     print.thres=T,
     main="IL6 ROC curve",
     col="#FF2E63",
     legacy.axes=T)

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
       title = "ROC curve")+
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
ggsave(filename = "03.lasso.ROC.png", width = 5, height = 5)
ggsave(filename = "03.lasso.ROC.pdf", width = 5, height = 5)

## 最佳截断值-1.554
confusion<-df.pred
confusion$pred<-ifelse(confusion$score>-1.554,'control','OSA')
table(confusion$group)
table(confusion$pred)
cs.en = data.frame(True=c('control','OSA','control','OSA'),
                   Predicted=c('OSA','OSA','control','control'),
                   Freq=c('9','25','8','0'))
cs.en
cs.en$Freq<-as.numeric(cs.en$Freq)
# 真阳性放在第一位
cs.en$True <- factor(cs.en$True, levels = c("control", "OSA"))
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
ggsave(filename = "04.lasso.Confusion.png", height = 6, width = 6, p1)
ggsave(filename = "04.lasso.Confusion.pdf", height = 6, width = 6, p1)

### 验证集 ---------
library(readxl)
library(tidyverse)
test<-read_xlsx('BJTC-158_qPCR.xlsx')
test<-column_to_rownames(test,var = 'sample')
test$sample<-ifelse(test$group=='Normal',0,1)
colnames(test)
library(pROC)
def_pred<-predict(mod, newdata=test,type = 'response')
pred<-data.frame(sample=rownames(test),pred=def_pred)
write.table(pred,'predict_test.xls',sep = '\t',row.names = F,quote = F)
roc<-multiclass.roc (as.ordered(test$sample) ,as.ordered(def_pred))
roc1<-roc(as.ordered(test$sample) ,as.ordered(def_pred))
round(auc(roc1),3)
round(ci(roc1),3)
pdf(file="AUC.pdf",width = 8,height = 6)
plot(roc1,print.auc=T, auc.polygon=T, grid=c(0.1, 0.2), grid.col=c("green","red"), 
     max.auc.polygon=T, auc.polygon.col="skyblue",print.thres=T)
dev.off()
