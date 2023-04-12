rm(list = ls())
setwd("/data/nas1/luchunlin/project/TY0307-11/")
if (! dir.exists("./04_model")){
  dir.create("./04_model")
}
setwd("./04_model")
library(magrittr)
library(ggplot2)
hub_gene <- read.delim2('../02_DEERS/DEERS.xls')
train.dat<-read.delim2("../00_rawdata/dat(GSE113079).xls", row.names = 1)  %>% lc.tableToNum
group = read.delim2("../00_rawdata/group(GSE113079).xls")
dat<-train.dat[hub_gene$symbol,group$sample]%>%t%>%as.data.frame()
dat<-merge(dat,group,by.x='row.names',by.y='sample')%>%tibble::column_to_rownames(var = 'Row.names')
dat$group<-factor(dat$group,levels = c('CAD','control'))
#LASSO----------
set.seed(1)
##4
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
coef.min = coef(res.lasso, s = "lambda.min")  ## lambda.min & lambda.1se 取一个
res.lasso$lambda.min
# 找出那些回归系数没有被惩罚为0的
active.min = which(coef.min@i != 0)
# 提取基因名称
lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i+1]
lasso_geneids <- lasso_geneids[-1]
lasso_geneids
# [1] "RCN2"     "HRC"      "DERL2"    "RNF183"   "CRH"      "TMED2"    "PPP1R15A" "IL1A"  
write.table(lasso_geneids,file = 'hubgene.xls',sep = '\t',row.names = F,quote = F)
colnames(dat)
#dat<-dat[,-c(3,4,6)]
set.seed(123)

res.lasso$lambda.min
dat <- dat[,c(lasso_geneids,'group')]
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
cs.en$True <- factor(cs.en$True, levels = c("control","CAD"))
cs.en$Predicted<-factor(cs.en$Predicted,levels = c('CAD','control'))
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

## 人工神经网络--------
set.seed(1002332)  
library(caret)
library(PRROC)
library(neuralnet)
input_data <- dat
folds <- createFolds(y=input_data$group,k=5)   ###分成5份
colnames(input_data)
class(input_data$group)

### LR模型
for  (i in 1:5){
  fold_test <- input_data[folds[[i]],]    #取fold i数据，建立测试集
  fold_train <- input_data[-folds[[i]],]  #建立训练集
  set.seed(10902893) ##108239
  fold_model <- neuralnet(
    group~RCN2+HRC+DERL2+RNF183+CRH+TMED2+PPP1R15A+IL1A,
    data=fold_train,
    threshold=0.3,  ##先将阈值从默认值 = 0.01 连续增加到 0.1、0.2、0.3 等
    hidden=6,  ##隐藏层数在输入层是输出层数   8
    err.fct = 'ce',
    linear.output = FALSE) ###人工神经网络
  net_predict<-compute(fold_model,fold_test)$net.result
  net_prediction<-c("control","CAD")[apply(net_predict,1,which.max)]
  roc_curve = roc.curve(net_prediction, weights.class0 =  fold_test$group == "CAD", curve = T)
  png(paste('0',i+4,'.roc_fold',i,'.png',sep=''), width = 7, height = 7, res = 600, units = "in", bg = "white",family='Times')
  par(pin = c(4,4), mar = c(6,6,6,1)) 
  plot(roc_curve, auc.main = T, legend = F, color = 'red', xlab = "1-Specificity", asp = 1 ,cex.axis=1.8,  ##坐标轴刻度文字的缩放倍数。类似cex。
       cex.lab=2.0,   ##坐标轴刻度文字的缩放倍数。类似cex。
       cex.main=2.0,   ##标题的缩放倍数
       main=paste("Fold-",i,'',sep=''),
       font.lab = 2, 
       font.main = 2, 
       font.sub =2)
  abline(0,1); dev.off()
  pdf(paste('0',i+4,'.roc_fold',i,'.pdf',sep=''), width = 7, height = 7, family='Times')
  par(pin = c(4,4), mar = c(6,6,6,1)) 
  plot(roc_curve, 
       auc.main = T, 
       legend = F, color ='red', xlab = "1-Specificity", asp = 1 ,cex.axis=1.8,  ##坐标轴刻度文字的缩放倍数。类似cex。
       cex.lab=2.0,   ##坐标轴刻度文字的缩放倍数。类似cex。
       cex.main=2.0,   ##标题的缩放倍数。
       main=paste("Fold-",i,sep=''),
       font.lab = 2, 
       font.main = 2,
       font.sub =2)
  abline(0,1); dev.off()
}



##可视化
i<-1
fold_test <- input_data[folds[[i]],]    #取fold i数据，建立测试集
fold_train <- input_data[-folds[[i]],]  #建立训练集
set.seed(10902893) ##108239
fold_model <- neuralnet(
  group~RCN2+HRC+DERL2+RNF183+CRH+TMED2+PPP1R15A+IL1A,
  data=fold_train,
  threshold=0.3,  ##先将阈值从默认值 = 0.01 连续增加到 0.1、0.2、0.3 等
  hidden=6,  ##隐藏层数在输入层是输出层数   8
  err.fct = 'ce',
  linear.output = FALSE) ###人工神经网络
plot(fold_model)

save(fold_model,file='model.RData')

#<span style="color:#009900;">#import the functin from Github</span>
library(devtools)
source_url('https://gist.githubusercontent.com/fawda123/7471137/raw/466c1474d0a505ff044412703516c34f1a4684a5/nnet_plot_update.r')
#<span style="color:#009900;">#plot each model</span>

pdf('10.nnet_plot.pdf', width = 9, height = 7)
par(pin = c(4,4), mar = c(1,1,1,1),family='Times',font=2) 
plot.nnet(fold_model)
dev.off()

png('10.nnet_plot.png', width = 9, height = 7,units='in',res=600)
par(pin = c(4,4), mar = c(1,1,1,1),family='Times',font=2) 
plot.nnet(fold_model)
dev.off()

fold1_results<-data.frame(row.names=rownames(fold_model$result.matrix), weights =round((fold_model$result.matrix)[,1],3))
#write.csv(fold1_results,'10.fold1_results.csv',quote=F)




# # 随机森林-LOOCV-------
# library(randomForest)
# set.seed(666)
# res.rf <- randomForest(group ~ ., data = dat,ntrees=1000,mtry = 2,importance = T)
# plot(res.rf, main = NULL)
# saveRDS(res.rf, "rf.models.rds")
# 
# #变量重要性评分
# importance(res.rf,type = 1)
# imp <- data.frame(symbol=res.rf$importance)%>%rownames_to_column(var = 'symbol')
# #重要性绘图
# varImpPlot(res.rf,main = 'importance')
# #使用测试集对构建好的随机森林进行测试
# irispred<-predict(res.rf,newdata = dat)
# table(irispred,dat$group)
# imp <- imp[order(imp$symbol.MeanDecreaseGini,decreasing = T),]
# library(boot)
# dat.cv <- dat[,imp$symbol]
# dat.cv<-merge(dat.cv,group,by.x='row.names',by.y='sample')%>%tibble::column_to_rownames(var = 'Row.names')
# table(dat.cv$group)
# dat.cv$group <- factor(ifelse(dat.cv$group=='CAD',1,0))
# auc <- data.frame(auc=rep(NA,24))
# cv.error <- data.frame(error=rep(NA,24))
# cnt <- 1
# while (cnt < 25) {
#   mod <- glm(group~.,family = 'binomial',data = dat.cv[c(1:cnt,25)])
#   library(pROC)
#   df.pred<-predict(mod, newdata=dat.cv,type = 'link')
#   df.pred = data.frame(sample=rownames(dat.cv),score=df.pred)
#   df.pred$group = factor(ifelse(dat.cv$group=='0','control','CAD'),levels = c('control','CAD'))
#   roc <- roc(df.pred$group, df.pred$score)
#   auc[cnt,] <- roc$auc
#   cv.error[cnt,] <- cv.glm(dat.cv,mod)$delta[1]
#   cnt = cnt + 1
# }
# 
