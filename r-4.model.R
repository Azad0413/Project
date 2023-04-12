rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/XA-0214-1/")
if (! dir.exists("./04_model")){
  dir.create("./04_model")
}
setwd("./04_model")

library(magrittr)
library(tidyverse)
library(lance)
library(ggplot2)
hub_gene <- read.delim2('../02_intersect/DEAOG.xls')
train.dat<-read.delim2("../00_rawdata/dat(GSE97537).xls", row.names = 1)  %>% lc.tableToNum
group = read.delim2("../00_rawdata/group(GSE97537).xls")
dat<-train.dat[hub_gene$symbol,group$sample]%>%t%>%as.data.frame()
dat<-merge(dat,group,by.x='row.names',by.y='sample')%>%tibble::column_to_rownames(var = 'Row.names')
table(group$group)
dat$group<-factor(dat$group,levels = c('CIRI','control'))
# 随机森林-------
library(randomForest)
set.seed(22)
#14
res.rf <- randomForest(group ~ ., data = dat,ntrees=500,mtry = 2,importance = T)
plot(res.rf, main = NULL)
# png("01.RF.ntree.png", width = 4, height = 4, bg = "white", units = "in", res = 300, family = "Times")
# plot(res.rf, main = NULL)
# dev.off()
# pdf("01.RF.ntree.pdf", width = 4, height = 4, family = "Times")
# plot(res.rf, main = NULL)
# dev.off()
saveRDS(res.rf, "rf.models.rds")

ntree = which.min(res.rf$err.rate[,1])
res.rf = randomForest::randomForest(group ~ ., data = dat, ntree = ntree)
#变量重要性评分
importance(res.rf,type = 1)
# mean decrease gini 计算每个变量对分类树每个节点上观测值的异质性的影响，从而比较变量的重要性。该值越大表示该变量的重要性越大
#重要性绘图
pdf(file = '01.RF.importance.pdf',w=5,h=6.5)
varImpPlot(res.rf,main = 'importance')
dev.off()
png(file = '01.RF.importance.png',w=400,h=500)
varImpPlot(res.rf,main = 'importance')
dev.off()
importance <- data.frame(res.rf$importance)%>%rownames_to_column(var='symbol')
importance <- importance[order(importance$MeanDecreaseGini,decreasing = T),]

importance <- importance[which(importance$MeanDecreaseGini>0),]
write.table(importance,'rf.importance.xls',sep = '\t',row.names = F,quote = F)

# XGBoost---------
# setwd('D:/Project_file/XA-0214-1/xgboost/')
# #install.packages('Ckmeans.1d.dp')
# library(Ckmeans.1d.dp)
# #BiocManager::install('lance')
# library(xgboost)
# library(magrittr)
# #library(lance)
# library(ggplot2)
# hub_gene <- read.delim2('DEAOG.xls')
# train.dat<-read.delim2("dat(GSE97537).xls", row.names = 1)
# group = read.delim2("group(GSE97537).xls")
# genes <- rownames(train.dat)
# train.dat <- as.data.frame(lapply(train.dat, as.numeric))
# rownames(train.dat) <- genes
# 
# dat<-train.dat[hub_gene$symbol,]%>%t%>%as.data.frame()
# dat<-merge(dat,group,by.x='row.names',by.y='sample')%>%tibble::column_to_rownames(var = 'Row.names')
# table(dat$group)
# dat$group<-factor(dat$group,levels = c('CIRI','control'))
# 
# #BiocManager::install('caret')
# library(caret)
# set.seed(123)
# folds <- createFolds(y=dat$group,k=5)
# ### XGB模型
# for  (i in 1:5){
#   fold_test <- dat[folds[[i]],]    #取fold i数据，建立测试集
#   fold_train <- dat[-folds[[i]],]  #建立训练集
#   train_matrix <- sparse.model.matrix(group ~.-1, data = fold_train)
#   # 训练集的数据预处理
#   # 将trainset的1-8列（自变量）转换为矩阵
#   traindata1<- data.matrix(fold_train [,c(1:38)])
#   # 利用Matrix函数，将sparse参数设置为TRUE，转化为稀疏矩阵
#   traindata2 <- Matrix(traindata1,sparse = T)
#   # 将因变量转换为numeric类型，-1是为了从0开始计数
#   train_y <- as.numeric(fold_train[,39])-1
#   # 将自变量和因变量拼接为list
#   traindata <- list(data=traindata2,label=train_y)
#   dtrain <- xgb.DMatrix(data = traindata$data, label = traindata$label)
#   # 测试集的数据预处理
#   # 将trainset的1-10列（自变量）转换为矩阵
#   testdata1<- data.matrix(fold_test [,c(1:38)])
#   # 利用Matrix函数，将sparse参数设置为TRUE，转化为稀疏矩阵
#   testdata2 <- Matrix(testdata1,sparse = T)
#   # 将因变量转换为numeric类型，-1是为了从0开始计数
#   test_y <- as.numeric(fold_test[,39])-1
#   # 将自变量和因变量拼接为list
#   testdata <- list(data=testdata2,label=test_y)
#   dtest <- xgb.DMatrix(data = testdata $data, label = testdata $label)
#   
#   set.seed(123456)
#   res_xgb <- xgboost(data = dtrain,max_depth=2, eta=0.2,
#                      objective='binary:logistic', nround=25)
#   importance <- xgb.importance(train_matrix@Dimnames[[2]], model = res_xgb)    ##特征重要度
#   head(importance)
#   # xgb.ggplot.importance(importance)
#   predicted=predict(res_xgb ,newdata = dtest)
#   fold_predict<-ifelse(predicted >0.5, "CIRI", "control") %>% as.factor()
#   index<-(confusionMatrix(fold_test$group, fold_predict)$overall)[c('Kappa','Accuracy')]
# }
# pdf('01.importance.pdf',w=6,h=6)
# xgb.ggplot.importance(importance)
# dev.off()
# png('01.importance.png',w=400,h=400)
# xgb.ggplot.importance(importance)
# dev.off()
# # median(importance$Importance)
# write.table(importance,file = '01.importance.xls',sep = '\t',row.names = F,quote = F)
# xgbgene <- importance
# # xgbgene <- importance[which(importance$Importance>median(importance$Importance)),]
# write.table(xgbgene,file = 'xgbgene.xls',sep = '\t',row.names = F,quote = F)
# 

## 取交集--------
xgbgene <- read.delim2('../04_model/XGBoost.xls')
hubgene <- data.frame(symbol=intersect(importance$symbol,xgbgene$Feature))
library(ggvenn)
mydata<-list('randomForest'=importance$symbol,'XGBoost'=xgbgene$Feature)
pdf('03.hubgene.pdf',w=6,h=6)
ggvenn(mydata,c('randomForest','XGBoost'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()
png('03.hubgene.png',w=500,h=500)
ggvenn(mydata,c('randomForest','XGBoost'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()


write.table(hubgene,file = 'hubgene.xls',sep = '\t',row.names = F,quote = F)
