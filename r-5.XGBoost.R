rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-386-10/")
if (! dir.exists("./05_XGBoost")){
  dir.create("./05_XGBoost")
}
setwd("./05_XGBoost")
#library(Ckmeans.1d.dp)
#BiocManager::install('lance')
library(xgboost)
library(magrittr)
#library(lance)
library(ggplot2)
hub_gene <- read.delim2('01.intersect.xls')
train.dat<-read.delim2("dat.fpkm.xls", row.names = 1)
group = read.delim2("group.xls")
genes <- rownames(train.dat)
train.dat <- as.data.frame(lapply(train.dat, as.numeric))
rownames(train.dat) <- genes

colnames(train.dat) <- gsub('.','-',colnames(train.dat),fixed = T)
dat<-train.dat[hub_gene$symbol,]%>%t%>%as.data.frame()
dat<-merge(dat,group,by.x='row.names',by.y='sample')%>%tibble::column_to_rownames(var = 'Row.names')
dat$group<-factor(dat$group,levels = c('Tumor','Normal'))
colnames(dat) <- gsub('-','_',colnames(dat),fixed = T)
colnames(dat) <- gsub(' ','_',colnames(dat),fixed = T)
library("Matrix")
#BiocManager::install('caret')
library(caret)
set.seed(1)
folds <- createFolds(y=dat$group,k=5)
### XGB模型
for  (i in 1:5){
  fold_test <- dat[folds[[i]],]    #取fold i数据，建立测试集
  fold_train <- dat[-folds[[i]],]  #建立训练集
  train_matrix <- sparse.model.matrix(group ~.-1, data = fold_train)
  # 训练集的数据预处理
  # 将trainset的1-8列（自变量）转换为矩阵
  traindata1<- data.matrix(fold_train [,c(1:156)])
  # 利用Matrix函数，将sparse参数设置为TRUE，转化为稀疏矩阵
  traindata2 <- Matrix(traindata1,sparse = T)
  # 将因变量转换为numeric类型，-1是为了从0开始计数
  train_y <- as.numeric(fold_train[,157])-1
  # 将自变量和因变量拼接为list
  traindata <- list(data=traindata2,label=train_y)
  dtrain <- xgb.DMatrix(data = traindata$data, label = traindata$label)
  # 测试集的数据预处理
  # 将trainset的1-10列（自变量）转换为矩阵
  testdata1<- data.matrix(fold_test [,c(1:156)])
  # 利用Matrix函数，将sparse参数设置为TRUE，转化为稀疏矩阵
  testdata2 <- Matrix(testdata1,sparse = T)
  # 将因变量转换为numeric类型，-1是为了从0开始计数
  test_y <- as.numeric(fold_test[,157])-1
  # 将自变量和因变量拼接为list
  testdata <- list(data=testdata2,label=test_y)
  dtest <- xgb.DMatrix(data = testdata $data, label = testdata $label)
  
  set.seed(123)
  res_xgb <- xgboost(data = dtrain,max_depth=6, eta=0.5,
                     objective='binary:logistic', nround=25)
  importance <- xgb.importance(train_matrix@Dimnames[[2]], model = res_xgb)    ##特征重要度
  head(importance)
  # xgb.ggplot.importance(importance)
  predicted=predict(res_xgb ,newdata = dtest)
  fold_predict<-ifelse(predicted >0.5, "Tumor", "Normal") %>% as.factor()
  index<-(confusionMatrix(fold_test$group, fold_predict)$overall)[c('Kappa','Accuracy')]
}
pdf('01.importance.pdf',w=6,h=7)
xgb.ggplot.importance(importance)
dev.off()
png('01.importance.png',w=400,h=500)
xgb.ggplot.importance(importance)
dev.off()
median(importance$Importance)
write.table(importance,file = '01.importance.xls',sep = '\t',row.names = F,quote = F)
xgbgene <- importance[which(importance$Importance>median(importance$Importance)),]
write.table(xgbgene,file = 'xgbgene.xls',sep = '\t',row.names = F,quote = F)
