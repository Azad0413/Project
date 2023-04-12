rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-370-8/")
if (! dir.exists("./06_XGBoost")){
  dir.create("./06_XGBoost")
}
setwd("./06_XGBoost")
library(magrittr)
library(ggplot2)
library(tidyverse)

hub_gene <- read.delim2('../03_DEGs/intersect.xls')
train.dat<-read.delim2("../03_DEGs/dat_final.xls", row.names = 1) 
group = read.delim2("../03_DEGs/group.xls")
genes <- rownames(train.dat)
train.dat <- as.data.frame(lapply(train.dat, as.numeric))
rownames(train.dat) <- genes
dat<-train.dat[hub_gene$.,]%>%t%>%as.data.frame()
dat<-merge(dat,group,by.x='row.names',by.y='sample')%>%tibble::column_to_rownames(var = 'Row.names')
dat$group<-factor(dat$group,levels = c('RB','control'))
library("Matrix")
#BiocManager::install('caret')
library(caret)
set.seed(123)
folds <- createFolds(y=dat$group,k=5)
### XGB模型
for  (i in 1:5){
  fold_test <- dat[folds[[i]],]    #取fold i数据，建立测试集
  fold_train <- dat[-folds[[i]],]  #建立训练集
  train_matrix <- sparse.model.matrix(group ~.-1, data = fold_train)
  # 训练集的数据预处理
  # 将trainset的1-8列（自变量）转换为矩阵
  traindata1<- data.matrix(fold_train [,c(1:10)])
  # 利用Matrix函数，将sparse参数设置为TRUE，转化为稀疏矩阵
  traindata2 <- Matrix(traindata1,sparse = T)
  # 将因变量转换为numeric类型，-1是为了从0开始计数
  train_y <- as.numeric(fold_train[,11])-1
  # 将自变量和因变量拼接为list
  traindata <- list(data=traindata2,label=train_y)
  dtrain <- xgb.DMatrix(data = traindata$data, label = traindata$label)
  # 测试集的数据预处理
  # 将trainset的1-10列（自变量）转换为矩阵
  testdata1<- data.matrix(fold_test [,c(1:10)])
  # 利用Matrix函数，将sparse参数设置为TRUE，转化为稀疏矩阵
  testdata2 <- Matrix(testdata1,sparse = T)
  # 将因变量转换为numeric类型，-1是为了从0开始计数
  test_y <- as.numeric(fold_test[,11])-1
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
  fold_predict<-ifelse(predicted >0.5, "control", "RB") %>% as.factor()
  index<-(confusionMatrix(fold_test$group, fold_predict)$overall)[c('Kappa','Accuracy')]
}
xgb.ggplot.importance(importance)
pdf('01.importance.pdf',w=6,h=6)
xgb.ggplot.importance(importance)
dev.off()
png('01.importance.png',w=400,h=400)
xgb.ggplot.importance(importance)
dev.off()

importance <- read.delim2('importance.xls')
lassogene <- read.csv('../05_lasso/lasso_genes.csv',header = F)
hubgene <- intersect(lassogene$V1,importance$Feature)%>%as.data.frame()
colnames(hubgene) <- 'symbol'
write.table(hubgene,file = 'hubgene.xls',sep = '\t',row.names = F,quote = F)

library(ggvenn)
mydata<-list(Lasso=lassogene$V1,XGBoost=importance$Feature)
pdf('02.venn.pdf',w=6,h=6)
ggvenn(mydata,c('Lasso','XGBoost'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       text_size = 6,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()
png('02.venn.png',w=400,h=400)
ggvenn(mydata,c('Lasso','XGBoost'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       text_size = 6,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()
