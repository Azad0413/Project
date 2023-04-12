rm(list = ls())
setwd('/data/nas1/luchunlin/project/BJTC-158(返修）/')

library(glmnet)
library(foreign)
library(rms)
library(pROC)

### 8个正常样本 34个OSAyangb
input = read.table("input.txt",sep="\t",header=T,check.names=F) 

input$sample <- factor(input$sample)
###两个数据矩阵，Y是结果，X是数据的变量,Y结果中 1表示发生OSA，0表示正常
y<-data.matrix(input[,1])

x<-as.matrix(input[,c(2:18)])

f1 =glmnet(x, y, family="binomial", nlambda=100, alpha=1)
#这里alpha=1为LASSO回归，如果等于0就是岭回归

print(f1)

###可以看到随着lambdas增加，自由度和残差减少，最小lambda为0.000021
pdf(file="lambda.pdf",width = 8,height = 6)
plot(f1,xvar="lambda", label=TRUE)
dev.off()
###横坐标为随着lambdas的对数，纵坐标为变量系数，
#可以看到随着lambdas增加变量系数不断减少，部分变量系数变为0(等于没有这个变量了)

###交叉验证（这部也可以不作）
predict(f1,newx=x[2:5,], type = "response")
###将Y的数据强制转换
y <- as.numeric(unlist(y))
cvfit=cv.glmnet(x,y)
pdf(file="cvfit.pdf",width = 8,height = 6)
plot(cvfit)
dev.off()
#我们这个图中有两条虚线，一个是均方误差最小时的λ值，
#一个是距离均方误差最小时一个标准误的λ值，
#有点拗口没关系，我们只要知道它是多少就可以了

cvfit$lambda.min#求出最小值
#[1] 0.06801549

# ##[1] 0.08991247
cvfit$lambda.1se#求出最小值一个标准误的λ值
#[1] 0.2077095

##我们得出这两个值后分别带进模型看一看
l.coef2<-coef(cvfit$glmnet.fit,s=0.06801549,exact= F)
l.coef1<-coef(cvfit$glmnet.fit,s=0.2077095,exact= F)

l.coef1##模型没有变量
l.coef2###还剩4个变量AREG 、ZFP36 、ATF3、DUSP1 
##14
# set.seed(14)
# dat <- data.frame(input[,c('AREG','ATF3','DUSP1','ZFP36')])
# dat$group<-ifelse(input$sample==0,'Normal','OSA')
# res.lasso <- cv.glmnet(as.matrix(dat[-ncol(dat)]), dat$group, family = "binomial", 
#                        type.measure = "auc")
# plot(res.lasso)
# plot(res.lasso$glmnet.fit, xvar = 'lambda')
# res.lasso = glmnet(as.matrix(dat[-ncol(dat)]), dat$group, family = "binomial",  lambda = res.lasso$lambda.min)
# df.pred = predict.glmnet(res.lasso, newx = as.matrix(dat[-ncol(dat)]), type = "link")[,1]
# df.pred = data.frame(sample = rownames(dat), score = df.pred)
# df.pred$group = input$sample
# #df.pred$group<-ifelse(df.pred$group==0,'Normal','OSA')
# 
# cs.en = confusion.glmnet(res.lasso, newx = as.matrix(dat[-ncol(dat)]), newy = dat$group) %>% as.matrix()
# cs.en = cs.en[2:1,2:1]
# cs.en
# cs.en <- cs.en %>% t %>% as.data.frame()
# 
# # 真阳性放在第一位
# cs.en$True <- factor(cs.en$True, levels = c('Normal','OSA'))
# p1 <- ggplot(cs.en, aes(x = Predicted, y = True, fill = Freq)) +
#   geom_tile(color = "grey") +
#   geom_text(aes(label = Freq), family = "Times", size = 6) +
#   scale_fill_gradient(low = "white", high = "royalblue") +
#   coord_fixed() +
#   labs(x="Predicted Classification", y="True Classification") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(size = 13, face = "bold", family = "Times"),
#         axis.text.y = element_text(size = 13, face = "bold", family = "Times"),
#         axis.title = element_text(size = 15, face = "bold", family = "Times"),
#         legend.position = "none",
#         panel.grid = element_blank())
# p1
## test------
# library(readxl)
# test<-read_xlsx('BJTC-158_qPCR.xlsx')
# test<-column_to_rownames(test,var = 'sample')
# df.pred = predict.glmnet(res.lasso, newx = as.matrix(test[-ncol(test)]), type = "link")[,1]
# df.pred = data.frame(sample = rownames(test), score = df.pred)
# df.pred$group = test$group
# cs.en = confusion.glmnet(res.lasso, newx = as.matrix(test[-ncol(test)]), newy = test$group) %>% as.matrix()
# cs.en
# cs.en = cs.en[2:1,2:1]
# cs.en <- cs.en %>% t %>% as.data.frame()
# 
# # 真阳性放在第一位
# cs.en$True <- factor(cs.en$True, levels = c('Normal','OSA'))
# p1 <- ggplot(cs.en, aes(x = Predicted, y = True, fill = Freq)) +
#   geom_tile(color = "grey") +
#   geom_text(aes(label = Freq), family = "Times", size = 6) +
#   scale_fill_gradient(low = "white", high = "royalblue") +
#   coord_fixed() +
#   labs(x="Predicted Classification", y="True Classification") +
#   theme_minimal() +
#   theme(axis.text.x = element_text(size = 13, face = "bold", family = "Times"),
#         axis.text.y = element_text(size = 13, face = "bold", family = "Times"),
#         axis.title = element_text(size = 15, face = "bold", family = "Times"),
#         legend.position = "none",
#         panel.grid = element_blank())
# p1
# 
# 
# library(pROC)
# lasso_roc <- roc(df.pred$group, df.pred$score)
# ggroc(lasso_roc,color = "red",
#       linetype = 1,
#       size = 1,
#       alpha = 1,
#       legacy.axes = T)+
#   geom_abline(intercept = 0,
#               slope = 1,
#               color = "grey",
#               size = 1,
#               linetype = 1)+
#   labs(x = "False Postive Rate(1 - Specificity)",
#        y = "True Positive Rate(Sensivity or Recall)",
#        title = "ROC curve")+
#   annotate("text",x = 0.70,y = 0.30,
#            label = paste("AUC =",signif(auc(lasso_roc),2)),
#            size = 5,family = "Times")+
#   theme_bw()+
#   theme(panel.background = element_rect(fill = "transparent"),
#         panel.grid = element_blank(),
#         axis.ticks.length = unit(0.4,"lines"),
#         axis.ticks = element_line(color = "black"),
#         axis.line = element_line(size = 0.5,colour = "black"),
#         axis.title = element_text(colour = "black",size = 15),
#         axis.text = element_text(colour = "black",size = 10),
#         plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
#         text = element_text(size = 8,color = "black",family = "Times"))

### lasso -logsitic----------
# rm(list = ls())
# 
# ### 8个正常样本 34个OSAyangb
# input = read.table("input.txt",sep="\t",header=T,check.names=F) 
# 
# input$sample <- factor(input$sample)
# ###两个数据矩阵，Y是结果，X是数据的变量,Y结果中 1表示发生OSA，0表示正常
# set.seed(1)
# y<-data.matrix(input[,1])
# 
# x<-as.matrix(input[,c(2:18)])
# 
# f1 =glmnet(x, y, family="binomial", nlambda=100, alpha=1)
# #这里alpha=1为LASSO回归，如果等于0就是岭回归
# 
# print(f1)
# 
# ###可以看到随着lambdas增加，自由度和残差减少，最小lambda为0.000021
# pdf(file="lambda.pdf",width = 8,height = 6)
# plot(f1,xvar="lambda", label=TRUE)
# dev.off()
# ###横坐标为随着lambdas的对数，纵坐标为变量系数，
# #可以看到随着lambdas增加变量系数不断减少，部分变量系数变为0(等于没有这个变量了)
# 
# ###交叉验证（这部也可以不作）
# predict(f1,newx=x[2:5,], type = "response")
# ###将Y的数据强制转换
# y <- as.numeric(unlist(y))
# cvfit=cv.glmnet(x,y)
# pdf(file="cvfit.pdf",width = 8,height = 6)
# plot(cvfit)
# dev.off()
# #我们这个图中有两条虚线，一个是均方误差最小时的λ值，
# #一个是距离均方误差最小时一个标准误的λ值，
# #有点拗口没关系，我们只要知道它是多少就可以了
# 
# cvfit$lambda.min#求出最小值
# #[1] 0.06801549
# 
# # ##[1] 0.08991247
# cvfit$lambda.1se#求出最小值一个标准误的λ值
# #[1] 0.2077095
# 
# ##我们得出这两个值后分别带进模型看一看
# l.coef2<-coef(cvfit$glmnet.fit,s=0.06801549,exact= F)
# l.coef1<-coef(cvfit$glmnet.fit,s=0.2077095,exact= F)
# 
# l.coef1##模型没有变量
# l.coef2###还剩4个变量AREG 、ZFP36 、ATF3、DUSP1 


mod<-glm(sample~AREG+ZFP36+ATF3+DUSP1,family="binomial",data= input)
summary(mod)

mod1<-lrm(sample~AREG+ZFP36+ATF3+DUSP1,data= input,x = T,y = T)
cal1 <- calibrate(mod1, cmethod='hare', method='boot', B=1000,data=vad)#建模组中绘制校准曲线
plot(cal1,xlim=c(0,1.0),ylim=c(0,1.0))#打印出校准曲线

###求OR和95%CI
exp(confint(mod))
#数据展示
#2.5 %       97.5 %
#  (Intercept) 545.23647619 3.047376e+15
#AREG          0.05940385 6.062758e+00
#ZFP36         0.09061539 2.620537e+01
#ATF3          0.02596174 1.436577e+01
#DUSP1         0.01461241 3.800547e+00
#

exp(coef(mod))
#
#(Intercept)         AREG        ZFP36         ATF3        DUSP1 
#3.753084e+07 5.948813e-01 1.638835e+00 6.464067e-01 2.931767e-01 
#
#查看预测结果
predicidata <- data.frame(input[,c('AREG','ATF3','DUSP1','ZFP36')])
predicidata$y=(predict(mod,predicidata,type="response"))
#predicidata$actual<-input$sample
# predicidata$pred<-ifelse(predicidata$y>0.5,1,0)
predicted = ifelse(predicidata$y >0.5, "OSA", "Normal") %>% as.factor()
input$group<-ifelse(input$sample=='0','Normal','OSA')
input$group<-factor(input$group,levels = c('Normal','OSA'))
cs.lg = confusionMatrix(input$group, predicted)$table
cs.lg



y.pred <- predict(mod,predicidata,type="response")
pdf(file="核密度估计图.pdf",width = 8,height = 6)
plot(density(y.pred),main='预测患病可能性的核密度估计图',xlab='预测患病可能性',y='密度')
dev.off()
predicidata<-as.matrix(predicidata)

library(ggplot2)
predicted.data <- data.frame(probability.of.hd=mod$fitted.values,sample=input$sample)
str(predicted.data)
predicted.data <- predicted.data[order(predicted.data$probability.of.hd, decreasing=FALSE),]
predicted.data$rank <- 1:nrow(predicted.data)
pdf(file="可视化结果.pdf",width = 8,height = 6)
ggplot(data=predicted.data, aes(x=rank, y=probability.of.hd)) +
  geom_point(aes(color=sample), alpha=1, shape=4, stroke=2) +
  xlab("Index") +
  ylab("Predicted probability of getting heart disease")
dev.off()


###AUC###
library(pROC)
def_pred<-predict(mod, newdata=input)##生成概率
roc<-multiclass.roc (as.ordered(input$sample) ,as.ordered(def_pred))#拟合ROC
roc1<-roc(as.ordered(input$sample) ,as.ordered(def_pred))


round(auc(roc1),3)##AUC
round(ci(roc1),3)##95%CI

plot(roc1)
pdf(file="AUC.pdf",width = 8,height = 6)
plot(roc1,print.auc=T, auc.polygon=T, grid=c(0.1, 0.2), grid.col=c("green","red"), 
     max.auc.polygon=T, auc.polygon.col="skyblue",print.thres=T)
dev.off()



ggroc1 <- ggroc(roc1, # roc()函数创建的对象
                alpha = 0.5, # 设置曲线透明度
                colour = "red",  # 设置曲线颜色
                linetype = 1, size = 1, # 设置曲线线型和大小 
                legacy.axes = TRUE)  # 坐标轴为"1-specificity"
##修改背景网格线
ggroc2 <- ggroc1+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  annotate("text",x=0.6,y=0.70,label = "AUC = 0.842",
           size = 5,fontface = "bold",family = "A")
ggroc2





###ggplot####ROC曲线绘制

library(pROC)

fit <- glm(sample ~AREG+ZFP36+ATF3+DUSP1,family="binomial",data= input)

roc1 <-roc(input$sample, as.vector(fitted.values(fit)), 
           percent=F,boot.n=1000, ci.alpha=0.9, 
           stratified=FALSE, plot=TRUE, grid=TRUE, 
           show.thres=TRUE, legacy.axes = TRUE, reuse.auc = TRUE,
           print.auc = TRUE, print.thres.col = "blue",
           ci=TRUE, ci.type="bars", print.thres.cex = 0.7,
           main = paste("ROC curve using","(N = ",nrow(input),")") )

plot(roc1$FP, roc1$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="False positive rate", ylab="True positive rate",
     main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)



###
score <-a
aupr=AUC(obs=score$link,pred=score$score,curve = "PR", simplif=TRUE, main = "PR curve")

##  混淆矩阵热图------

z = c(30,4,4,4)

z = matrix(z, ncol=2) 
colnames(z) = c("OSA","Normal") 
rownames(z) = c("OSA","Normal") 

image(z[,ncol(z):1], axes=FALSE) 


axis(2, at = seq(0, 1, length=length(colnames(z))), labels=colnames(z)) 

head(z)

library(ComplexHeatmap)

pdf(file="2.pdf",width = 8,height = 5)
Heatmap(z,name = "Random forest ",
        col = c("turquoise","blue","dodgerblue"),
        cluster_rows=FALSE,cluster_columns = FALSE,
        row_names_side = c("left"),
        column_names_side = c("top"),
        column_names_rot = 0)
dev.off()


### 验证集 ---------
library(readxl)
test<-read_xlsx('BJTC-158_qPCR.xlsx')
test<-column_to_rownames(test,var = 'sample')
test$sample<-ifelse(test$group=='Normal',0,1)

colnames(test)
predicidata <- data.frame(test[,c("ZFP36","AREG","DUSP1","ATF3")])
predicidata$y=(predict(mod,predicidata,type="response"))
predicted = ifelse(predicidata$y >0.5, "OSA", "Normal") %>% as.factor()
predicted
test$group<-factor(test$group,levels=c('Normal','OSA'))
cs.lg = confusionMatrix(test$group, predicted)$table
cs.lg
library(pROC)
def_pred<-predict(mod, newdata=test,type = 'response')##生成概率
def_pred
write.table(def_pred,'predict_test.xls',sep = '\t',row.names = F,quote = F)
roc<-multiclass.roc (as.ordered(test$sample) ,as.ordered(def_pred))#拟合ROC
roc1<-roc(as.ordered(test$sample) ,as.ordered(def_pred))
round(auc(roc1),3)##AUC
round(ci(roc1),3)##95%CI

plot(roc1)
pdf(file="AUC.pdf",width = 8,height = 6)
plot(roc1,print.auc=T, auc.polygon=T, grid=c(0.1, 0.2), grid.col=c("green","red"), 
     max.auc.polygon=T, auc.polygon.col="skyblue",print.thres=T)
dev.off()



ggroc1 <- ggroc(roc1, # roc()函数创建的对象
                alpha = 0.5, # 设置曲线透明度
                colour = "red",  # 设置曲线颜色
                linetype = 1, size = 1, # 设置曲线线型和大小 
                legacy.axes = TRUE)  # 坐标轴为"1-specificity"
##修改背景网格线
ggroc2 <- ggroc1+
  theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())+
  annotate("text",x=0.6,y=0.70,label = "AUC = 0.842",
           size = 5,fontface = "bold",family = "A")
ggroc2





