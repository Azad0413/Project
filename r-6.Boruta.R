rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-386-10/")
if (! dir.exists("./06_Boruta")){
  dir.create("./06_Boruta")
}
setwd("./06_Boruta")

#install.packages("Boruta")
library(Boruta)
hub_gene <- read.delim2('../04_DEm7Glnc/01.intersect.xls')
train.dat<-read.delim2("../00_rawdata/dat.fpkm.xls", row.names = 1)%>%lc.tableToNum()
group = read.delim2("../01_WGCNA/group.xls")
colnames(train.dat) <- gsub('.','-',colnames(train.dat),fixed = T)

dat<-train.dat[hub_gene$symbol,]%>%t%>%as.data.frame()
dat<-merge(dat,group,by.x='row.names',by.y='sample')%>%tibble::column_to_rownames(var = 'Row.names')
dat$group<-factor(dat$group,levels = c('Tumor','Normal'))
colnames(dat) <- gsub('-','_',colnames(dat),fixed = T)
colnames(dat) <- gsub(' ','_',colnames(dat),fixed = T)

set.seed(123)
#boruta.train <- Boruta(x=dat[-ncol(dat)],y=dat$group, doTrace = 2,pValue = 0.05,mcAdj = T,maxRuns = 300)
boruta.train <- Boruta(group~., data = dat, doTrace = 2)
print(boruta.train)
# 38 attributes confirmed important: AC005082.12, AC096574.5, AC144831.1, AL161668.5, AP000696.2 and 33 more;
# 1 attributes confirmed unimportant: RP1_265C24.8;
# 3 tentative attributes left: RP11_37L2.1, RP11_77I22.2, RP5_943J3.2;
png(file = '01.importance.png',w=600,h=400)
plot(boruta.train, xlab = "", xaxt = "n")
lz<-lapply(1:ncol(boruta.train$ImpHistory),function(i)
  
  boruta.train$ImpHistory[is.finite(boruta.train$ImpHistory[,i]),i])

names(lz) <- colnames(boruta.train$ImpHistory)

Labels <- sort(sapply(lz,median))

axis(side = 1,las=2,labels = names(Labels),
       
       at = 1:ncol(boruta.train$ImpHistory), cex.axis = 0.6)
dev.off()
pdf(file = '01.importance.pdf',w=7,h=5)
plot(boruta.train, xlab = "", xaxt = "n")
lz<-lapply(1:ncol(boruta.train$ImpHistory),function(i)
  
  boruta.train$ImpHistory[is.finite(boruta.train$ImpHistory[,i]),i])

names(lz) <- colnames(boruta.train$ImpHistory)

Labels <- sort(sapply(lz,median))

axis(side = 1,las=2,labels = names(Labels),
     
     at = 1:ncol(boruta.train$ImpHistory), cex.axis = 0.6)
dev.off()


# 蓝色的盒状图对应一个阴影属性的最小、平均和最大Z分数。红色、黄色和绿色的盒状图分别代表拒绝、暂定和确认属性的Z分数。
#现在我们对实验性属性进行判定。实验性属性将通过比较属性的Z分数中位数和最佳阴影属性的Z分数中位数被归类为确认或拒绝。
library(dplyr)
boruta.imp <- function(x){
  imp <- reshape2::melt(x$ImpHistory, na.rm=T)[,-1]
  colnames(imp) <- c("Variable","Importance")
  imp <- imp[is.finite(imp$Importance),]
  
  variableGrp <- data.frame(Variable=names(x$finalDecision), 
                            finalDecision=x$finalDecision)
  
  showGrp <- data.frame(Variable=c("shadowMax", "shadowMean", "shadowMin"),
                        finalDecision=c("shadowMax", "shadowMean", "shadowMin"))
  
  variableGrp <- rbind(variableGrp, showGrp)
  
  boruta.variable.imp <- merge(imp, variableGrp, all.x=T)
  
  sortedVariable <- boruta.variable.imp %>% group_by(Variable) %>% 
    summarise(median=median(Importance)) %>% arrange(median)
  sortedVariable <- as.vector(sortedVariable$Variable)
  
  
  boruta.variable.imp$Variable <- factor(boruta.variable.imp$Variable, levels=sortedVariable)
  
  invisible(boruta.variable.imp)
}


final.boruta <- TentativeRoughFix(boruta.train)


print(final.boruta)
getSelectedAttributes(final.boruta, withTentative = F)
boruta.df <- attStats(final.boruta)
class(boruta.df)
print(boruta.df)

#提取重要的变量和可能重要的变量

boruta.finalVarsWithTentative <- data.frame(Item=getSelectedAttributes(final.boruta, withTentative = T), 
                                            Type="Boruta_with_tentative")

# 交叉验证选择参数并拟合模型,定义一个函数生成一些列用来测试的mtry (一系列不大于总变量数的数值)。

generateTestVariableSet <- function(num_toal_variable){
  max_power <- ceiling(log10(num_toal_variable))
  tmp_subset <- c(unlist(sapply(1:max_power, function(x) (1:10)^x, simplify = F)), ceiling(max_power/3))
  #return(tmp_subset)
  base::unique(sort(tmp_subset[tmp_subset<num_toal_variable]))
}
# generateTestVariableSet(78)
# 选择关键特征变量相关的数据
# 提取训练集的特征变量子集
boruta_train_data <- dat[, boruta.finalVarsWithTentative$Item]
boruta_mtry <- generateTestVariableSet(ncol(boruta_train_data))
# 使用 Caret 进行调参和建模

library(caret)
# Create model with default parameters
trControl <- trainControl(method="repeatedcv", number=10, repeats=5)


# 根据经验或感觉设置一些待查询的参数和参数值
tuneGrid <- expand.grid(mtry=boruta_mtry)
train_data_group <- as.factor(dat$group)
borutaConfirmed_rf_default <- train(x=boruta_train_data, y = train_data_group, method="rf", 
                                    tuneGrid = tuneGrid, # 
                                    metric="rf", #metric='Kappa';Accuracy;rf
                                    trControl=trControl)
borutaConfirmed_rf_default

# 可视化不同参数的准确性分布

plot(borutaConfirmed_rf_default)

# 可视化重要的变量
#pdf(file = paste0("02.BORUTA_DEGs.pdf"),width = 5,height = 4)
#a <- dev.cur()   #记录pdf设备
#png(file = paste0("02.BORUTA_DEGs.png"),width = 5, height=4, units="in", res=300) 
dev.control("enable")
par(mar = c(5,10,3,5));#下、左、上、右
dotPlot(varImp(borutaConfirmed_rf_default))
importance <- varImp(borutaConfirmed_rf_default)
importance <- importance$importance%>%as.data.frame()%>%rownames_to_column(var = 'symbol')
importance <- importance[order(importance$Overall,decreasing = T),]
write.table(importance,file = '01.importance.xls',sep = '\t',row.names = F,quote = F)
borutagene <- importance[c(1:10),]%>%select('symbol')
write.table(borutagene,file = '02.brtgene.xls',sep = '\t',row.names = F,quote = F)

#dev.copy(which = a)  #复制来自png设备的图片到pdf
#dev.off()
#dev.off()
# 提取最终选择的模型，并绘制 ROC 曲线评估模型
# borutaConfirmed <- borutaConfirmed_rf_default$finalModel


# # 提取confirmed基因
# temdat <- boruta.variable.imp %>% as.data.frame()
# temdat <- temdat[,c(1,3)] %>% distinct()
# table(temdat$finalDecision)
# temdat <- temdat[grep("Confirmed",temdat$finalDecision),]
# BORUTA_gene <- temdat$Variable %>% as.character()

borutagene <- boruta.df[which(boruta.df$decision=='Confirmed'),]
write.table(boruta.df,file = '01.boruta.df.xls',sep = '\t',row.names = T,quote = F)
write.table(borutagene,file = '02.brtgene.xls',sep = '\t',row.names = T,quote = F)

## 
