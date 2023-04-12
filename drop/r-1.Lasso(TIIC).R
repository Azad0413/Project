rm(list = ls())
setwd("/data/nas1/luchunlin/project/JNZK-214-8/")
if (! dir.exists("./02_Lasso(TIIC)")){
  dir.create("./02_Lasso(TIIC)")
}
setwd("./02_Lasso(TIIC)")
dat<-read.delim2('../01_ssGSEA/ssgsea_result.xls',row.names = 1)%>%lc.tableToNum()
#dat<-read.delim2('../01_CIBERSORT/cibersort_result.xls',row.names = 1)%>%lc.tableToNum()
colnames(dat)<-gsub('.','-',colnames(dat),fixed = T)
train_data<-t(dat)%>%as.data.frame()
survival<-read.delim2('/data/nas1/luchunlin/TCGA_survival/TCGA-COAD.survival.tsv')
train_data$sample<-rownames(train_data)
train_data<-merge(survival,train_data,by="sample")
train_data<-column_to_rownames(train_data,var = 'sample')%>%select(-'X_PATIENT')
setwd('../00_rawdata/')
write.table(train_data,file = 'train_data(TIIC).xls',sep = '\t',row.names = T,quote = F)
setwd('../02_Lasso(TIIC)/')
library(glmnet)
x_all <- subset(train_data, select = -c(OS, OS.time))
y_all <- subset(train_data, select = c(OS, OS.time))
# 拟合模型
library(survival)
fit <- glmnet(as.matrix(x_all), Surv(y_all$OS.time,y_all$OS), 
              family = "cox") 
#dev.new()
png(filename = "01.lasso_model.png", height = 400, width = 500)
plot(fit, xvar = "lambda",label = TRUE, las=1)
dev.off()
pdf(file = "01.lasso_model.pdf", height = 5)
plot(fit, xvar = "lambda",label = TRUE, las=1)
dev.off()
# 交叉验证拟合模型
set.seed(4)
#2 4
cvfit = cv.glmnet(as.matrix(x_all),
                  Surv(y_all$OS.time,y_all$OS),nfold=10,
                  family = "cox") 

png(filename = "02.lasso_verify.png", height = 400, width = 500)
plot(cvfit, las =1)
dev.off()
pdf(file = "02.lasso_verify.pdf", height = 5)
plot(cvfit, las =1)
dev.off()
# 提取指定lambda时特征的系数
coef.min = coef(cvfit, s = "lambda.min")  ## lambda.min & lambda.1se 取一个
cvfit$lambda.min
## 0.005360113
active.min = which(coef.min@i != 0)
coef.min
lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i+1]
lasso_geneids

write(lasso_geneids, "lasso_genes.csv")

