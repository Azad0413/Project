rm(list = ls())
setwd("/data/nas1/luchunlin/project/JNZK-204(modify)/")
if (! dir.exists("./07_cindex")){
  dir.create("./07_cindex")
}
setwd("./07_cindex")

##构建COX模型，绘制列线图---------
train_phenotype<-read.delim2('../06_clinical/phenotype.xls')
train_phenotype$OS<-as.numeric(train_phenotype$OS)
train_phenotype$OS.time<-as.numeric(train_phenotype$OS.time)
train_phenotype2<-train_phenotype

table(train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('a','',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('b','',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('TX',NA,train_phenotype2$T.stage)
train_phenotype2[train_phenotype2==''] <- NA
train_phenotype2$T.stage <- gsub('T1','T1/2',train_phenotype2$T.stage)
train_phenotype2$T.stage <- gsub('T2','T1/2',train_phenotype2$T.stage)
train_phenotype2$T.stage <- gsub('T3','T3/4',train_phenotype2$T.stage)
train_phenotype2$T.stage <- gsub('T4','T3/4',train_phenotype2$T.stage)

table(train_phenotype2$Grade)
table(train_phenotype2$ATI)

table(train_phenotype2$Gender)

library(tidyverse)
library(lance)
colnames(train_phenotype2)
colnames(train_phenotype2)<-c('id','T.stage','Grade','Gender','ATI','OS','OS.time')
risk<-read.delim2('../06_clinical/risk.xls')%>%lc.tableToNum()
sub_risk <- subset(risk, select = c(sample, risk))
colnames(sub_risk) <- c('id','riskScore')
sub_risk$id <- gsub('.','-',sub_risk$id,fixed = T)
train_risk_clinical <- merge(train_phenotype2,
                             sub_risk,
                             by = "id")
rownames(train_risk_clinical) <- train_risk_clinical$id
train_risk_clinical = subset(train_risk_clinical, select = -c(id))
dim(train_risk_clinical)
colnames(train_risk_clinical)
library(survival)
multi_cov<-c('riskScore',"T.stage")
cox_data_prog <- as.formula(paste0('Surv(OS.time, OS)~',
                                   paste(multi_cov,
                                         sep = '',
                                         collapse = '+')))
cox_more_prog <- coxph(cox_data_prog,
                       data = as.data.frame(train_risk_clinical))

#C指数即一致性指数，用来评价模型的预测能力。c指数是指所有病人对子中预测结果与实际结果一致的对子所占的比例。
C_index <- cox_more_prog$concordance['concordance']
if(C_index >= 0.9){
  print("High accuracy")
}else{
  if(C_index < 0.9 & C_index >= 0.7){
    print("Medium accuracy")
  }else{
    print("Low accuracy")
  }
}
sum.surv<-summary(cox_more_prog)
c_index<-sum.surv$concordance
sum.surv$conf.int
sum.surv$coefficients
sum.surv$concordance
c_index
# C      se(C) 
# 0.71053211 0.03151625 

LOW <- 0.71053211-0.03151625
High <- 0.71053211+0.03151625
