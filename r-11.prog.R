rm(list = ls())
setwd("/data/nas1/luchunlin/project/JNZK-214-8/")
if (! dir.exists("./11_progmodel")){
  dir.create("./11_progmodel")
}
setwd("./11_progmodel")
##构建COX模型，绘制列线图---------
train_phenotype<-read.delim2('../10_clinical/phenotype.xls')
train_phenotype$OS<-as.numeric(train_phenotype$OS)
train_phenotype$OS.time<-as.numeric(train_phenotype$OS.time)
train_phenotype2<-train_phenotype
table(train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('a','',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('b','',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('Tis',NA,train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('T1','1/2',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('T2','1/2',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('T3','3',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('T4','4',train_phenotype2$T.stage)
table(train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('a','',train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('b','',train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('c','',train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('N0','0',train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('N1','1',train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('N2','2',train_phenotype2$N.stage)
table(train_phenotype2$M.stage)
train_phenotype2$M.stage<-gsub('a','',train_phenotype2$M.stage)
train_phenotype2$M.stage<-gsub('b','',train_phenotype2$M.stage)
train_phenotype2$M.stage<-gsub('MX',NA,train_phenotype2$M.stage)
train_phenotype2[train_phenotype2=='']<-NA 
train_phenotype2$M.stage<-gsub('M0','0',train_phenotype2$M.stage)
train_phenotype2$M.stage<-gsub('M1','1',train_phenotype2$M.stage)
table(train_phenotype2$age)
train_phenotype2$age<-ifelse(train_phenotype2$age>60,'>60','<=60')
table(train_phenotype2$gender)
#train_phenotype2$gender <- ifelse(train_phenotype2$gender=='male',1,0)
colnames(train_phenotype2)
colnames(train_phenotype2)<-c('id','Age','Gender','T.stage','N.stage','M.stage','OS','OS.time')
risk<-read.delim2('../09_risk/risk.xls')%>%lc.tableToNum()
sub_risk <- subset(risk, select = c(id, riskScore))
train_risk_clinical <- merge(train_phenotype2,
                             sub_risk,
                             by = "id")
rownames(train_risk_clinical) <- train_risk_clinical$id
train_risk_clinical = subset(train_risk_clinical, select = -c(id))
dim(train_risk_clinical)
colnames(train_risk_clinical)
library(survival)
multi_cov<-c('riskScore','Age','Gender',"T.stage","N.stage","M.stage")
cox_data_prog <- as.formula(paste0('Surv(OS.time, OS)~',
                                   paste(multi_cov,
                                         sep = '',
                                         collapse = '+')))
cox_more_prog <- coxph(cox_data_prog,
                       data = as.data.frame(train_risk_clinical))
# Nomogram
library(rms)
ddist <- datadist(train_risk_clinical)
options(datadist='ddist')

# 构建COX模型，绘制列线图

res.cox <- psm(cox_data_prog,
               data = train_risk_clinical, dist = 'lognormal')
surv <- Survival(res.cox) # 构建生存概率函数
function(x) surv(365, x) # 1年事件发生概率
function(x) surv(1095, x) # 3年事件发生概率
function(x) surv(1825, x) # 5年事件发生概率

nom.cox <- nomogram(res.cox,
                    fun = list(function(x) surv(365, x),
                               function(x) surv(1095, x),
                               function(x) surv(1825, x)),
                    funlabel=c("1-year Survival Probability", "3-year Survival Probability", "5-year Survival Probability"),
                    maxscale = 10,
                    fun.at = c(0.01,seq(0.1,0.9,by=0.2),0.95,0.99),
                    lp=F)



plot(nom.cox, cex.axis  = 1.5, cex.var = 1.7)
png(filename = "01.nomogram_line_points.png", height = 700, width = 1200)
plot(nom.cox, cex.axis  = 1.5, cex.var = 1.7)
dev.off()
pdf(file = "01.nomogram_line_points.pdf", height = 9, width = 17)
plot(nom.cox, cex.axis  = 1.5, cex.var = 1.7)
dev.off()
##07-4 构建校准曲线---------
coxm_1 <- cph(cox_data_prog,
              data=train_risk_clinical,
              surv=T,
              x=T,
              y=T,
              time.inc = 365)
cal_1<-calibrate(coxm_1,u=365,cmethod='KM',m=100)

##绘制3年生存期校曲线
##time.in 和 u 要是一样的，都是要评价的时间节点
coxm_3 <- cph(cox_data_prog,
              data=train_risk_clinical,
              surv=T,
              x=T,y=T,
              time.inc = 3*365)
cal_3 <-calibrate(coxm_3,u=3*365,cmethod='KM',m=100,B=1000)

coxm_5 <- cph(cox_data_prog,
              data=train_risk_clinical,
              surv=T,
              x=T,y=T,
              time.inc = 5*365)
cal_5 <-calibrate(coxm_5,u=5*365,cmethod='KM',m=100,B=100)

par(mar=c(7,4,4,3),cex=1.5)
plot(cal_1,
     subtitles = F,
     lwd=2,lty=1, ##设置线条形状和尺寸
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 1-5 year Survival Probability',#便签
     ylab='Actual 1-5 year Survival Probability',#标签
     col="#00468b",#设置一个颜色
     xlim = c(0,1),ylim = c(0,1)) ##x轴和y轴范围
plot(cal_3,
     add = T,
     subtitles = F,
     lwd=2,lty=1,  ##设置线条宽度和线条类型
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 1-5 year Survival Probability',#便签
     ylab='Actual 1-5 year Survival Probability',#标签
     col="#ed0000",#设置一个颜色
     xlim = c(0,1),ylim = c(0,1)) ##x轴和y轴范围
plot(cal_5,
     add = T,
     subtitles = F,
     lwd=2,lty=1, ##设置线条形状和尺寸
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 1-5 year Survival Probability',#便签
     ylab='Actual 1-year Survival Probability',#标签
     col="#42b540",#设置一个颜色
     xlim = c(0,1),ylim = c(0,1)) ##x轴和y轴范围

#加上图例
legend("bottomright", legend=c("1-year", "3-year", "5-year"), 
       col=c("#00468b", "#ed0000", "#42b540"), 
       lwd=2)
#调整对角线
abline(0,1,lty=5,lwd=2,col="grey")


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
c_index
