rm(list = ls())
setwd("/data/nas1/luchunlin/project/JNZK-255-2/")
if (! dir.exists("./15_progvalue")){
  dir.create("./15_progvalue")
}
setwd("./15_progvalue")
library(tidyverse)
survival<-read.delim2("../04_univariate_cox/survival.xls")
phenotype<-read_tsv(file = '/data/nas1/luchunlin/TCGA_phenotype/TCGA-BLCA.GDC_phenotype.tsv.gz')
train_phenotype<-data.frame(sample=phenotype$submitter_id.samples,
                            Age=phenotype$age_at_initial_pathologic_diagnosis,
                            Gender=phenotype$gender.demographic,
                            Grade=phenotype$neoplasm_histologic_grade,
                            Stage=phenotype$tumor_stage.diagnoses,
                            T.stage=phenotype$pathologic_T,
                            N.stage=phenotype$pathologic_N,
                            M.stage=phenotype$pathologic_M)

train_phenotype<-merge(train_phenotype,survival,by='sample')
# write.table(train_phenotype,file = 'phenotype.xls',sep = '\t',row.names = F,quote = F)

train_phenotype$OS<-as.numeric(train_phenotype$OS)
train_phenotype$OS.time<-as.numeric(train_phenotype$OS.time)
train_phenotype2<-train_phenotype

table(train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('a','',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('b','',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('TX',NA,train_phenotype2$T.stage)
train_phenotype2[train_phenotype2==''] <- NA
train_phenotype2$T.stage<-gsub('T0','0/1/2',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('T1','0/1/2',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('T2','0/1/2',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('T3','3',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('T4','4',train_phenotype2$T.stage)


table(train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('NX',NA,train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('N2','2/3',train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('N3','2/3',train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('N0','0',train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('N1','1',train_phenotype2$N.stage)

# train_phenotype2$N.stage[train_phenotype2$N.stage==''] <- NA


table(train_phenotype2$M.stage)
train_phenotype2$M.stage<-gsub('MX',NA,train_phenotype2$M.stage)
train_phenotype2$M.stage<-gsub('M0','0',train_phenotype2$M.stage)
train_phenotype2$M.stage<-gsub('M1','1',train_phenotype2$M.stage)

# table(train_phenotype2$Stage)
# train_phenotype2$Stage<-gsub('not reported',NA,train_phenotype2$Stage)
# # train_phenotype2$Stage<-gsub('a','',train_phenotype2$Stage)
# # train_phenotype2$Stage<-gsub('b','',train_phenotype2$Stage)
# # train_phenotype2$Stage<-gsub('c','',train_phenotype2$Stage)
# train_phenotype2$Stage<-gsub('stage iv','4',train_phenotype2$Stage)
# train_phenotype2$Stage<-gsub('stage iii','3',train_phenotype2$Stage)
# train_phenotype2$Stage<-gsub('stage ii','1/2',train_phenotype2$Stage)
# train_phenotype2$Stage<-gsub('stage i','1/2',train_phenotype2$Stage)

table(train_phenotype2$Grade)
train_phenotype2$Grade<-ifelse(train_phenotype2$Grade=='High Grade',1,0)

table(train_phenotype2$Age)
median(train_phenotype2$Age)
train_phenotype2$Age <- ifelse(train_phenotype2$Age>70,1,0)

table(train_phenotype2$Gender)
train_phenotype2$Gender <- ifelse(train_phenotype2$Gender=='female',1,0)

colnames(train_phenotype2)
colnames(train_phenotype2)[1]<-'id'
risk<-read.delim2('../06_risk/risk.xls')
sub_risk <- subset(risk, select = c(id, riskScore))

train_risk_clinical <- merge(train_phenotype2,
                             sub_risk,
                             by = "id")
rownames(train_risk_clinical) <- train_risk_clinical$id
train_risk_clinical = subset(train_risk_clinical, select = -c(id))
dim(train_risk_clinical)
colnames_train <- colnames(train_risk_clinical)
covariates_train <- colnames_train[-which(colnames_train %in% c("OS", "OS.time"))]
class(train_risk_clinical$riskScore)
train_risk_clinical$riskScore <- as.numeric(train_risk_clinical$riskScore)

train_risk_clinical$T.stage<-factor(train_risk_clinical$T.stage)
train_risk_clinical$N.stage<-factor(train_risk_clinical$N.stage)
train_risk_clinical$M.stage<-factor(train_risk_clinical$M.stage)
train_risk_clinical$Grade <- factor(train_risk_clinical$Grade)
# train_risk_clinical$Stage <- factor(train_risk_clinical$Stage)

library(survival)
res.risk = coxph(Surv(time = OS.time, event = OS) ~ riskScore, data = train_risk_clinical) %>% summary
res.risk = c(res.risk$conf.int[-2], res.risk$coefficients[5])
res.age = coxph(Surv(time = OS.time, event = OS) ~ Age, data = train_risk_clinical) %>% summary
res.age = c(res.age$conf.int[-2], res.age$coefficients[5])

res.gender = coxph(Surv(time = OS.time, event = OS) ~ Gender, data = train_risk_clinical) %>% summary
res.gender = c(res.gender$conf.int[-2], res.gender$coefficients[5])
# 
# res.stage = coxph(Surv(time = OS.time, event = OS) ~ Stage, data = train_risk_clinical) %>% summary
# res.stage = cbind(res.stage$conf.int[,-2], res.stage$coefficients[,5])

res.T_stage = coxph(Surv(time = OS.time, event = OS) ~ T.stage, data = train_risk_clinical) %>% summary
res.T_stage = cbind(res.T_stage$conf.int[,-2], res.T_stage$coefficients[,5])

res.M_stage = coxph(Surv(time = OS.time, event = OS) ~ M.stage, data = train_risk_clinical) %>% summary
res.M_stage = c(res.M_stage$conf.int[-2], res.M_stage$coefficients[5])

res.N_stage = coxph(Surv(time = OS.time, event = OS) ~ N.stage, data = train_risk_clinical) %>% summary
res.N_stage = cbind(res.N_stage$conf.int[,-2], res.N_stage$coefficients[,5])

res.grade = coxph(Surv(time = OS.time, event = OS) ~ Grade, data = train_risk_clinical) %>% summary
res.grade = c(res.grade$conf.int[-2], res.grade$coefficients[5])

res.ref = c(1,1,1,NA)

res = rbind(res.risk,res.age,res.gender,res.ref,res.T_stage,res.ref,res.N_stage,res.M_stage,res.grade) %>% as.data.frame()
rownames(res)
res$Indicators = c("riskScore","Age",'Gender',"T0/1/2(Reference)",'T3','T4',"N0(Reference)",'N1','N2/3','M1 vs.M0','High Grade vs.Low Grade')
colnames(res) = c("hr","low","up","pv","Indicator")
res$p = signif(res$pv, 2) %>% paste0("p = ", .)
res$p[is.na(res$pv)] = NA
res$Indicator = factor(res$Indicator, levels = rev(res$Indicator))
rownames(res) <- res$Indicator
res2 <- data.frame(p.value=res$pv,
                   HR=res$hr,
                   HR.95L=res$low,
                   HR.95H=res$up,
                   Indicator=res$Indicator)
rownames(res2) <- res2$Indicator
write.table(res2, file = "univariate_cox_prog_forest.xls", sep = "\t", quote = F)
res2 <- subset(res2, select = -c(Indicator))
library(tidyr)
hz <- paste(round(res2$HR,3),
            "(",round(res2$HR.95L,3),
            "-",round(res2$HR.95H,3),")",sep = "")
hz
hz[c(4,7)] <- ""

tabletext <- cbind(c(NA,rownames(res2)),
                   c("P value",ifelse(res2$p.value<0.0001,
                                      "< 0.0001",
                                      round(res2$p.value,4))),
                   c("Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1
library(forestplot)
pdf(file = '01.uncoxforest.pdf',height = 6, width = 13, onefile = F)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,rep(FALSE, 7)),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,res2$HR),
           lower=c(NA,res2$HR.95L), #95%置信区间下限
           upper=c(NA,res2$HR.95H), #95%置信区间上限
           boxsize=0.1,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(1.2,"cm"), #固定行高
           graphwidth = unit(.6,"npc"), #图在表中的宽度比例
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),
           # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
           #                 "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
           #                 "59" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
           # mar=unit(rep(0.01, times = 4), "cm"),#图形页边距
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Hazard Ratio",
           grid = T,
           title = "Univariate") # 垂直于x轴的网格线，对应每个刻度
dev.off()

png(file = '01.uncoxforest.png',height = 500, width = 1000)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,rep(FALSE, 7)),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,res2$HR),
           lower=c(NA,res2$HR.95L), #95%置信区间下限
           upper=c(NA,res2$HR.95H), #95%置信区间上限
           boxsize=0.1,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(0,1,2,3,4,5,6,7,8,9,10,11,12,13), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(1.2,"cm"), #固定行高
           graphwidth = unit(.6,"npc"), #图在表中的宽度比例
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),
           # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
           #                 "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
           #                 "59" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
           # mar=unit(rep(0.01, times = 4), "cm"),#图形页边距
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Hazard Ratio",
           grid = T,
           title = "Univariate") # 垂直于x轴的网格线，对应每个刻度
dev.off()
## 07-2 多因素Cox----------
res.mul = coxph(Surv(time = OS.time, event = OS) ~ riskScore +Age+T.stage+N.stage+M.stage, data = train_risk_clinical)%>% summary
res.mul = cbind(res.mul$conf.int[,-2], res.mul$coefficients[,5]) %>% as.data.frame()
rownames(res.mul)

res.mul = rbind(res.mul[c(1:2),], res.ref, res.mul[c(3:4),],res.ref,res.mul[c(5:7),])
res.mul$Indicators = c("riskScore","Age","T0/1/2(Reference)",'T3','T4',"N0(Reference)",'N1','N2/3','M1 vs.M0')
colnames(res.mul) = c("hr","low","up","pv","Indicator")
res.mul$p = signif(res.mul$pv, 2) %>% paste0("p = ", .)
res.mul$p[is.na(res.mul$pv)] = NA
res.mul$Indicator = factor(res.mul$Indicator, levels = rev(res.mul$Indicator))
rownames(res.mul) <- res.mul$Indicator

multi_res <- data.frame(p.value=res.mul$pv,
                        HR=res.mul$hr,
                        HR.95L=res.mul$low,
                        HR.95H=res.mul$up,
                        Indicator=res.mul$Indicator)
rownames(multi_res) <- multi_res$Indicator
multi_res
write.csv(multi_res,
          file = "multivariate_cox_prog_result.csv",
          quote = F,
          row.names = T)
multi_res <- subset(multi_res, select = -c(Indicator))

library(tidyr)
hz <- paste(round(multi_res$HR,3),
            "(",round(multi_res$HR.95L,3),
            "-",round(multi_res$HR.95H,3),")",sep = "")
hz

hz[c(3,6)] <- ""
tabletext <- cbind(c(NA,rownames(multi_res)),
                   c("P value",ifelse(multi_res$p.value<0.0001,
                                      "< 0.0001",
                                      round(multi_res$p.value,4))),
                   c("Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1
pdf(file = '02.mulcoxforest.pdf',height = 5, width = 12, onefile = F)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,multi_res$HR),
           lower=c(NA,multi_res$HR.95L), #95%置信区间下限
           upper=c(NA,multi_res$HR.95H), #95%置信区间上限
           boxsize=0.1,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(0,1,2,3,4,5,6), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(1.2,"cm"), #固定行高
           graphwidth = unit(.6,"npc"), #图在表中的宽度比例
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),
           # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
           #                 "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
           #                 "59" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
           # mar=unit(rep(0.01, times = 4), "cm"),#图形页边距
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Hazard Ratio",
           grid = T,
           title = "Multivariate") # 垂直于x轴的网格线，对应每个刻度
dev.off()

png(file = '02.mulcoxforest.png',height = 400, width = 800)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE,FALSE),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,multi_res$HR),
           lower=c(NA,multi_res$HR.95L), #95%置信区间下限
           upper=c(NA,multi_res$HR.95H), #95%置信区间上限
           boxsize=0.1,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(0,1,2,3,4,5,6), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(1.2,"cm"), #固定行高
           graphwidth = unit(.6,"npc"), #图在表中的宽度比例
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),
           # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
           #                 "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
           #                 "59" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
           # mar=unit(rep(0.01, times = 4), "cm"),#图形页边距
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Hazard Ratio",
           grid = T,
           title = "Multivariate") # 垂直于x轴的网格线，对应每个刻度
dev.off()
##07-3 构建COX模型，绘制列线图---------
multi_cov<-c('riskScore',"N.stage")

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
#train_risk_clinical$Age <- ifelse(train_risk_clinical$Age==1,'>60','<=60')
res.cox <- psm(cox_data_prog,
               data = train_risk_clinical, dist = 'lognormal')
surv <- Survival(res.cox) # 构建生存概率函数
function(x) surv(365, x) # 1年事件发生概率
#function(x) surv(730, x) # 2年事件发生概率
function(x) surv(1095, x) # 3年事件发生概率
function(x) surv(1825, x) # 5年事件发生概率
# function(x) surv(2555, x) # 7年事件发生概率
nom.cox <- nomogram(res.cox,
                    fun = list(function(x) surv(365, x),
                               function(x) surv(1095, x),
                               function(x) surv(1825, x)),
                    funlabel=c("1-year Survival Probability", "3-year Survival Probability", "5-year Survival Probability"),
                    maxscale = 10,
                    fun.at = c(0.01,seq(0.1,0.9,by=0.2),0.95,0.99),
                    lp=F)
plot(nom.cox, cex.axis  = 1.5, cex.var = 1.7)
png(filename = "03.nomogram_line_points.png", height = 600, width = 1200)
plot(nom.cox, cex.axis  = 1.5, cex.var = 1.7)
dev.off()
pdf(file = "03.nomogram_line_points.pdf", height = 8, width = 17)
plot(nom.cox, cex.axis  = 1.5, cex.var = 1.7)
dev.off()
##07-4 构建校准曲线---------
coxm_1 <- cph(cox_data_prog,
              data=train_risk_clinical,
              surv=T,
              x=T,
              y=T,
              time.inc = 365*1)
cal_1<-calibrate(coxm_1,u=365*1,cmethod='KM',m=100,B=200)

##绘制3年生存期校曲线
##time.in 和 u 要是一样的，都是要评价的时间节点
coxm_3 <- cph(cox_data_prog,
              data=train_risk_clinical,
              surv=T,
              x=T,y=T,
              time.inc = 3*365)
cal_3 <-calibrate(coxm_3,u=3*365,cmethod='KM',m=100,B=200)

coxm_5 <- cph(cox_data_prog,
              data=train_risk_clinical,
              surv=T,
              x=T,y=T,
              time.inc = 5*365)
cal_5 <-calibrate(coxm_5,u=5*365,cmethod='KM',m=100,B=1000)
pdf(file = '04.calibrate.pdf',w=9,h=8)
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
     ylab='Actual 1-5 year Survival Probability',#标签
     col="#42b540",#设置一个颜色
     xlim = c(0,1),ylim = c(0,1)) ##x轴和y轴范围

#加上图例
legend("bottomright", legend=c("1-year", "3-year", "5-year"), 
       col=c("#00468b", "#ed0000", "#42b540"), 
       lwd=2)
#调整对角线
abline(0,1,lty=5,lwd=2,col="grey")
dev.off()

png(file = '04.calibrate.png',w=700,h=600)
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
     ylab='Actual 1-5 year Survival Probability',#标签
     col="#42b540",#设置一个颜色
     xlim = c(0,1),ylim = c(0,1)) ##x轴和y轴范围

#加上图例
legend("bottomright", legend=c("1-year", "3-year", "5-year"), 
       col=c("#00468b", "#ed0000", "#42b540"), 
       lwd=2)
#调整对角线
abline(0,1,lty=5,lwd=2,col="grey")
dev.off()

##DCA曲线-------
# ## DCA-----
library(rmda)
dca.dat <- train_risk_clinical
# dca.dat$Age <- ifelse(dca.dat$Age=='>60',1,0)
complex<-decision_curve(OS ~riskScore+N.stage,data =dca.dat,family = binomial(link ='logit'),
                        thresholds = seq(0,1, by = 0.01),
                        confidence.intervals= 0.95,
                        study.design = 'case-control',
                        population.prevalence= 0.3
)
riskScore.dca <- decision_curve(OS ~riskScore,data =dca.dat,family = binomial(link ='logit'),
                                thresholds = seq(0,1, by = 0.01),
                                confidence.intervals= 0.95,
                                study.design = 'case-control',
                                population.prevalence= 0.3
)

nstage.dca <- decision_curve(OS ~N.stage,data =dca.dat,family = binomial(link ='logit'),
                            thresholds = seq(0,1, by = 0.01),
                            confidence.intervals= 0.95,
                            study.design = 'case-control',
                            population.prevalence= 0.3
)

dca.list <- list(riskScore.dca,nstage.dca,complex)


pdf("05.DCA.pdf", width = 5, height = 5)
plot_decision_curve(dca.list,
                    curve.names=c('riskScore','N.Stage','Nomogram'),
                    cost.benefit.axis =FALSE,col= c("#FDC086","#7FC97F","#BEAED4","#74CCBE","#ED7474"),
                    confidence.intervals=FALSE,
                    standardize = FALSE)
dev.off()
png("05.DCA.png", width = 500, height = 500)
plot_decision_curve(dca.list,
                    curve.names=c('riskScore','N.Stage','Nomogram'),
                    cost.benefit.axis =FALSE,col= c("#FDC086","#7FC97F","#BEAED4","#74CCBE","#ED7474"),
                    confidence.intervals=FALSE,
                    standardize = FALSE)
dev.off()




