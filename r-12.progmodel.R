rm(list = ls())
# 09 预后模型的构建与评价--------
setwd("/data/nas1/luchunlin/project/BJTC-308/")
if (! dir.exists("./12_prog_model")){
  dir.create("./12_prog_model")
}
setwd("./12_prog_model")

## 07-1 单因素Cox----------
train_phenotype<-read.delim2("/data/nas1/luchunlin/project/BJTC-308/07_clinical_index/phenotype.xls")
train_phenotype$OS<-as.numeric(train_phenotype$OS)
train_phenotype$OS.time<-as.numeric(train_phenotype$OS.time)
train_phenotype2<-train_phenotype
table(train_phenotype2$TNM.stage)
train_phenotype2$TNM.stage<-gsub('tnm stage: III','3',train_phenotype2$TNM.stage,fixed = T)
train_phenotype2$TNM.stage<-gsub('tnm stage: II','1/2',train_phenotype2$TNM.stage,fixed = T)
train_phenotype2$TNM.stage<-gsub('tnm stage: I','1/2',train_phenotype2$TNM.stage,fixed = T)

table(train_phenotype2$T.stage)

train_phenotype2$T.stage<-gsub('t stage: T1','1',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('t stage: T2','2',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('t stage: T3','3/4',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('t stage: T4','4/4',train_phenotype2$T.stage)

table(train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('n stage: N0','0',train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('n stage: N1','1',train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('n stage: N2','2/3',train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('n stage: N3','2/3',train_phenotype2$N.stage)

train_phenotype2$gender<-ifelse(train_phenotype2$gender=='female',0,1)

colnames(train_phenotype2)
colnames(train_phenotype2)<-c('id','Gender','T.stage','N.stage','TNM.stage','OS','OS.time')
risk<-read.delim2('/data/nas1/luchunlin/project/BJTC-308/05_risk/risk.xls')%>%lc.tableToNum()
sub_risk <- subset(risk, select = c(id, riskScore))

train_risk_clinical <- merge(train_phenotype2,
                             sub_risk,
                             by = "id")
rownames(train_risk_clinical) <- train_risk_clinical$id
train_risk_clinical = subset(train_risk_clinical, select = -c(id))
dim(train_risk_clinical)
colnames_train <- colnames(train_risk_clinical)
covariates_train <- colnames_train[-which(colnames_train %in% c("OS", "OS.time"))]

train_risk_clinical$TNM.stage<-factor(train_risk_clinical$TNM.stage)
train_risk_clinical$T.stage<-factor(train_risk_clinical$T.stage)
train_risk_clinical$N.stage<-factor(train_risk_clinical$N.stage)
train_risk_clinical$Gender<-factor(train_risk_clinical$Gender)

library(survival)
res.risk = coxph(Surv(time = OS.time, event = OS) ~ riskScore, data = train_risk_clinical) %>% summary
res.risk = c(res.risk$conf.int[-2], res.risk$coefficients[5])
res.gender = coxph(Surv(time = OS.time, event = OS) ~ Gender, data = train_risk_clinical) %>% summary
res.gender = c(res.gender$conf.int[-2], res.gender$coefficients[5])

res.T_stage = coxph(Surv(time = OS.time, event = OS) ~ T.stage, data = train_risk_clinical) %>% summary
res.T_stage = cbind(res.T_stage$conf.int[,-2], res.T_stage$coefficients[,5])

res.N_stage = coxph(Surv(time = OS.time, event = OS) ~ N.stage, data = train_risk_clinical) %>% summary
res.N_stage = cbind(res.N_stage$conf.int[,-2], res.N_stage$coefficients[,5])

res.stage = coxph(Surv(time = OS.time, event = OS) ~ TNM.stage, data = train_risk_clinical) %>% summary

res.stage = c(res.stage$conf.int[-2], res.stage$coefficients[5])
res.ref = c(1,1,1,NA)

res = rbind(res.risk,res.gender,res.ref,res.T_stage,res.ref,res.N_stage,res.stage) %>% as.data.frame()
rownames(res)
res$Indicators = c("riskScore","Gender","T1(Reference)",'T2','T3',"T4",'N1(Reference)','N1','N2/3','stage3 vs.stage1/2')
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
hz[c(3,7)] <- ""

tabletext <- cbind(c(NA,rownames(res2)),
                   c("P value",ifelse(res2$p.value<0.0001,
                                      "< 0.0001",
                                      round(res2$p.value,4))),
                   c("Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1
library(forestplot)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,rep(FALSE, 7)),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,res2$HR),
           lower=c(NA,res2$HR.95L), #95%置信区间下限
           upper=c(NA,res2$HR.95H), #95%置信区间上限
           boxsize=0.1,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(0,5,10,15,20,25,30,35,40,45,50,55,60), #横坐标刻度
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
## 07-2 多因素Cox----------
res.mul = coxph(Surv(time = OS.time, event = OS) ~ riskScore +TNM.stage, data = train_risk_clinical)%>% summary
res.mul = cbind(res.mul$conf.int[,-2], res.mul$coefficients[,5]) %>% as.data.frame()
rownames(res.mul)
res.mul = rbind(res.mul[1,],res.mul[2,])
res.mul$Indicators = c("riskScore","stage3 vs.stage1/2")
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

#hz[c(2,5)] <- ""
tabletext <- cbind(c(NA,rownames(multi_res)),
                   c("P value",ifelse(multi_res$p.value<0.0001,
                                      "< 0.0001",
                                      round(multi_res$p.value,4))),
                   c("Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1
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
           xticks = c(0,1,2,3,4,5,6,7,8,9,10), #横坐标刻度
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
##07-3 构建COX模型，绘制列线图---------
multi_cov<-c('riskScore',"TNM.stage")
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
function(x) surv(730, x) # 2年事件发生概率
function(x) surv(1095, x) # 3年事件发生概率

nom.cox <- nomogram(res.cox,
                    fun = list(function(x) surv(365, x),
                               function(x) surv(730, x),
                               function(x) surv(1095, x)),
                    funlabel=c("1-year Survival Probability", "2-year Survival Probability", "3-year Survival Probability"),
                    maxscale = 10,
                    fun.at = c(0.01,seq(0.1,0.9,by=0.2),0.95,0.99),
                    lp=F)



plot(nom.cox, cex.axis  = 1.5, cex.var = 1.7)
png(filename = "nomogram_line_points.png", height = 600, width = 1200)
plot(nom.cox, cex.axis  = 1.5, cex.var = 1.7)
dev.off()
pdf(file = "nomogram_line_points.pdf", height = 8, width = 17)
plot(nom.cox, cex.axis  = 1.5, cex.var = 1.7)
dev.off()
##07-4 构建校准曲线---------
coxm_1 <- cph(cox_data_prog,
              data=train_risk_clinical,
              surv=T,
              x=T,
              y=T,
              time.inc = 365)
cal_1<-calibrate(coxm_1,u=365,cmethod='KM',m=30)

##绘制3年生存期校曲线
##time.in 和 u 要是一样的，都是要评价的时间节点
coxm_3 <- cph(cox_data_prog,
              data=train_risk_clinical,
              surv=T,
              x=T,y=T,
              time.inc = 3*365)
cal_3 <-calibrate(coxm_3,u=3*365,cmethod='KM',m=30)

coxm_2 <- cph(cox_data_prog,
              data=train_risk_clinical,
              surv=T,
              x=T,y=T,
              time.inc = 2*365)
cal_2 <-calibrate(coxm_2,u=2*365,cmethod='KM',m=30)

par(mar=c(7,4,4,3),cex=1.5)
plot(cal_1,
     subtitles = F,
     lwd=2,lty=1, ##设置线条形状和尺寸
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 1-3 year Progression-free Interval',#便签
     ylab='Actual 1-3 year Progression-free Interval(Proportion)',#标签
     col="#00468b",#设置一个颜色
     xlim = c(0,1),ylim = c(0,1)) ##x轴和y轴范围
plot(cal_2,
     add = T,
     subtitles = F,
     lwd=2,lty=1,  ##设置线条宽度和线条类型
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 1-3 year Progression-free Interval',#便签
     ylab='Actual 1-3 year Progression-free Interval(Proportion)',#标签
     col="#ed0000",#设置一个颜色
     xlim = c(0,1),ylim = c(0,1)) ##x轴和y轴范围
plot(cal_3,
     add = T,
     subtitles = F,
     lwd=2,lty=1, ##设置线条形状和尺寸
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 1-3 year Progression-free Interval',#便签
     ylab='Actual 1-3year Progression-free Interval(Proportion)',#标签
     col="#42b540",#设置一个颜色
     xlim = c(0,1),ylim = c(0,1)) ##x轴和y轴范围

#加上图例
legend("bottomright", legend=c("1-year", "2-year", "3-year"), 
       col=c("#00468b", "#ed0000", "#42b540"), 
       lwd=2)
#调整对角线
abline(0,1,lty=5,lwd=2,col="grey")
# ## 07-5 KM----------
# res.mul = coxph(Surv(time = OS.time, event = OS)~ riskScore +TNM.stage, data = train_risk_clinical)
# train_risk_clinical$riskScore_prog = predict(res.mul, newdata = train_risk_clinical, type = "lp")
# 
# train_risk_clinical$risk_prog = ifelse(train_risk_clinical$riskScore_prog > median(train_risk_clinical$riskScore_prog, na.rm = T), "high", "low")
# train_risk_clinical$risk_prog = factor(train_risk_clinical$risk_prog, levels = c("high", "low"), labels = c("High risk", "Low risk"))
# surv.fit = survfit(Surv(time = OS.time,event = OS) ~ risk_prog, data = train_risk_clinical)
# prog_survival_median <- ggsurvplot(surv.fit,
#                                    pval = TRUE, 
#                                    conf.int = F,
#                                    legend.labs=c("High risk","Low risk" ),
#                                    legend.title="Risk score",
#                                    title="Prognosis Model KM",
#                                    font.main = c(15,"bold"),
#                                    risk.table = TRUE, 
#                                    risk.table.col = "strata", 
#                                    linetype = "strata", 
#                                    surv.median.line = "hv",
#                                    ggtheme = theme_bw(), 
#                                    palette = c("#A73030FF", "#0073C2FF"))
# prog_survival_median
# 
# ## 07-6 ROC--------
# library(survivalROC)
# library(tidyverse)
# # 开始验证
# train_risk_clinical2 <- train_risk_clinical
# 
# library(survival)
# library(survminer)
# 
# riskscore <- function(survival_cancer_df,
#                       candidate_genes_for_cox,
#                       cox_report){
#   library("dplyr")
#   risk_score_table <- survival_cancer_df[, candidate_genes_for_cox]
#   for (each_sig_gene in colnames(risk_score_table)){
#     risk_score_table$each_sig_gene <- risk_score_table[,each_sig_gene]*
#       (summary(cox_report)$coefficients[each_sig_gene,1])
#   }
#   risk_score_table <- cbind(risk_score_table,
#                             "total_risk_score"=exp(rowSums(risk_score_table))) %>%
#     cbind(survival_cancer_df[,c("OS.time", "OS")])
#   risk_score_table <- risk_score_table[,c("OS.time",
#                                           "OS",
#                                           candidate_genes_for_cox,
#                                           "total_risk_score")]
# }
# 
# 
# # candidate_genes_for_cox2 <- c(rownames(cox_table)[cox_table[,3]>0.05])
# 
# risk_score_table_multi_cox2<-train_risk_clinical2[,-10]
# multi_ROC_out <- function(time_vector, risk_score_table){
#   library(survivalROC)
#   single_ROC <- function(single_time){
#     for_ROC <- survivalROC(Stime=risk_score_table$OS.time,
#                            status=risk_score_table$OS,
#                            marker=risk_score_table$riskScore_prog,
#                            predict.time=single_time,method = 'KM')
#     data.frame('True_positive'=for_ROC$TP, 'False_positive'=for_ROC$FP,
#                'Cut_values'=for_ROC$cut.values, 'Time'=rep(single_time, length(for_ROC$TP)),
#                'AUC'=rep(for_ROC$AUC, length(for_ROC$TP)))
#   }
#   multi_ROC_list <- lapply(time_vector, single_ROC)
#   do.call(rbind, multi_ROC_list)
# }
# 
# for_multi_ROC <- multi_ROC_out(time_vector = c(365*seq(1,5,2)), 
#                                risk_score_table = train_risk_clinical2)
# for_multi_ROC$Time <- factor(for_multi_ROC$Time)
# 
# # 画ROC曲线 
# library(scales)
# library(geomROC)
# library(plotROC)
# auc_y1 <- round(for_multi_ROC[which(for_multi_ROC$Time==365),5][1],2)
# auc_y3 <- round(for_multi_ROC[which(for_multi_ROC$Time==1095),5][1],2)
# auc_y5 <- round(for_multi_ROC[which(for_multi_ROC$Time==1825),5][1],2)
# 
# ROC <- ggplot(for_multi_ROC, aes(x=False_positive,
#                                  y=True_positive, 
#                                  label=Cut_values, 
#                                  color=Time)) + 
#   scale_color_manual(breaks = c("365", "1095", "1825"),
#                      labels = c("1 years", "3 years", "5 years"),
#                      values = c("#4682B4", "#FF4040", "#20B2AA")) +
#   geom_roc(labels = F, stat = 'identity') + 
#   style_roc() + 
#   geom_abline(slope = 1, intercept = 0, color = 'gray', linetype=2) +
#   theme_bw() +
#   labs(title = "ROC for Prognosis Model") +
#   theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
#         panel.grid = element_blank(),
#         plot.title = element_text(size = 15, hjust = 0.5, face = "bold")) +
#   # annotate('text', x=.75, y=.25, label=paste('AUC of 1 years =', round(auc_y1,2))) + 
#   # annotate('text', x=.75, y=.15, label=paste('AUC of 2 years =', round(auc_y3,2))) + 
#   # annotate('text', x=.75, y=.05, label=paste('AUC of 3 years =', round(auc_y5,2))) +
#   annotate("text", x=0.75, y=c(0.25, 0.15, 0.05),
#            label = c(paste('AUC of 1 years =', format(auc_y1,nsmall=2)),
#                      paste('AUC of 3 years =', format(auc_y3,nsmall=2)),
#                      paste('AUC of 5 years =', format(auc_y5,nsmall=2))))
# ROC
# ggsave('ROC for Prognosis Model.png', ROC,width = 5, height = 4)
# ggsave('ROC for Prognosis Model.pdf', ROC,width = 5, height = 4)
