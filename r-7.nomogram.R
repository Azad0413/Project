rm(list = ls())
setwd("/data/nas1/luchunlin/project/TY0307-11/")
if (! dir.exists("./07_nomogram")){
  dir.create("./07_nomogram")
}
setwd("./07_nomogram")
library(magrittr)
library(ggplot2)
library(car)
library(rms)
library(lance)
library(pROC)
hubgene <- read.delim2('../04_model/hubgene.xls')
dat<-read.delim2("../00_rawdata/dat(GSE113079).xls", row.names = 1)%>% lc.tableToNum
group = read.delim2("../00_rawdata/group(GSE113079).xls")
hub_exp <- t(dat[hubgene$x,])%>%as.data.frame()
hub_exp$group <- group$group
##基于多因素Logistic回归分析结果构建列线图
ddist <- datadist(hub_exp)
options(datadist='ddist')
hubgene
lrm <-lrm(group ~ RCN2+HRC+DERL2+RNF183+CRH+TMED2+PPP1R15A+IL1A, data=hub_exp, x=TRUE, y=TRUE,maxit=1000)
print(lrm)
#nomogram计算部分，此处的f_lrm及对应的多因素logistic回归函数。
pdf("01.nomogram.pdf",width=8,height=6)
nomogram <- nomogram(lrm,fun=function(x)1/(1+exp(-x)), ##逻辑回归计算公式
                     fun.at = c(0.01,0.5,0.99),#风险轴刻度
                     funlabel = "Prob of death", #风险轴便签
                     lp=F,  ##是否显示系数轴
                     conf.int = F, ##每个得分的置信度区间，用横线表示,横线越长置信度越
                     abbrev = F#是否用简称代表因子变量
)
#绘制nomogram
plot(nomogram)
dev.off()
png("01.nomogram.png",width=600,height=400)
nomogram <- nomogram(lrm,fun=function(x)1/(1+exp(-x)), ##逻辑回归计算公式
                     fun.at = c(0.01,0.5,0.99),#风险轴刻度
                     funlabel = "Prob of death", #风险轴便签
                     lp=F,  ##是否显示系数轴
                     conf.int = F, ##每个得分的置信度区间，用横线表示,横线越长置信度越
                     abbrev = F#是否用简称代表因子变量
)
#绘制nomogram
plot(nomogram)
dev.off()

cal1 <- calibrate(lrm, cmethod='KM', method='boot', B=30)#建模组中绘制校准曲线
pdf("02.calibrate.pdf",width=8,height=6, family = "Times")
par(mar = c(6,5,2,2))
plot(cal1, lwd=2, lty=1, 
     cex.lab=1.5, cex.axis=1.2, cex.main=1.5, cex.sub=1.2, 
     xlim=c(0, 1), ylim= c(0, 1), 
     xlab="Nomogram-Predicted Probability of death risk", 
     ylab="Actual death (proportion)", 
     col=c("#00468BFF", "#ED0000FF", "#42B540FF"),
     legend=FALSE)
lines(cal1[, c(1:3)], type ="l", lwd=2, pch=16, col=c("#00468BFF"))
abline(0, 1, lty=3, lwd=2) 
legend(x=.6, y=.4, legend=c("Apparent", "Bias-corrected", "Ideal"), 
       lty=c(1, 1, 2), lwd = 2, col=c("#00468BFF", "black", "black"), bty="n")
#text(x = 0.2, y = 0.8, paste0("Hosmer-Lemeshow "))
#text(x = 0.34, y = 0.8, as.expression(bquote(italic('p')==.(pval))))
dev.off()
png("02.calibrate.png",width=8,height=6, units = "in", res = 600, family = "Times")
par(mar = c(6,5,2,2))
plot(cal1, lwd=2, lty=1, 
     cex.lab=1.5, cex.axis=1.2, cex.main=1.5, cex.sub=1.2, 
     xlim=c(0, 1), ylim= c(0, 1), 
     xlab="Nomogram-Predicted Probability of death risk", 
     ylab="Actual death (proportion)", 
     col=c("#00468BFF", "#ED0000FF", "#42B540FF"),
     legend=FALSE)
lines(cal1[, c(1:3)], type ="l", lwd=2, pch=16, col=c("#00468BFF"))
abline(0, 1, lty=3, lwd=2) 
legend(x=.6, y=.4, legend=c("Apparent", "Bias-corrected", "Ideal"), 
       lty=c(1, 1, 2), lwd = 2, col=c("#00468BFF", "black", "black"), bty="n")
#text(x = 0.2, y = 0.8, paste0("Hosmer-Lemeshow "))
#text(x = 0.34, y = 0.8, as.expression(bquote(italic('p')==.(pval))))
dev.off()


colnames(hub_exp)
library(rmda)
# pdf("03.DCA.pdf",width=8,height=8)
# png("03.DCA.png",width=800,height=800)
hub_exp$group<-ifelse(hub_exp$group=='CAD',1,0)
colnames(hub_exp)
complex<-decision_curve(group ~ RCN2+HRC+DERL2+RNF183+CRH+TMED2+PPP1R15A+IL1A,data = hub_exp,family = binomial(link ='logit'),
                        thresholds = seq(0,1, by = 0.01),
                        confidence.intervals= 0.95,
                        study.design = 'case-control',
                        population.prevalence= 0.3
)
# RCN2.dca <- decision_curve(group ~ RCN2,data = hub_exp,family = binomial(link ='logit'),
#                            thresholds = seq(0,1, by = 0.01),
#                            confidence.intervals= 0.95,
#                            study.design = 'case-control',
#                            population.prevalence= 0.3
# )
# HRC.dca <- decision_curve(group ~ HRC,data = hub_exp,family = binomial(link ='logit'),
#                thresholds = seq(0,1, by = 0.01),
#                confidence.intervals= 0.95,
#                study.design = 'case-control',
#                population.prevalence= 0.3
# )
# DERL2.dca <- decision_curve(group ~ DERL2,data = hub_exp,family = binomial(link ='logit'),
#                            thresholds = seq(0,1, by = 0.01),
#                            confidence.intervals= 0.95,
#                            study.design = 'case-control',
#                            population.prevalence= 0.3
# )
# RNF183.dca <- decision_curve(group ~ RNF183,data = hub_exp,family = binomial(link ='logit'),
#                             thresholds = seq(0,1, by = 0.01),
#                             confidence.intervals= 0.95,
#                             study.design = 'case-control',
#                             population.prevalence= 0.3
# )
# CRH.dca <- decision_curve(group ~ CRH,data = hub_exp,family = binomial(link ='logit'),
#                              thresholds = seq(0,1, by = 0.01),
#                              confidence.intervals= 0.95,
#                              study.design = 'case-control',
#                              population.prevalence= 0.3
# )
# TMED2.dca <- decision_curve(group ~ TMED2,data = hub_exp,family = binomial(link ='logit'),
#                           thresholds = seq(0,1, by = 0.01),
#                           confidence.intervals= 0.95,
#                           study.design = 'case-control',
#                           population.prevalence= 0.3
# )
# PPP1R15A.dca <- decision_curve(group ~ PPP1R15A,data = hub_exp,family = binomial(link ='logit'),
#                             thresholds = seq(0,1, by = 0.01),
#                             confidence.intervals= 0.95,
#                             study.design = 'case-control',
#                             population.prevalence= 0.3
# )
# IL1A.dca <- decision_curve(group ~ IL1A,data = hub_exp,family = binomial(link ='logit'),
#                                thresholds = seq(0,1, by = 0.01),
#                                confidence.intervals= 0.95,
#                                study.design = 'case-control',
#                                population.prevalence= 0.3
# )
# dca.list <- list(RCN2.dca,HRC.dca,DERL2.dca,RNF183.dca,CRH.dca,TMED2.dca,PPP1R15A.dca,IL1A.dca,complex)


# plot_decision_curve(dca.list,
#                     curve.names=c('RCN2','HRC','DERL2','RNF183','CRH','TMED2','PPP1R15A','IL1A','nomogram'),
#                     cost.benefit.axis =FALSE,col= c("#FDC086","#7FC97F","#BEAED4","#74CCBE","#ED7474",'orange','blue','darkgreen','red'),
#                     confidence.intervals=FALSE,
#                     standardize = FALSE)
pdf("03.DCA.pdf",width=6,height=5)
plot_decision_curve(complex,
                    curve.names=c('nomogram'),
                    cost.benefit.axis =FALSE,col= c('blue'),
                    confidence.intervals=FALSE,
                    standardize = FALSE,legend.position = 'right')
dev.off()
png("03.DCA.png",width=500,height=400)
plot_decision_curve(complex,
                    curve.names=c('nomogram'),
                    cost.benefit.axis =FALSE,col= c('blue'),
                    confidence.intervals=FALSE,
                    standardize = FALSE,legend.position = 'right')
dev.off()
##临床影响曲线CIC-------
pdf(file = '04.CIC.pdf', width = 6, height = 5)
plot_clinical_impact(complex,population.size = 1000,cost.benefit.axis = T,
                     n.cost.benefits= 8,col = c('red','blue'),
                     confidence.intervals= T,ylim=c(0,1000),
                     legend.position= "topright")
dev.off()
#红色曲线（Numberhigh risk）表示，在各个阈概率下，被simple或complex模型划分为阳性（高风险）的人数；
#蓝色曲线（Number high risk with outcome）为各个阈概率下真阳性的人数。意义一目了然。
png(file = '04.CIC.png', width = 400, height = 350)
plot_clinical_impact(complex,population.size = 1000,cost.benefit.axis = T,
                     n.cost.benefits= 8,col = c('red','blue'),
                     confidence.intervals= T,ylim=c(0,1000),
                     legend.position= "topright")
dev.off()
