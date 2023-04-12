rm(list = ls())
setwd("/data/nas1/luchunlin/project/HF-0106-2/")
if (! dir.exists("./11_nomogram")){
  dir.create("./11_nomogram")
}
setwd("./11_nomogram")
library(magrittr)
library(ggplot2)
library(car)
# BiocManager::install('rms')
library(rms)
library(lance)
library(pROC)
# hubgene <- read.delim2('../07_features/features.xls')
hubgene <- c('MAP2K1','TAF1','HP')
dat<-read.delim2("../00_rawdata/dat(GSE169568).xls")#%>% lc.tableToNum
gene <- rownames(dat)
dat <-  as.data.frame(lapply(dat,as.numeric))
rownames(dat) <- gene
group = read.delim2("../00_rawdata/group(GSE169568).xls")
hub_exp <- t(dat[hubgene,])%>%as.data.frame()
hub_exp$group <- group$group
##基于多因素Logistic回归分析结果构建列线图
ddist <- datadist(hub_exp)
options(datadist='ddist')
hubgene
lrm <-lrm(group ~ MAP2K1+TAF1+HP, data=hub_exp, x=TRUE, y=TRUE,maxit=1000)
print(lrm)
#nomogram计算部分，此处的f_lrm及对应的多因素logistic回归函数。
pdf("01.nomogram.pdf",width=8,height=6)
nomogram <- nomogram(lrm,fun=function(x)1/(1+exp(-x)), ##逻辑回归计算公式
                     fun.at = c(0.01,0.5,0.99),#风险轴刻度
                     funlabel = "Prob of disease", #风险轴便签
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
                     funlabel = "Prob of disease", #风险轴便签
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
     xlab="Nomogram-Predicted Probability of disease risk", 
     ylab="Actual disease (proportion)", 
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
     xlab="Nomogram-Predicted Probability of disease risk", 
     ylab="Actual disease (proportion)", 
     col=c("#00468BFF", "#ED0000FF", "#42B540FF"),
     legend=FALSE)
lines(cal1[, c(1:3)], type ="l", lwd=2, pch=16, col=c("#00468BFF"))
abline(0, 1, lty=3, lwd=2) 
legend(x=.6, y=.4, legend=c("Apparent", "Bias-corrected", "Ideal"), 
       lty=c(1, 1, 2), lwd = 2, col=c("#00468BFF", "black", "black"), bty="n")
#text(x = 0.2, y = 0.8, paste0("Hosmer-Lemeshow "))
#text(x = 0.34, y = 0.8, as.expression(bquote(italic('p')==.(pval))))
dev.off()

