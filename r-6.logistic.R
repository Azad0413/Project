rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-302/")
if (! dir.exists("./06_logistic")){
  dir.create("./06_logistic")
}
setwd("./06_logistic")
library(magrittr)
library(ggplot2)
library(car)
library(rms)
library(lance)
library(pROC)
hub_gene <- read.delim2('/data/nas1/luchunlin/project/BJTC-302/10_validation/hubgene.final.xls')
train.dat<-read.delim2("/data/nas1/luchunlin/project/BJTC-302/00_rawdata/dat.xls", row.names = 1)  %>% lc.tableToNum
sample<-read.delim2('/data/nas1/luchunlin/project/BJTC-302/00_rawdata/group.xls')
sample<-sample[which(sample$group=='Burn'),]
group = read.delim2("/data/nas1/luchunlin/project/BJTC-302/00_rawdata/survival.xls")
group<-group[group$sample%in%sample$sample,]
group<-group[order(group$OS),]
group$OS<-ifelse(group$OS=='Survivor','Alive','Dead')
colnames(group)<-c('sample','group')
dat<-train.dat[hub_gene$hubgene,group$sample]%>%t%>%as.data.frame()
dat<-merge(dat,group,by.x='row.names',by.y='sample')%>%tibble::column_to_rownames(var = 'Row.names')
dat$group<-factor(dat$group,levels = c('Alive','Dead'))
colnames(dat)

# Lasso-logistic
##3
set.seed(4)
library(glmnet)

res.lasso <- cv.glmnet(as.matrix(dat[-ncol(dat)]), dat$group, family = "binomial")
plot(res.lasso)
plot(res.lasso$glmnet.fit, xvar = 'lambda')
ggsave("01.lasso.CV.png", plot(res.lasso), width = 6, height = 6, dpi = 300, units = "in", bg = "white")
ggsave("01.lasso.CV.pdf", plot(res.lasso), width = 6, height = 6, dpi = 300, units = "in", bg = "white")
ggsave("02.lasso.Coef.png", plot(res.lasso$glmnet.fit, xvar = 'lambda'), width = 6, height = 6, dpi = 300, units = "in", bg = "white")
ggsave("02.lasso.Coef.pdf", plot(res.lasso$glmnet.fit, xvar = 'lambda'), width = 6, height = 6, dpi = 300, units = "in", bg = "white")
l.coef<-coef(res.lasso$glmnet.fit,s=res.lasso$lambda.min,exact= F)
l.coef
res.lasso$lambda.min
## 0.01101454
colnames(dat)
#dat<-dat[,-c(3,4,6)]
set.seed(4)
res.lasso = glmnet(as.matrix(dat[-ncol(dat)]), dat$group, family = "binomial",  lambda = res.lasso$lambda.min)
saveRDS(list(res.lasso = res.lasso), "models.rds")
df.pred = predict.glmnet(res.lasso, newx = as.matrix(dat[-ncol(dat)]), type = "link")[,1]
df.pred = data.frame(sample = rownames(dat), score = df.pred)
df.pred$group = group$group[match(df.pred$sample, group$sample)]

# cs.en = confusion.glmnet(res.lasso, newx = as.matrix(dat[-ncol(dat)]), newy = dat$group) %>% as.matrix()
# cs.en = cs.en[2:1,2:1]
# cs.en
# cs.en <- cs.en %>% t %>% as.data.frame()
# # 真阳性放在第一位
# cs.en$True <- factor(cs.en$True, levels = c("Dead", "Burn"))
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
# ggsave(filename = "03.lasso.Confusion.png", height = 6, width = 6, p1)
# ggsave(filename = "03.lasso.Confusion.pdf", height = 6, width = 6, p1)

# library(PRROC)
# # roc.en = roc.curve(df.pred$score, weights.class0 = df.pred$group == "control", curve = T)
# # png("04.lasso.ROC.png", width = 4, height = 4, res = 300, units = "in", bg = "white")
# # plot(roc.en, auc.main = T, legend = F, color = F, xlab = "1-Specificity", asp = 1)
# # abline(0,1); dev.off()
# # pdf("04.lasso.ROC.pdf", width = 4, height = 4)
# # plot(roc.en, auc.main = T, legend = F, color = F, xlab = "1-Specificity", asp = 1)
# # abline(0,1); dev.off()
# 
# pr.en = pr.curve(df.pred$score, weights.class0 = df.pred$group == "Dead", curve = T)
# png("05.lasso.PR.png", width = 4, height = 4, res = 300, units = "in", bg = "white")
# plot(pr.en, auc.main = T, legend = F, color = 'red', asp = 1)
# dev.off()
# pdf("05.lasso.PR.pdf", width = 4, height = 4)
# plot(pr.en, auc.main = T, legend = F, color = 'red', asp = 1)
# dev.off()

library(pROC)
lasso_roc <- roc(df.pred$group, df.pred$score)
ggroc(lasso_roc,color = "red",
      linetype = 1,
      size = 1,
      alpha = 1,
      legacy.axes = T)+
  geom_abline(intercept = 0,
              slope = 1,
              color = "grey",
              size = 1,
              linetype = 1)+
  labs(x = "False Postive Rate(1 - Specificity)",
       y = "True Positive Rate(Sensivity or Recall)",
       title = "ROC curve")+
  annotate("text",x = 0.70,y = 0.30,
           label = paste("AUC =",signif(auc(lasso_roc),2)),
           size = 5,family = "Times")+
  theme_bw()+
  theme(panel.background = element_rect(fill = "transparent"),
        panel.grid = element_blank(),
        axis.ticks.length = unit(0.4,"lines"),
        axis.ticks = element_line(color = "black"),
        axis.line = element_line(size = 0.5,colour = "black"),
        axis.title = element_text(colour = "black",size = 15),
        axis.text = element_text(colour = "black",size = 10),
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        text = element_text(size = 8,color = "black",family = "Times"))
ggsave(filename = "03.lasso.ROC.png", width = 5, height = 5)
ggsave(filename = "03.lasso.ROC.pdf", width = 5, height = 5)


# library(ROCR)
# pred<-prediction(df.pred$score, df.pred$group)
# prc<-performance(pred,"prec", "rec")
# auc_ROCR <- performance(pred, measure = "auc")
# auc_ROCR <- auc_ROCR@y.values[[1]]
# png("05.lasso.PR.png", width = 4, height = 4, res = 300, units = "in", bg = "white", family = "Times")
# plot(prc,
#      xlim=c(0,1), ylim=c(0,1),
#      col='red',
#      #colorize=T, 
#      xlab= "recall",ylab="precision",
#      main='PR curve',
#      lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
# # abline(1,-1)
# text(x = 0.2, y = 0.3, label = paste0("AUC = ", signif(auc_ROCR,2)))
# dev.off()
# pdf("05.lasso.PR.pdf", width = 4, height = 4, family = "Times")
# plot(prc,
#      xlim=c(0,1), ylim=c(0,1),
#      col='red',
#      #colorize=T, 
#      xlab= "recall",ylab="precision",
#      main='PR curve',
#      lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
# # abline(1,-1)
# text(x = 0.2, y = 0.3, label = paste0("AUC = ", signif(auc_ROCR,2)))
# dev.off()

# https://www.icode9.com/content-4-1399138.html
# library(dcurves)
# df.pred$prob = (df.pred$score - min(df.pred$score))/(max(df.pred$score) - min(df.pred$score))
# dca(group ~ prob, data = df.pred,
#     label = list(prob = "Lasso Logistic")) %>% plot(smooth = T) +
#   scale_color_manual(values = c("#7FC97F","#BEAED4","#FDC086")) +
#   theme(aspect.ratio = 1, legend.position = "top",
#         text = element_text(size = 12, family = "Times"),
#         panel.grid = element_blank())
# ggsave("06.lasso.DCA.png", width = 5, height = 4, dpi = 300)
# ggsave("06.lasso.DCA.pdf", width = 5, height = 4, dpi = 300)


# colnames(dat)
# library(rmda)
# pdf("06.DCA.pdf",width=8,height=8)
# png("06.DCA.png",width=800,height=800)
# dat$group<-ifelse(dat$group=='PD',1,0)
# colnames(dat)
# complex<-decision_curve(group ~ DVL2+DNMT1+ABL1+RAF1+NOTCH1+RELA+PDGFRB,data = dat,family = binomial(link ='logit'),
#                         thresholds = seq(0,1, by = 0.01),
#                         confidence.intervals= 0.95,
#                         study.design = 'case-control',
#                         population.prevalence= 0.3
# )
# 
# 
# plot_decision_curve(complex,
#                     curve.names=c('logistic'),
#                     cost.benefit.axis =FALSE,col= c('black'),
#                     confidence.intervals=FALSE,
#                     standardize = FALSE)
# dev.off()
# 
# 

##基于多因素Logistic回归分析结果构建列线图
ddist <- datadist(dat)
options(datadist='ddist')
lrm <-lrm(group ~ S100A8+ITGAM, data=dat, x=TRUE, y=TRUE,maxit=1000)
print(lrm)
#nomogram计算部分，此处的f_lrm及对应的多因素logistic回归函数。
pdf("04.nomogram.pdf",width=6,height=4)
nomogram <- nomogram(lrm,fun=function(x)1/(1+exp(-x)), ##逻辑回归计算公式
                     fun.at = c(0.01,0.1,0.3,0.5,0.8,0.9,0.99),#风险轴刻度
                     funlabel = "Prob of death", #风险轴便签
                     lp=F,  ##是否显示系数轴
                     conf.int = F, ##每个得分的置信度区间，用横线表示,横线越长置信度越
                     abbrev = F#是否用简称代表因子变量
)
#绘制nomogram
plot(nomogram)
dev.off()
png("04.nomogram.png",width=600,height=400)
nomogram <- nomogram(lrm,fun=function(x)1/(1+exp(-x)), ##逻辑回归计算公式
                     fun.at = c(0.01,0.1,0.3,0.5,0.8,0.9,0.99),#风险轴刻度
                     funlabel = "Prob of death", #风险轴便签
                     lp=F,  ##是否显示系数轴
                     conf.int = F, ##每个得分的置信度区间，用横线表示,横线越长置信度越
                     abbrev = F#是否用简称代表因子变量
)
#绘制nomogram
plot(nomogram)
dev.off()

pdf("05.calibrate.pdf",width=8,height=6, family = "Times")
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
png("05.calibrate.png",width=8,height=6, units = "in", res = 600, family = "Times")
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


# 
# ### 单因素logistic------
# res.SPI1 <- glm(group ~ SPI1, data= dat, family = binomial())
# 
# res.S100A9 <- glm(group ~ S100A9, data= dat, family = binomial())
# res.S100A8 <- glm(group ~ S100A8, data= dat, family = binomial())
# res.TLR8 <- glm(group ~ TLR8, data= dat, family = binomial())
# res.S100A12 <- glm(group ~ S100A12, data= dat, family = binomial())
# res.ITGAM <- glm(group ~ ITGAM, data= dat, family = binomial())
# res.FCGR3B <- glm(group ~ FCGR3B, data= dat, family = binomial())
# res.CYBB <- glm(group ~ CYBB, data= dat, family = binomial())
# summary(res.CYBB)
# 
# ### 多因素logistic------
# mod <-glm(group ~ S100A8+ITGAM, data=dat, family = binomial())
# summary(mod)
# # lrm <-lrm(group ~ SPI1+TLR8+S100A9+S100A12+S100A8+FCGR3B+ITGAM+CYBB, data=dat, x=TRUE, y=TRUE,maxit=1000)
# # print(lrm)
# 
# ##基于多因素Logistic回归分析结果构建列线图
# ddist <- datadist(dat)
# options(datadist='ddist')
# #nomogram计算部分，此处的f_lrm及对应的多因素logistic回归函数。
# pdf(file=paste(output_dir, "\\nomogram.pdf", sep = ""),width=10,height=5) 
# nomogram <- nomogram(lrm,fun=function(x)1/(1+exp(-x)), ##逻辑回归计算公式
#                      fun.at = c(0.01,0.1,0.3,0.5,0.8,0.9,0.99),#风险轴刻度
#                      funlabel = "Prob of death", #风险轴便签
#                      lp=F,  ##是否显示系数轴
#                      conf.int = F, ##每个得分的置信度区间，用横线表示,横线越长置信度越
#                      abbrev = F#是否用简称代表因子变量
# )
# #绘制nomogram
# plot(nomogram)
# dev.off()
# 
# ###AUC###
# colnames(dat)
# library(pROC)
# pridictdata<-dat[,c('S100A8','ITGAM')]
# df.pred<-predict(mod, newdata=pridictdata,type = 'link')
# df.pred = data.frame(sample = rownames(dat), score = df.pred)
# df.pred$group<-factor(group$group,levels = c('Dead','Alive'))
# roc <- roc(df.pred$group, df.pred$score,)
# plot(roc,
#      print.auc=T,
#      print.auc.x=0.4,print.auc.y=0.5,
#      #auc.polygon=T,
#      #auc.polygon.con="#fff7f7",
#      grid=c(0.5,0.2),
#      grid.col=c("black","black"),
#      print.thres=T,
#      main="ROC curve",
#      col="#FF2E63",
#      legacy.axes=T)
