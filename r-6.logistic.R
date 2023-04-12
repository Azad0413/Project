rm(list = ls())
# Logistic----------
setwd("/data/nas1/luchunlin/project/BJTC-317")
if (! dir.exists("./06_Logistic")){
  dir.create("./06_Logistic")
}
setwd("./06_Logistic")
library(magrittr)
library(ggplot2)
hub_gene <- data.frame(sybmol=c('DVL2','DNMT1','ABL1','RAF1','NOTCH1','RELA','PDGFRB'))
write.table(hub_gene,file = 'hubgene.xls',sep = '\t',row.names = F,quote = F)
train.dat<-read.delim2("/data/nas1/luchunlin/project/BJTC-317/00_rawdata/dat.final.xls", row.names = 1)  %>% lc.tableToNum
group = read.delim2("/data/nas1/luchunlin/project/BJTC-317/00_rawdata/group.xls")
dat<-train.dat[hub_gene$sybmol,group$sample]%>%t%>%as.data.frame()
dat<-log2(dat+0.0001)
dat<-merge(dat,group,by.x='row.names',by.y='sample')%>%tibble::column_to_rownames(var = 'Row.names')
dat$group<-factor(dat$group,levels = c('PD','control'))

# Lasso-logistic
set.seed(123)
library(glmnet)
res.lasso <- cv.glmnet(as.matrix(dat[-ncol(dat)]), dat$group, family = "binomial", 
                       type.measure = "auc")

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
set.seed(123)
res.lasso = glmnet(as.matrix(dat[-ncol(dat)]), dat$group, family = "binomial",  lambda = res.lasso$lambda.min)
saveRDS(list(res.lasso = res.lasso), "models.rds")

df.pred = predict.glmnet(res.lasso, newx = as.matrix(dat[-ncol(dat)]), type = "link")[,1]
df.pred = data.frame(sample = rownames(dat), score = df.pred)
df.pred$group = group$group[match(df.pred$sample, group$sample)]

cs.en = confusion.glmnet(res.lasso, newx = as.matrix(dat[-ncol(dat)]), newy = dat$group) %>% as.matrix()
cs.en = cs.en[2:1,2:1]
cs.en
cs.en <- cs.en %>% t %>% as.data.frame()
# 真阳性放在第一位
cs.en$True <- factor(cs.en$True, levels = c("control", "PD"))
p1 <- ggplot(cs.en, aes(x = Predicted, y = True, fill = Freq)) +
  geom_tile(color = "grey") +
  geom_text(aes(label = Freq), family = "Times", size = 6) +
  scale_fill_gradient(low = "white", high = "royalblue") +
  coord_fixed() +
  labs(x="Predicted Classification", y="True Classification") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 13, face = "bold", family = "Times"),
        axis.text.y = element_text(size = 13, face = "bold", family = "Times"),
        axis.title = element_text(size = 15, face = "bold", family = "Times"),
        legend.position = "none",
        panel.grid = element_blank())
p1
ggsave(filename = "03.lasso.Confusion.png", height = 6, width = 6, p1)
ggsave(filename = "03.lasso.Confusion.pdf", height = 6, width = 6, p1)

library(PRROC)
# roc.en = roc.curve(df.pred$score, weights.class0 = df.pred$group == "control", curve = T)
# png("04.lasso.ROC.png", width = 4, height = 4, res = 300, units = "in", bg = "white")
# plot(roc.en, auc.main = T, legend = F, color = F, xlab = "1-Specificity", asp = 1)
# abline(0,1); dev.off()
# pdf("04.lasso.ROC.pdf", width = 4, height = 4)
# plot(roc.en, auc.main = T, legend = F, color = F, xlab = "1-Specificity", asp = 1)
# abline(0,1); dev.off()

pr.en = pr.curve(df.pred$score, weights.class0 = df.pred$group == "control", curve = T)
png("05.lasso.PR.png", width = 4, height = 4, res = 300, units = "in", bg = "white")
plot(pr.en, auc.main = T, legend = F, color = 'red', asp = 1)
dev.off()
pdf("05.lasso.PR.pdf", width = 4, height = 4)
plot(pr.en, auc.main = T, legend = F, color = 'red', asp = 1)
dev.off()

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
ggsave(filename = "04.lasso.ROC.png", width = 5, height = 5)
ggsave(filename = "04.lasso.ROC.pdf", width = 5, height = 5)


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
# 



colnames(dat)
library(rmda)
pdf("06.DCA.pdf",width=8,height=8)
png("06.DCA.png",width=800,height=800)
dat$group<-ifelse(dat$group=='PD',1,0)
colnames(dat)
complex<-decision_curve(group ~ DVL2+DNMT1+ABL1+RAF1+NOTCH1+RELA+PDGFRB,data = dat,family = binomial(link ='logit'),
                        thresholds = seq(0,1, by = 0.01),
                        confidence.intervals= 0.95,
                        study.design = 'case-control',
                        population.prevalence= 0.3
)


plot_decision_curve(complex,
                    curve.names=c('logistic'),
                    cost.benefit.axis =FALSE,col= c('black'),
                    confidence.intervals=FALSE,
                    standardize = FALSE)
dev.off()
