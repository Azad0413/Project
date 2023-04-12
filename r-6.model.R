rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-385-10/")
if (! dir.exists("./06_model")){
  dir.create("./06_model")
}
setwd("./06_model")
library(magrittr)
library(tidyverse)
library(lance)
library(ggplot2)
hub_gene <- read.delim2('../03_intersect/intersect.xls')
train.dat<-read.delim2("../01_DEGs/dat_final.xls", row.names = 1)  %>% lc.tableToNum
group = read.delim2("../01_DEGs/group.xls")
dat<-train.dat[hub_gene$symbol,group$sample]%>%t%>%as.data.frame()
dat<-merge(dat,group,by.x='row.names',by.y='sample')%>%tibble::column_to_rownames(var = 'Row.names')
dat$group<-factor(dat$group,levels = c('RB','control'))
# 随机森林
library(randomForest)
set.seed(11)
res.rf <- randomForest(group ~ ., data = dat,ntrees=500,mtry = 2,importance = T)
plot(res.rf, main = NULL)
png("01.RF.ntree.png", width = 4, height = 4, bg = "white", units = "in", res = 300, family = "Times")
plot(res.rf, main = NULL)
dev.off()
pdf("01.RF.ntree.pdf", width = 4, height = 4, family = "Times")
plot(res.rf, main = NULL)
dev.off()
saveRDS(res.rf, "rf.models.rds")
##混淆矩阵----
ntree = which.min(res.rf$err.rate[,1])
res.rf = randomForest::randomForest(group ~ ., data = dat, ntree = ntree)
df.rf = data.frame(sample = rownames(res.rf$votes), vote = res.rf$votes[,2])
df.rf$group = group$group[match(df.rf$sample, group$sample)]
data <- res.rf$confusion[2:1,2:1] %>% t %>% as.data.frame() %>% tibble::rownames_to_column(var = "true_label")
colnames(data)[2:3] <- c("pred_patient", "pred_normal")
data <- tidyr::gather(data, pred_label, value, -1)
data$true_label<-factor(data$true_label,levels = c('control','RB'))
#data$pred_label<-factor(data$pred_label,levels = c('RB','control'))
p1 <- ggplot(data, aes(x = pred_label, y = true_label, fill = value)) +
  geom_tile() +
  geom_text(aes(label = value), family = "Times", size = 6) +
  scale_x_discrete(breaks = c("pred_normal", "pred_patient"),
                   labels = c("RB","control")) +
  scale_fill_gradient(low = "white", high = "royalblue") +
  coord_fixed() +
  labs(x="Predicted Classification", y="True Classification") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 13, face = "bold", family = "Times"),
        axis.text.y = element_text(size = 13, face = "bold", family = "Times"),
        axis.title = element_text(size = 15, face = "bold", family = "Times"),
        legend.position = "none")
p1
ggsave(filename = "02.RF.Confusion.png", height = 6, width = 6, p1)
ggsave(filename = "02.RF.Confusion.pdf", height = 6, width = 6, p1)

#根据随机森林中的不同数来绘制误差率

#变量重要性评分
importance(res.rf,type = 1)
#重要性绘图
pdf(file = '03.RF.importance.pdf',w=5,h=6.5)
varImpPlot(res.rf,main = 'importance')
dev.off()
png(file = '03.RF.importance.png',w=400,h=550)
varImpPlot(res.rf,main = 'importance')
dev.off()
importance <- data.frame(res.rf$importance)%>%rownames_to_column(var='symbol')
importance <- importance[order(importance$MeanDecreaseGini,decreasing = T),]
importance <- importance[c(1:10),]
write.table(importance,'importance.xls',sep = '\t',row.names = F,quote = F)
#使用测试集对构建好的随机森林进行测试
irispred<-predict(res.rf,newdata = dat)
table(irispred,dat$group)
#数据点的边距为正确归类的比例减去被归到其他类的最大比例  #观测值被判断正确的概率图
pdf('04.RF.probality.pdf',w=6,h=6)
plot(margin(res.rf,dat$group),main='Probality graph of the observed value being judged correct',
     ylab='Probality',
     xlab='Predict the number of samples')
dev.off()
png('04.RF.probality.png',w=600,h=600)
plot(margin(res.rf,dat$group),main='Probality graph of the observed value being judged correct',
     ylab='Probality',
     xlab='Predict the number of samples')
dev.off()

# ## PR-----
# library(PRROC)
# pr.en = pr.curve(df.rf$vote, weights.class0 = df.rf$group == "control", curve = T)
# png("04.RF.PR.png", width = 4, height = 4, res = 300, units = "in", bg = "white")
# plot(pr.en, auc.main = T, legend = F, color = 'red', asp = 1)
# dev.off()
# pdf("04.RF.PR.pdf", width = 4, height = 4)
# plot(pr.en, auc.main = T, legend = F, color = 'red', asp = 1)
# dev.off()

## roc------
library(pROC)
rf_roc <- roc(df.rf$group, df.rf$vote)
ggroc(rf_roc,color = "red",
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
       title = "RF ROC curve")+
  annotate("text",x = 0.70,y = 0.30,
           label = paste("AUC =",signif(auc(rf_roc),2)),
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
ggsave(filename = "04.RF.ROC.png", width = 5, height = 5)
ggsave(filename = "04.RF.ROC.pdf", width = 5, height = 5)

# ## DCA-----
# library(dcurves)
# df.rf$prob = (df.rf$vote - min(df.rf$vote))/(max(df.rf$vote) - min(df.rf$vote))
# dca(group ~ prob, data = df.rf,
#     label = list(prob = "RF")) %>% plot(smooth = T) +
#   scale_color_manual(values = c("#7FC97F","#BEAED4","#FDC086")) +
#   theme(aspect.ratio = 1, legend.position = "top",
#         text = element_text(size = 12, family = "Times"),
#         panel.grid = element_blank())
# ggsave("06.RF.DCA.png", width = 5, height = 4, dpi = 300)
# ggsave("06.RF.DCA.pdf", width = 5, height = 4, dpi = 300)

# Lasso-logistic---------
set.seed(123)
library(glmnet)
res.lasso <- cv.glmnet(as.matrix(dat[-ncol(dat)]), dat$group, family = "binomial", 
                       type.measure = "auc")

plot(res.lasso)
plot(res.lasso$glmnet.fit, xvar = 'lambda')
ggsave("01.lasso.CV.png", plot(res.lasso), width = 6, height = 5, dpi = 300, units = "in", bg = "white")
ggsave("01.lasso.CV.pdf", plot(res.lasso), width = 6, height = 5, dpi = 300, units = "in", bg = "white")
ggsave("02.lasso.Coef.png", plot(res.lasso$glmnet.fit, xvar = 'lambda'), width = 6, height = 5, dpi = 300, units = "in", bg = "white")
ggsave("02.lasso.Coef.pdf", plot(res.lasso$glmnet.fit, xvar = 'lambda'), width = 6, height = 5, dpi = 300, units = "in", bg = "white")
l.coef<-coef(res.lasso$glmnet.fit,s=res.lasso$lambda.min,exact= F)
l.coef
res.lasso$lambda.min
#0.06493903
##PDE8B         0.1790931
##ESRRB         1.5516442
##SPRY2         0.5833901
lasso_geneids <- l.coef@Dimnames[[1]][l.coef@i+1]
lasso_geneids <- lasso_geneids[-1]
lasso_geneids
write(lasso_geneids, "lasso_genes.csv")
set.seed(123)
res.lasso = glmnet(as.matrix(dat[-ncol(dat)]), dat$group, family = "binomial",  lambda = res.lasso$lambda.min)
saveRDS(list(res.lasso = res.lasso), "lasso.models.rds")

df.pred = predict.glmnet(res.lasso, newx = as.matrix(dat[-ncol(dat)]), type = "link")[,1]
df.pred = data.frame(sample = rownames(dat), score = df.pred)
df.pred$group = group$group[match(df.pred$sample, group$sample)]

cs.en = confusion.glmnet(res.lasso, newx = as.matrix(dat[-ncol(dat)]), newy = dat$group) %>% as.matrix()
cs.en = cs.en[2:1,2:1]
cs.en
cs.en <- cs.en %>% t %>% as.data.frame()
# 真阳性放在第一位
cs.en$True <- factor(cs.en$True, levels = c("control","RB"))
cs.en$Predicted<-factor(cs.en$Predicted,levels = c('RB','control'))
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

# library(PRROC)
# plot(pr.en, auc.main = T, legend = F, color = 'red', asp = 1)
# pr.en = pr.curve(df.pred$score, weights.class0 = df.pred$group == "control", curve = T)
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
       title = "lasso ROC curve")+
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


hub_gene <- intersect(lasso_geneids,importance$symbol)%>%as.data.frame()
colnames(hub_gene) <- 'symbol'
write.table(hub_gene,file = 'hubgene.xls',sep = '\t',row.names = F,quote = F)
library(ggvenn)
mydata<-list('Lasso-logistic'=lasso_geneids,RandomForest=importance$symbol)
pdf('06.venn.pdf',w=6,h=6)
ggvenn(mydata,c('Lasso-logistic','RandomForest'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       text_size = 6,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()
png('06.venn.png',w=400,h=400)
ggvenn(mydata,c('Lasso-logistic','RandomForest'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       text_size = 6,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()
