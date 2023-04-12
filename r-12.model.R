rm(list = ls())
# 01 获取数据集--------------
setwd("/data/nas1/luchunlin/project/BJTC-356/")
if (! dir.exists("./11_model")){
  dir.create("./11_model")
}
setwd("./11_model")
library(magrittr)
library(ggplot2)
hub_gene <- read.delim2('../06_quadrant/final.gene.xls')
train.dat<-read.delim2("../00_rawdata/dat.fpkm.xls", row.names = 1)  %>% lc.tableToNum
group = read.delim2("../01_GSVA/group.xls")
dat<-train.dat[hub_gene$symbol,group$sample]%>%t%>%as.data.frame()
dat <- log2(dat+1)
dat<-merge(dat,group,by.x='row.names',by.y='sample')%>%tibble::column_to_rownames(var = 'Row.names')
table(group$group)
dat$group<-factor(dat$group,levels = c('Persistent','Paroxysmal'))
colnames(dat) <- gsub('-','_',colnames(dat),fixed = T)

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
#0.00156685
colnames(dat)
#dat<-dat[,-c(3,4,6)]
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
cs.en$True <- factor(cs.en$True, levels = c("Paroxysmal","Persistent"))
cs.en$Predicted<-factor(cs.en$Predicted,levels = c('Persistent','Paroxysmal'))
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
ggsave(filename = "04.lasso.ROC.png", width = 5, height = 5)
ggsave(filename = "04.lasso.ROC.pdf", width = 5, height = 5)

##随机森林-------
library(randomForest)
set.seed(36)
##6
res.rf <- randomForest(group ~ ., data = dat,ntrees=1000,mtry = 2,importance = T)
plot(res.rf, main = NULL)
png("05.RF.ntree.png", width = 4, height = 4, bg = "white", units = "in", res = 300, family = "Times")
plot(res.rf, main = NULL)
dev.off()
pdf("05.RF.ntree.pdf", width = 4, height = 4, family = "Times")
plot(res.rf, main = NULL)
dev.off()
saveRDS(res.rf, "rf.models.rds")
##混淆矩阵----
ntree = which.min(res.rf$err.rate[,1])
res.rf = randomForest::randomForest(group ~ ., data = dat, ntree = ntree)
df.rf = data.frame(sample = rownames(res.rf$votes), vote = res.rf$votes[,2])
df.rf$group = group$group[match(df.rf$sample, group$sample)]
data <- res.rf$confusion[1:2,1:2] %>% t %>% as.data.frame() %>% tibble::rownames_to_column(var = "true_label")
res.rf$confusion
data <- tidyr::gather(data, pred_label, value, -1)
data$true_label<-factor(data$true_label,levels = c('Paroxysmal','Persistent'))
data$pred_label <- factor(data$pred_label,levels = c('Persistent','Paroxysmal'))

p1 <- ggplot(data, aes(x = pred_label, y = true_label, fill = value)) +
  geom_tile() +
  geom_text(aes(label = value), family = "Times", size = 6) +
  scale_x_discrete(breaks = c('Persistent','Paroxysmal'),
                   labels = c("Persistent","Paroxysmal")) +
  scale_fill_gradient(low = "white", high = "royalblue") +
  coord_fixed() +
  labs(x="Predicted Classification", y="True Classification") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 13, face = "bold", family = "Times"),
        axis.text.y = element_text(size = 13, face = "bold", family = "Times"),
        axis.title = element_text(size = 15, face = "bold", family = "Times"),
        legend.position = "none")
p1
ggsave(filename = "06.RF.Confusion.png", height = 6, width = 6, p1)
ggsave(filename = "06.RF.Confusion.pdf", height = 6, width = 6, p1)

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
ggsave(filename = "07.RF.ROC.png", width = 5, height = 5)
ggsave(filename = "07.RF.ROC.pdf", width = 5, height = 5)


