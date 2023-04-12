rm(list = ls())
setwd("/data/nas1/luchunlin/project/NN-0118-2/")
if (! dir.exists("./08_randomforest")){
  dir.create("./08_randomforest")
}
setwd("./08_randomforest")

library(magrittr)
library(ggplot2)
hub_gene <- read.delim2('../06_PPI/hubgene.xls')
train.dat<-read.delim2("../00_rawdata/dat.fpkm.xls", row.names = 1)  %>% lc.tableToNum
colnames(train.dat)<-gsub('.','-',colnames(train.dat),fixed = T)
train.dat <- log2(train.dat+1)
group = read.delim2("../01_DEGs(TCGA)/group.xls")
dat<-train.dat[hub_gene$symbol,group$sample]%>%t%>%as.data.frame()
dat<-merge(dat,group,by.x='row.names',by.y='sample')%>%tibble::column_to_rownames(var = 'Row.names')
table(group$group)
dat$group<-factor(dat$group,levels = c('Tumor','Normal'))
# 随机森林
library(randomForest)
set.seed(123)
##38
res.rf <- randomForest(group ~ ., data = dat,ntrees=1000,mtry = 2,importance = T)
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
data$true_label<-factor(data$true_label,levels = c('Normal','Tumor'))
#data$pred_label<-factor(data$pred_label,levels = c('asthma','control'))
p1 <- ggplot(data, aes(x = pred_label, y = true_label, fill = value)) +
  geom_tile() +
  geom_text(aes(label = value), family = "Times", size = 6) +
  scale_x_discrete(breaks = c("pred_normal", "pred_patient"),
                   labels = c("control", "asthma")) +
  scale_fill_gradient(low = "white", high = "royalblue") +
  coord_fixed() +
  labs(x="Predicted Classification", y="True Classification") +
  theme_minimal() +
  theme(axis.text.x = element_text(size = 13, face = "bold", family = "Times"),
        axis.text.y = element_text(size = 13, face = "bold", family = "Times"),
        axis.title = element_text(size = 15, face = "bold", family = "Times"),
        legend.position = "none")
p1
# ggsave(filename = "02.RF.Confusion.png", height = 6, width = 6, p1)
# ggsave(filename = "02.RF.Confusion.pdf", height = 6, width = 6, p1)
#根据随机森林中的不同数来绘制误差率

#变量重要性评分
importance(res.rf,type = 1)
#重要性绘图
pdf(file = '02.importance.pdf',w=5,h=6)
varImpPlot(res.rf,main = 'importance')
dev.off()
png(file = '02.importance.png',w=400,h=500)
varImpPlot(res.rf,main = 'importance')
dev.off()

