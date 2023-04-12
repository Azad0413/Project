rm(list = ls())
# Logistic----------
setwd("/data/nas1/luchunlin/project/BJTC-302")
if (! dir.exists("./07_verify")){
  dir.create("./07_verify")
}
setwd("./07_verify")

### 验证集-------
dat<-read.delim2('/data/nas1/luchunlin/project/BJTC-302/10_validation/dat_va.xls',row.names = 1)%>%lc.tableToNum()
group<-read.delim2('/data/nas1/luchunlin/project/BJTC-302/10_validation/group.va.xls')
table(group$group)
dat.va<-dat[,group$sample]
hub_gene<-read.delim2('/data/nas1/luchunlin/project/BJTC-302/10_validation/hubgene.final.xls')
test.dat<-dat.va[hub_gene$hubgene,]%>%t%>%as.data.frame()
colnames(test.dat)
#test.dat<-test.dat[,-c(3,4,6)]
models <- readRDS("../06_logistic/models.rds")
df.pred <- predict(models$res.lasso, newx = as.matrix(test.dat), type = "link") %>% as.data.frame %>% 
  tibble::rownames_to_column(var = "sample")
colnames(df.pred)[2] <- "score"
df.pred$group = group$group[match(df.pred$sample, group$sample)]
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
       title = "ROC curve (GSE37069)")+
  annotate("text",x = 0.70,y = 0.30,
           label = paste("AUC =", signif(auc(lasso_roc),2)),
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
ggsave(filename = "ROC.verify.png", width = 5, height = 5)
ggsave(filename = "ROC.verify.pdf", width = 5, height = 5)
