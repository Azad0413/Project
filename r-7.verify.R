rm(list = ls())
# Logistic----------
setwd("/data/nas1/luchunlin/project/BJTC-317")
if (! dir.exists("./07_verify")){
  dir.create("./07_verify")
}
setwd("./07_verify")

### 验证集1-------
gset<-getGEO("GSE20163",destdir = '.',GSEMatrix = T,getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
a=gset[[1]]
gpl<-getGEO("GPL96",destdir = '.')
gpl<-Table(gpl)    
colnames(gpl)
probe2symbol<-gpl %>%
  dplyr::select('ID','Gene Symbol')%>%
  filter('Gene Symbol'!='')%>%
  separate('Gene Symbol',c('symbol','drop'),sep = '///')%>%
  dplyr::select(-drop)
dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symbol$ID<-as.character(probe2symbol$ID)
dat<-dat %>%
  inner_join(probe2symbol,by='ID')%>% 
  dplyr::select(-ID)%>%     ## 去除多余信息
  dplyr::select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
pd<-pData(a)
group<-data.frame(sample=pd$geo_accession,
                  group=pd$source_name_ch1)
table(group$group)
group$group<-ifelse(group$group=='SN, control','control','PD')
table(group$group)
dat.va<-dat[,group$sample]
write.table(dat.va,file = 'dat_va1.xls',sep = '\t',row.names = T,quote = F)
write.table(group,file = 'group_va1.xls',sep = '\t',row.names = F,quote = F)
hub_gene <- read.delim2('/data/nas1/luchunlin/project/BJTC-317/06_Logistic/hubgene.xls')
test.dat<-dat.va[hub_gene$sybmol,]%>%t%>%as.data.frame()
colnames(test.dat)
#test.dat<-test.dat[,-c(3,4,6)]
models <- readRDS("../06_Logistic/models.rds")
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
       title = "ROC curve (GSE20163)")+
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
ggsave(filename = "ROC.verify1.png", width = 5, height = 5)
ggsave(filename = "ROC.verify1.pdf", width = 5, height = 5)
# 验证集2--------

gset2<-getGEO("GSE7621",destdir = '.',GSEMatrix = T,getGPL = F)
expr2<-as.data.frame(exprs(gset2[[1]]))
a2=gset2[[1]]
gpl2<-getGEO("GPL570",destdir = '.')
gpl2<-Table(gpl2)    
colnames(gpl2)
probe2symbol2<-gpl2 %>%
  dplyr::select('ID','Gene Symbol')%>%
  filter('Gene Symbol'!='')%>%
  separate('Gene Symbol',c('symbol','drop'),sep = '///')%>%
  dplyr::select(-drop)
dat2<-expr2
dat2$ID<-rownames(dat2)
dat2$ID<-as.character(dat2$ID)
probe2symbol2$ID<-as.character(probe2symbol2$ID)
dat2<-dat2 %>%
  inner_join(probe2symbol2,by='ID')%>% 
  dplyr::select(-ID)%>%     ## 去除多余信息
  dplyr::select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
pd2<-pData(a2)
group2<-data.frame(sample=pd2$geo_accession,
                  group=pd2$characteristics_ch1)
table(group2$group)
group2$group<-ifelse(group2$group=='Old Control','control','PD')
table(group2$group)
dat.va2<-dat2[,group2$sample]
write.table(dat.va,file = 'dat_va2.xls',sep = '\t',row.names = T,quote = F)
write.table(group,file = 'group_va2.xls',sep = '\t',row.names = F,quote = F)
## 
test.dat2<-dat.va2[hub_gene$sybmol,]%>%t%>%as.data.frame()
colnames(test.dat2)
#test.dat2<-test.dat2[,-c(3,4,6)]
models <- readRDS("../06_Logistic/models.rds")
df.pred2 <- predict(models$res.lasso, newx = as.matrix(test.dat2), type = "link") %>% as.data.frame %>% 
  tibble::rownames_to_column(var = "sample")
colnames(df.pred2)[2] <- "score"
df.pred2$group = group2$group[match(df.pred2$sample, group2$sample)]
library(pROC)
lasso_roc2 <- roc(df.pred2$group, df.pred2$score)
ggroc(lasso_roc2,color = "red",
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
       title = "ROC curve (GSE7621)")+
  annotate("text",x = 0.70,y = 0.30,
           label = paste("AUC =", signif(auc(lasso_roc2),2)),
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
ggsave(filename = "ROC.verify2.png", width = 5, height = 5)
ggsave(filename = "ROC.verify2.pdf", width = 5, height = 5)

## 验证集3 弃-----
# gset<-getGEO("GSE99039",destdir = '.',GSEMatrix = T,getGPL = F)
# expr<-as.data.frame(exprs(gset[[1]]))
# a=gset[[1]]
# gpl<-getGEO("GPL570",destdir = '.')
# gpl<-Table(gpl)    
# colnames(gpl)
# probe2symbol<-gpl %>%
#   dplyr::select('ID','Gene Symbol')%>%
#   filter('Gene Symbol'!='')%>%
#   separate('Gene Symbol',c('symbol','drop'),sep = '///')%>%
#   dplyr::select(-drop)
# dat<-expr
# dat$ID<-rownames(dat)
# dat$ID<-as.character(dat$ID)
# probe2symbol$ID<-as.character(probe2symbol$ID)
# dat<-dat %>%
#   inner_join(probe2symbol,by='ID')%>% 
#   dplyr::select(-ID)%>%     ## 去除多余信息
#   dplyr::select(symbol,everything())%>%     ## 重新排列
#   mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
#   arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
#   distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
#   dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
#   tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
# pd<-pData(a)
# group<-data.frame(sample=pd$geo_accession,
#                   group=pd$characteristics_ch1.2)
# table(group$group)
# group<-group[which(group$group==c('disease label: CONTROL','disease label: IPD')),]
# group$group<-ifelse(group$group=='disease label: CONTROL','control','PD')
# table(group$group)
# dat.va<-dat[,group$sample]
# write.table(dat.va,file = 'dat_va.xls',sep = '\t',row.names = T,quote = F)
# write.table(group,file = 'group_va.xls',sep = '\t',row.names = F,quote = F)
# ## 
# hub_gene <- read.delim2('/data/nas1/luchunlin/project/BJTC-317/06_Logistic/hubgene.xls')
# test.dat<-dat.va[hub_gene$sybmol,]%>%t%>%as.data.frame()
# test.dat<-test.dat[,-c(3,4,6)]
# models <- readRDS("../06_Logistic/models.rds")
# df.pred <- predict(models$res.lasso, newx = as.matrix(test.dat), type = "link") %>% as.data.frame %>% 
#   tibble::rownames_to_column(var = "sample")
# colnames(df.pred)[2] <- "score"
# df.pred$group = group$group[match(df.pred$sample, group$sample)]
# library(pROC)
# lasso_roc <- roc(df.pred$group, df.pred$score)
# ggroc(lasso_roc,color = "red",
#       linetype = 1,
#       size = 1,
#       alpha = 1,
#       legacy.axes = T)+
#   geom_abline(intercept = 0,
#               slope = 1,
#               color = "grey",
#               size = 1,
#               linetype = 1)+
#   labs(x = "False Postive Rate(1 - Specificity)",
#        y = "True Positive Rate(Sensivity or Recall)",
#        title = "ROC curve (GSE99039)")+
#   annotate("text",x = 0.70,y = 0.30,
#            label = paste("AUC =", signif(auc(lasso_roc),2)),
#            size = 5,family = "Times")+
#   theme_bw()+
#   theme(panel.background = element_rect(fill = "transparent"),
#         panel.grid = element_blank(),
#         axis.ticks.length = unit(0.4,"lines"),
#         axis.ticks = element_line(color = "black"),
#         axis.line = element_line(size = 0.5,colour = "black"),
#         axis.title = element_text(colour = "black",size = 15),
#         axis.text = element_text(colour = "black",size = 10),
#         plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
#         text = element_text(size = 8,color = "black",family = "Times"))
# ggsave(filename = "01.lasso.ROC.png", width = 5, height = 5)
# ggsave(filename = "01.lasso.ROC.pdf", width = 5, height = 5)
