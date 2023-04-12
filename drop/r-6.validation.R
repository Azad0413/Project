rm(list = ls())
setwd("/data/nas1/luchunlin/project/TY0307-11/")
if (! dir.exists("./06_validation")){
  dir.create("./06_validation")
}
setwd("./06_validation")
hubgene <- read.delim2('../04_model/hubgene.xls')
# hubgene <- read.delim2('../02_DEERS/DEERS.xls')
dat<-read.delim2("../00_rawdata/dat(GSE42148).xls", row.names = 1)%>% lc.tableToNum
group = read.delim2("../00_rawdata/group(GSE42148).xls")
control.sample<-group$sample[which(group$group=='control')]
## 
# LOC442546
hub_exp<-dat[hubgene$x,]
#unmap <- dat['LOC442546',]
#rownames(unmap) <- 'NUPR1'
#hub_exp <- rbind(unmap,hub_exp)
hub_exp2<-hub_exp
hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))

## 样本分组
hub_exp2$Group<-ifelse(hub_exp2$sample%in%control.sample,'control','CAD')

##分面图形
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(rstatix)
stat.test<-hub_exp2%>%
  group_by(Symbol)%>%
  wilcox_test(expr ~ Group)%>%
  adjust_pvalue(method = 'fdr')
stat.test$p<-ifelse(stat.test$p<0.001,"***",
                        ifelse(stat.test$p<0.05,"**",
                               ifelse(stat.test$p<0.05,"*",'ns')))
exp_plot <- ggplot(hub_exp2,aes(x = Group, y = expr, color = Group)) +
  #geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#4682B4","#CD3700"), name = "Group")+
  labs(title="Expression", x="", y = "",size=20) +
  stat_pvalue_manual(stat.test,
                     y.position = c(10,6,10,6.5,10,10),
                     size = 3.2,
                     family = "Times",
                     label = "p",
                     #parse = T,
                     face = "bold")+
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=15),
        axis.text.x=element_text(angle=0,hjust=0.5,colour="black",face="bold",size=12), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill='none')+
  facet_wrap(~Symbol,scales = "free",nrow = 2) 
exp_plot
ggsave('01.expression.pdf',exp_plot,width = 7,height = 4)
ggsave('01.expression.png',exp_plot,width = 7,height = 4)

test.dat<-hub_exp%>%t%>%as.data.frame()
colnames(test.dat)
test.dat$group<-factor(group$group)
## LASSO-----
models <- readRDS("../04_model/lasso.models.rds")

df.pred <- predict(models$res.lasso, newx = as.matrix(test.dat[-ncol(test.dat)]), type = "link") %>% as.data.frame %>%
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
       title = "Validation ROC curve")+
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
ggsave(filename = "02.ROC.verify.png", width = 5, height = 5)
ggsave(filename = "02.ROC.verify.pdf", width = 5, height = 5)



group$group = factor(group$group, levels = c("CAD", "control"))
library(ROCR)
library(ggplot2)
hub_exp2<-t(hub_exp)
## 绘制ROC曲线
library(pROC)
for (i in c(1:24)) {
  roc<-roc(group$group,hub_exp2[,i],levels=c("CAD", "control"))
#  png(paste0(i, ".", colnames(hub_exp2)[i],".png"),width = 300,height = 300)
  plot(roc,
       print.auc=T,
       print.auc.x=0.4,print.auc.y=0.5,
       #auc.polygon=T,
       #auc.polygon.con="#fff7f7",
       grid=c(0.5,0.2),
       grid.col=c("black","black"),
       #print.thres=T,
       main=colnames(hub_exp2)[i],
       col="#FF2E63",
       legacy.axes=T)
  # dev.off()
  # pdf(paste0(i, ".", colnames(hub_exp2)[i],".pdf"),width = 4,height = 4)
  # plot(roc,
  #      print.auc=T,
  #      print.auc.x=0.4,print.auc.y=0.5,
  #      #auc.polygon=T,
  #      #auc.polygon.con="#fff7f7",
  #      grid=c(0.5,0.2),
  #      grid.col=c("black","black"),
  #      #print.thres=T,
  #      main=colnames(hub_exp2)[i],
  #      col="#FF2E63",
  #      legacy.axes=T)
  # dev.off()
  i<-i+1
}

## PCA-------
cluster_exp<-t(hub_exp)
pca<-prcomp(cluster_exp,center = TRUE,scale. = TRUE)
df<-pca$x  ## 提取PC score
df <- as.data.frame(df)
summ <- summary(pca)
summ
# 提取主成分的方差贡献率,生成坐标轴标题
summ <- summary(pca)
summ
xlab <- paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")
ylab <- paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")
# 绘制PCA得分图
table(group$group)
library(ggplot2)
p.pca <- ggplot(data = df,aes(x = PC1,y = PC2,color = group$group))+
  stat_ellipse(aes(fill = group$group),
               type = "norm",geom = "polygon",alpha = 0.25,color = NA)+ # 添加置信椭圆
  geom_point(size = 3.5)+
  labs(x = xlab,y = ylab,color = "Group",title = "PCA Scores Plot")+
  guides(fill = "none")+
  theme_bw()+
  scale_fill_manual(values = c("purple","orange","pink"))+
  scale_colour_manual(values = c("purple","orange","pink"))+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text = element_text(size = 11),axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),legend.title = element_text(size = 13),
        plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))
p.pca
ggsave('01.pca.pdf',p.pca,w=6,h=5)
ggsave('01.pca.png',p.pca,w=6,h=5)


