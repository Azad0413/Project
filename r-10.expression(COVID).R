rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/BJTC-420-1/")
if (! dir.exists("./10_expression(COVID)")){
  dir.create("./10_expression(COVID)")
}
setwd("./10_expression(COVID)")
dat<-read.delim2("../00_rawdata/dat(GSE171110).xls", row.names = 1)%>% lc.tableToNum
group <- read.delim2('../00_rawdata/group(GSE171110).xls')
table(group$group)
hubgene <- read.delim2('../05_PPI/hubgene.xls')

diff <- read.delim2('../01_DEGs(COVID)/DEG_sig.xls',row.names = 1)
#1 train-----
hub_exp<-dat[hubgene$symbol,]
hub_exp <- edgeR::cpm(hub_exp)%>%as.data.frame()
hub_exp2<-log2(hub_exp+1)
hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))

control.sample <- group$sample[which(group$group=='control')]
## 样本分组
hub_exp2$Group<-ifelse(hub_exp2$sample%in%control.sample,'control','COVID 19')

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

df<-diff[stat.test$Symbol,]
stat.test$p.adj<-df$padj%>%as.numeric()%>%round(digits = 3)
stat.test$p.adj<-ifelse(stat.test$p.adj<0.001,"***",
                        ifelse(stat.test$p.adj<0.05,"**",
                               ifelse(stat.test$p.adj<0.05,"*",'ns')))
hub_exp2$Group <- factor(hub_exp2$Group,levels = c('COVID 19','control'))
exp_plot <- ggplot(hub_exp2,aes(x = Group, y = expr, color = Group)) +
  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.4,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#CD3700","#4682B4"), name = "Group")+
  labs(title="Expression(train)", x="", y = "",size=20) +
  stat_pvalue_manual(stat.test,
                     y.position = c(20,19,18,19,21,19,19,21,20),
                     size = 3.2,
                     family = "Times",
                     label = "p.adj",
                     #parse = T,
                     face = "bold")+
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=15),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=12), 
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
ggsave('01.expression(train).pdf',exp_plot,width = 8,height = 7)
ggsave('01.expression(train).png',exp_plot,width = 8,height = 7)


## 训练集
dat.va<-read.delim2('../00_rawdata/dat(GSE152418).xls',row.names = 1)%>%lc.tableToNum()
dat.va <- na.omit(dat.va)
#colnames(dat.va) <- gsub('X','',colnames(dat.va),fixed = T)
group.va <- read.delim2('../00_rawdata/group(GSE152418).xls')
table(group.va$group)
control.sample <- group.va$sample[which(group.va$group=='control')]
#1 validation-----
hub_exp<-dat.va[hubgene$symbol,]
# hub_exp <- edgeR::cpm(hub_exp)%>%as.data.frame()
hub_exp2<-log2(hub_exp+1)
# hub_exp2 <- hub_exp
hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))

## 样本分组
hub_exp2$Group<-ifelse(hub_exp2$sample%in%control.sample,'control','COVID 19')

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

write.table(stat.test,file = 'wilcox.validation.xls',sep = '\t',row.names = F,quote = F)
stat.test$p.adj <- round(stat.test$p.adj,digits = 3)
stat.test$p.adj<-ifelse(stat.test$p.adj<0.001,"***",
                        ifelse(stat.test$p.adj<0.05,"**",
                               ifelse(stat.test$p.adj<0.05,"*",'ns')))
hub_exp2$Group <- factor(hub_exp2$Group,levels = c('COVID 19','control'))
exp_plot <- ggplot(hub_exp2,aes(x = Group, y = expr, color = Group)) +
  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.4,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#CD3700","#4682B4"), name = "Group")+
  labs(title="Expression(validation)", x="", y = "",size=20) +
  stat_pvalue_manual(stat.test,
                     y.position = c(20,19,18,19,21,19,19,21,20),
                     size = 3.2,
                     family = "Times",
                     label = "p.adj",
                     #parse = T,
                     face = "bold")+
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=15),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=12), 
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
ggsave('02.expression(validation).pdf',exp_plot,width = 8,height = 7)
ggsave('02.expression(validation).png',exp_plot,width = 8,height = 7)

