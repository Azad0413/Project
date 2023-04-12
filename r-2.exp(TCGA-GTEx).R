rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/BJTC-299_modIfied/")
if (! dir.exists("./02_exp(TCGA-GTEx)")){
  dir.create("./02_exp(TCGA-GTEx)")
}
setwd("./02_exp(TCGA-GTEx)")

dat<-read.csv('../00_rawdata/dat.merge.xls',sep = '\t')
dat <- edgeR::cpm(dat)
group <- read.delim2('../00_rawdata/group.merge.xls')
table(group$group)
hubgene <- data.frame(symbol=c('IGF1','NGF','GCLM','PYCR1','EFEMP1','APOC3','IFNB1'))
control.sample <- group$sample[which(group$group=='normal')]
#1 train-----
# hub_exp <- dat[rownames(dat)%in%hubgene$symbol,]
hub_exp<-dat[hubgene$symbol,]
hub_exp2<-log2(hub_exp+1)%>%as.data.frame()
hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))

## 样本分组
## 样本分组
hub_exp2$Group<-ifelse(hub_exp2$sample%in%control.sample,'Normal','Tumor')

hub_exp2$Group <- factor(hub_exp2$Group,levels = c('Tumor','Normal'))
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
write.table(stat.test,file = 'wilcox.test.xls',sep = '\t',row.names = F,quote = F)
exp_plot <- ggplot(hub_exp2,aes(x = Group, y = expr, color = Group)) +
  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.3,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#CD3700","#4682B4"), name = "Group")+
  labs(title="Expression(TCGA-GTEx)", x="", y = "",size=20) +
  stat_pvalue_manual(stat.test,
                     y.position = c(8,13,12,2.5,9,5.5,12),
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
ggsave('01.expression(TCGA-GTEx).pdf',exp_plot,width = 8,height = 6)
ggsave('01.expression(TCGA-GTEx).png',exp_plot,width = 8,height = 6)
