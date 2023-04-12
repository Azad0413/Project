rm(list = ls())
setwd("/data/nas1/luchunlin/project/CD-0601-2/")
if (! dir.exists("./12_expression")){
  dir.create("./12_expression")
}
setwd("./12_expression")
library(tidyverse)
library(lance)
dat<-read.delim2("../00_rawdata/dat(GSE65682).xls", row.names = 1)%>% lc.tableToNum
group <- read.delim2('../00_rawdata/group(GSE65682).xls')
table(group$group)
hubgene <- read.delim2('../07_features/features.xls')
diff <- read.delim2('../02_DEGs/DEG_sig(GSE65682).xls',row.names = 1)
#1 train-----
hub_exp<-dat[hubgene$symbol,]
hub_exp2<-hub_exp
hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))

control.sample <- group$sample[which(group$group=='control')]
## 样本分组
hub_exp2$Group<-ifelse(hub_exp2$sample%in%control.sample,'control','Sepsis')
hub_exp2$Group <- factor(hub_exp2$Group,levels = c('Sepsis','control'))

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
stat.test$p.adj<-df$adj.P.Val%>%as.numeric()%>%round(digits = 3)
stat.test$p.adj<-ifelse(stat.test$p.adj<0.001,"***",
                        ifelse(stat.test$p.adj<0.05,"**",
                               ifelse(stat.test$p.adj<0.05,"*",'ns')))
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
  labs(title="Expression(Train)", x="", y = "",size=20) +
  stat_pvalue_manual(stat.test,
                     y.position = c(7,10,9,8,8,12,8),
                     size = 3.2,
                     family = "Times",
                     label = "p.adj",
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
ggsave('01.expression(train).pdf',exp_plot,width = 7,height = 6)
ggsave('01.expression(train).png',exp_plot,width = 7,height = 6)


## 训练集
dat.va<-read.delim2('../00_rawdata/dat(GSE28750).xls',row.names = 1)%>%lc.tableToNum()
#colnames(dat.va) <- gsub('X','',colnames(dat.va),fixed = T)
group.va <- read.delim2('../00_rawdata/group(GSE28750).xls')
table(group.va$group)
control.sample <- group.va$sample[which(group.va$group=='control')]

hub_exp<-dat.va[hubgene$symbol,]
hub_exp2<-hub_exp
hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))

## 样本分组
table(group$group)
hub_exp2$Group<-ifelse(hub_exp2$sample%in%control.sample,'control','Sepsis')
hub_exp2$Group <- factor(hub_exp2$Group,levels = c('Sepsis','control'))


##分面图形
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(rstatix)
stat.test<-hub_exp2%>%
  group_by(Symbol)%>%
  t_test(expr ~ Group)%>%
  adjust_pvalue(method = 'fdr')
write.table(stat.test,file = 'GSE28750.test.xls',sep = '\t',row.names = F,quote = F)
stat.test$p<-ifelse(stat.test$p<0.001,"***",
                    ifelse(stat.test$p<0.05,"**",
                           ifelse(stat.test$p<0.05,"*",'ns')))
exp_plot <- ggplot(hub_exp2,aes(x = Group, y = expr, color = Group)) +
  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#4682B4","#CD3700"), name = "Group")+
  labs(title="Expression(Validation)", x="", y = "",size=20) +
  stat_pvalue_manual(stat.test,
                     y.position = c(6.5,10,9,10,8,12,11),
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
ggsave('01.expression(GSE28750).pdf',exp_plot,width = 7,height = 6)
ggsave('01.expression(GSE28750).png',exp_plot,width = 7,height = 6)

