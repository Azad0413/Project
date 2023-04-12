
rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-327/")
if (! dir.exists("./04_expression")){
  dir.create("./04_expression")
}
setwd("./04_expression")
library(lance)
library(tidyverse)
# AIH.PBC------
##将基因的矩阵提取出来
hubgene<-read.delim2('/data/nas1/luchunlin/project/BJTC-327/02_DEIRG/DEIRG.xls')
dat1<-read.delim2('/data/nas1/luchunlin/project/BJTC-327/00_rawdata/AIH.PBC.fpkm.xls',row.names = 1)%>%lc.tableToNum()
df1<-read.delim2('/data/nas1/luchunlin/project/BJTC-327/01_DEGs/DEG_sig(AIH vs.PBC).xls',row.names = 1)
dat1<-na.omit(dat1)
dat1<-log2(dat1+1)
group1<-read.delim2("/data/nas1/luchunlin/project/BJTC-327/01_DEGs/group1.xls")
table(group1$group)
control.sample<-group1$sample[which(group1$group=='AIH')]
hub_exp1<-dat1[hubgene$.,]
hub_exp2<-hub_exp1
hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))

## 样本分组
hub_exp2$Group<-ifelse(hub_exp2$sample%in%control.sample,'AIH','PBC')
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
df1<-df1[stat.test$Symbol,]
stat.test$p<-df1$pvalue%>%as.numeric()%>%round(digits = 3)
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
  labs(title="", x="", y = "",size=20) +
  stat_pvalue_manual(stat.test,
                     y.position = c(9,4,1,0.5,7,7,5,4,7,7,6,7.5,1.5,1.5,4,5,2,1.5),
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
        panel.grid.minor = element_blank())+facet_wrap(~Symbol,scales = "free",nrow = 3) 
exp_plot
ggsave('expression(AIH vs.PBC).pdf',exp_plot,width = 10,height = 6)
ggsave('expression(AIH vs.PBC).png',exp_plot,width = 10,height = 6)

### HIH HBV----
##将基因的矩阵提取出来
dat2<-read.delim2('/data/nas1/luchunlin/project/BJTC-327/00_rawdata/AIH.HBV.fpkm.xls',row.names = 1)%>%lc.tableToNum()
df2<-read.delim2('/data/nas1/luchunlin/project/BJTC-327/01_DEGs/DEG_sig(AIH vs.HBV).xls',row.names = 1)
dat2<-na.omit(dat2)
dat2<-log2(dat2+1)
group2<-read.delim2("/data/nas1/luchunlin/project/BJTC-327/01_DEGs/group2.xls")
table(group2$group)
control.sample<-group2$sample[which(group2$group=='AIH')]
hub_exp3<-dat2[hubgene$.,]
hub_exp4<-hub_exp3
hub_exp4$Symbol<-rownames(hub_exp4)
hub_exp4<-gather(hub_exp4,
                 key = sample,
                 value = expr,
                 -c("Symbol"))

## 样本分组
hub_exp4$Group<-ifelse(hub_exp4$sample%in%control.sample,'AIH','HBV')
##分面图形
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(rstatix)
stat.test2<-hub_exp4%>%
  group_by(Symbol)%>%
  t_test(expr ~ Group)%>%
  adjust_pvalue(method = 'fdr')
df2<-df2[stat.test2$Symbol,]
stat.test2$p<-df2$pvalue%>%as.numeric()%>%round(digits = 3)
stat.test2$p<-ifelse(stat.test2$p<0.001,"***",
                    ifelse(stat.test2$p<0.05,"**",
                           ifelse(stat.test2$p<0.05,"*",'ns')))
exp_plot <- ggplot(hub_exp4,aes(x = Group, y = expr, color = Group)) +
  #geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#4682B4","#CD3700"), name = "Group")+
  labs(title="", x="", y = "",size=20) +
  stat_pvalue_manual(stat.test2,
                     y.position = c(9,4,1,0.5,7,7,5,4,7,7,6,7.5,1.5,1.5,4,5,2,1.5),
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
        panel.grid.minor = element_blank())+facet_wrap(~Symbol,scales = "free",nrow = 3) 
exp_plot
ggsave('expression(AIH vs.HBV).pdf',exp_plot,width = 10,height = 6)
ggsave('expression(AIH vs.HBV).png',exp_plot,width = 10,height = 6)
