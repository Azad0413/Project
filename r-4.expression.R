rm(list = ls())
# 01 获取数据集--------------
setwd("/data/nas1/luchunlin/project/GY0315/")
if (! dir.exists("./04_expression")){
  dir.create("./04_expression")
}
setwd("./04_expression")
library(tidyverse)
library(lance)
dat<-read.delim2("../00_rawdata/dat.xls", row.names = 1)%>% lc.tableToNum
group = read.delim2("../00_rawdata/group.xls")
table(group$group)
group$group<-gsub(' ','_',group$group,fixed = T)
control.sample<-group$sample[which(group$group=='Low_BMD')]
df<-read.delim2('../01_DEGs/DEG_sig.xls',row.names = 1)
hubgene<-read.delim2('../07_validation(roc)/hub.final.xls')
hub_exp<-dat[hubgene$symbol,]
hub_exp2<-hub_exp
hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))

## 样本分组
hub_exp2$Group<-ifelse(hub_exp2$sample%in%control.sample,'Low BMD','High BMD')

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
df<-df[stat.test$Symbol,]
stat.test$p.adj<-df$adj.P.Val%>%as.numeric()%>%round(digits = 3)
stat.test$p.adj<-ifelse(stat.test$p.adj<0.001,"***",
                        ifelse(stat.test$p.adj<0.05,"**",
                               ifelse(stat.test$p.adj<0.05,"*",'ns')))
exp_plot <- ggplot(hub_exp2,aes(x = Group, y = expr, fill = Group)) +
  geom_violin(trim=F) + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA,
               fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#CD3700","#4682B4"), name = "Group")+
  labs(title="Expression level for Hubgene", x="", y = "",size=20) +
  # stat_pvalue_manual(stat.test,
  #                    y.position = c(9,11,9,11,11,9,9,13,13,13),
  #                    size = 3.2,
  #                    family = "Times",
  #                    label = "p.adj",
  #                    #parse = T,
  #                    face = "bold")+
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=15),
        axis.text.x=element_text(angle=0,hjust=0.5,colour="black",face="bold",size=8), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill='none')+facet_wrap(~Symbol,scales = "free",nrow = 3) 
exp_plot
ggsave('expression.pdf',exp_plot,width = 6,height = 6)
ggsave('expression.png',exp_plot,width = 6,height = 6)
