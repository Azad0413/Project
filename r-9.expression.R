rm(list = ls())
setwd("/data/nas1/luchunlin/project/NN-0118-2/")
if (! dir.exists("./09_expression")){
  dir.create("./09_expression")
}
setwd("./09_expression")
library(tidyverse)
library(lance)
dat<-read.delim2('../00_rawdata/dat.fpkm.xls',row.names = 1)%>%lc.tableToNum()
colnames(dat)<-gsub('.','-',colnames(dat),fixed = T)
group <- read.delim2('../01_DEGs(TCGA)/group.xls')
table(group$group)
hubgene <- 'CCNB1'
diff <- read.delim2('../01_DEGs(TCGA)/DEG_sig.xls',row.names = 1)
#1 train-----
hub_exp<-dat[hubgene,]
hub_exp2<-log2(hub_exp+1)
hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))

control.sample <- group$sample[which(group$group=='Normal')]
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

df<-diff[stat.test$Symbol,]
stat.test$p.adj<-df$padj%>%as.numeric()%>%round(digits = 3)
stat.test$p.adj<-ifelse(stat.test$p.adj<0.001,"***",
                        ifelse(stat.test$p.adj<0.05,"**",
                               ifelse(stat.test$p.adj<0.05,"*",'ns')))
exp_plot <- ggplot(hub_exp2,aes(x = Group, y = expr, fill = Group)) +
  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.2,
               position=position_dodge(0.9),
               outlier.shape = NA,fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#CD3700","#4682B4"), name = "Group")+
  labs(title="CCNB1 Expression(Train)", x="", y = "Expression Level",size=20) +
  # stat_pvalue_manual(stat.test,
  #                    y.position = c(10),
  #                    size = 3.2,
  #                    family = "Times",
  #                    label = "p.adj",
  #                    #parse = T,
  #                    face = "bold")+
  stat_compare_means(data = hub_exp2,
                     mapping = aes(group = Group),
                     label ="p.signif",
                     method = 'wilcox.test',label.x = 1.5) +
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
   guides(fill='none')#+
  # facet_wrap(~Symbol,scales = "free",nrow = 2) 
exp_plot
ggsave('01.expression(train).pdf',exp_plot,width = 5,height = 5)
ggsave('01.expression(train).png',exp_plot,width = 5,height = 5)

## 验证集------
dat.va<-read.delim2('../00_rawdata/dat(GSE19188).xls',row.names = 1)%>%lc.tableToNum()
dat.va <- na.omit(dat.va)
#colnames(dat.va) <- gsub('X','',colnames(dat.va),fixed = T)
group.va <- read.delim2('../00_rawdata/group(GSE19188).xls')
table(group.va$group)
control.sample <- group.va$sample[which(group.va$group=='Normal')]

hub_exp<-dat.va[hubgene,]
hub_exp2<-hub_exp
hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))

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
  t_test(expr ~ Group)%>%
  adjust_pvalue(method = 'fdr')
write.table(stat.test,file = 'GSE19188.test.xls',sep = '\t',row.names = F,quote = F)
stat.test$p<-ifelse(stat.test$p<0.001,"***",
                    ifelse(stat.test$p<0.05,"**",
                           ifelse(stat.test$p<0.05,"*",'ns')))
exp_plot <- ggplot(hub_exp2,aes(x = Group, y = expr, fill = Group)) +
  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.2,
               position=position_dodge(0.9),
               outlier.shape = NA,fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#CD3700","#4682B4"), name = "Group")+
  labs(title="CCNB1 Expression(Test)", x="", y = "Expression Level",size=20) +
  # stat_pvalue_manual(stat.test,
  #                    y.position = c(10),
  #                    size = 3.2,
  #                    family = "Times",
  #                    label = "p.adj",
  #                    #parse = T,
  #                    face = "bold")+
  stat_compare_means(data = hub_exp2,
                     mapping = aes(group = Group),
                     label ="p.signif",
                     method = 'wilcox.test',label.x = 1.5) +
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
  guides(fill='none')#+
# facet_wrap(~Symbol,scales = "free",nrow = 2) 
exp_plot
ggsave('02.expression(GSE19188).pdf',exp_plot,width = 5,height = 5)
ggsave('02.expression(GSE19188).png',exp_plot,width = 5,height = 5)

