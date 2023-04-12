rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-386-10/")
if (! dir.exists("./13_clinical(exp.va)")){
  dir.create("./13_clinical(exp.va)")
}
setwd("./13_clinical(exp.va)")

## 01.验证集-------
dat<-read.delim2('../00_rawdata/dat.va.xls')%>%lc.tableToNum()
dat <- na.omit(dat)
dat <- t(dat)%>%as.data.frame()
phenotype<-read.delim2('../00_rawdata/Merge_clinical.txt')
validation_phenotype<-data.frame(sample=phenotype$icgc_sample_id,
                                 age=phenotype$donor_age_at_diagnosis,
                                 gender=phenotype$donor_sex,
                                 T.stage=phenotype$donor_tumour_stage_at_diagnosis,
                                 N.stage=phenotype$donor_tumour_stage_at_diagnosis,
                                 M.stage=phenotype$donor_tumour_stage_at_diagnosis)
dat$sample <- rownames(dat)
test.data<<-merge(validation_phenotype,dat,by='sample')%>%column_to_rownames(var = 'sample')
hub.exp <- t(test.data[,c(6:12)])%>%as.data.frame()%>%lc.tableToNum()
pd <- test.data[,c(1:5)]
pd2 <- pd
table(pd2$age)
pd2$age<-ifelse(pd2$age>60,'>60','<=60')

table(pd2$T.stage)
pd2$T.stage <- substr(pd2$T.stage,1,2)
pd2$T.stage<-gsub('T3','T3/4',pd2$T.stage)
pd2$T.stage<-gsub('T4','T3/4',pd2$T.stage)
pd2$T.stage<-gsub('T1','T1/2',pd2$T.stage)
pd2$T.stage<-gsub('T2','T1/2',pd2$T.stage)

table(pd2$N.stage)
pd2$N.stage <- substr(pd2$N.stage,3,4)
pd2$N.stage<-gsub('NX',NA,pd2$N.stage)

table(pd2$M.stage)
pd2$M.stage <- substr(pd2$M.stage,5,6)
pd2$M.stage<-gsub('MX',NA,pd2$M.stage)

colnames(pd2)
colnames(pd2)<-c('Age','Gender','T.stage','N.stage','M.stage')

hub_exp2<-log2(hub.exp+1)%>%as.data.frame()
hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))
## age-----
## 样本分组
age.dat <- data.frame(sample=rownames(pd2),group=pd2$Age)
age.dat <- na.omit(age.dat)
control.sample <- rownames(pd2)[which(pd2$Age=='>60')]
hub_exp2 <- hub_exp2[hub_exp2$sample%in%age.dat$sample,]
hub_exp2$Group<-ifelse(hub_exp2$sample%in%control.sample,'>60','<=60')
##分面图形
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(rstatix)
library(lance)
stat.test<-hub_exp2%>%lc.tableToNum()%>%
  group_by(Symbol)%>%
  wilcox_test(expr ~ Group)%>%
  adjust_pvalue(method = 'fdr')
write.table(stat.test,file = 'wilcox.age.xls',sep = '\t',row.names = F,quote = F)
hub_exp2<-lc.tableToNum(hub_exp2)
age_plot <- ggplot(hub_exp2, aes(x=Symbol, 
                                 y=expr,
                                 fill=Group)) +
  #  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#A73030FF", "#0073C2FF"), name = "Group")+
  labs(title="Expression between Differett Age(Validation)", x="", y = "expression level",size=20) +
  stat_compare_means(data = hub_exp2,
                     mapping = aes(group = Group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F,
                     label.y = 4.5) +
  #  geom_signif(comparisons = my_comparisons,
  #              test = t.test,
  #              map_signif_level = T)+
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=8),
        legend.title = element_text(face = "bold", size = 10),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
age_plot
ggsave('01.Age(Validation).pdf',age_plot,width = 7,height = 5)
ggsave('01.Age(Validation).png',age_plot,width = 7,height = 5)


## Gender-----
## 样本分组
Gender.dat <- data.frame(sample=rownames(pd2),group=pd2$Gender)
Gender.dat <- na.omit(Gender.dat)
table(pd2$Gender)
control.sample <- rownames(pd2)[which(pd2$Gender=='male')]
hub_exp2 <- hub_exp2[hub_exp2$sample%in%Gender.dat$sample,]
hub_exp2$Group<-ifelse(hub_exp2$sample%in%control.sample,'male','female')
##分面图形
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(rstatix)
library(lance)
stat.test<-hub_exp2%>%lc.tableToNum()%>%
  group_by(Symbol)%>%
  wilcox_test(expr ~ Group)%>%
  adjust_pvalue(method = 'fdr')
write.table(stat.test,file = 'wilcox.Gender.xls',sep = '\t',row.names = F,quote = F)
hub_exp2<-lc.tableToNum(hub_exp2)
Gender_plot <- ggplot(hub_exp2, aes(x=Symbol, 
                                    y=expr,
                                    fill=Group)) +
  #  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#A73030FF", "#0073C2FF"), name = "Group")+
  labs(title="Expression between Differentt Gender(Validation)", x="", y = "expression level",size=20) +
  stat_compare_means(data = hub_exp2,
                     mapping = aes(group = Group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F,
                     label.y = 4.5) +
  #  geom_signif(comparisons = my_comparisons,
  #              test = t.test,
  #              map_signif_level = T)+
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=8),
        legend.title = element_text(face = "bold", size = 10),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
Gender_plot
ggsave('02.Gender(Validation).pdf',Gender_plot,width = 7,height = 5)
ggsave('02.Gender(Validation).png',Gender_plot,width = 7,height = 5)


## T.stage-----
## 样本分组
T.stage.dat <- data.frame(sample=rownames(pd2),group=pd2$T.stage)
T.stage.dat <- na.omit(T.stage.dat)
table(pd2$T.stage)
control.sample <- rownames(pd2)[which(pd2$T.stage=='T1/2')]
hub_exp2 <- hub_exp2[hub_exp2$sample%in%T.stage.dat$sample,]
hub_exp2$Group<-ifelse(hub_exp2$sample%in%control.sample,'T1/2','T3/4')
##分面图形
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(rstatix)
library(lance)
stat.test<-hub_exp2%>%lc.tableToNum()%>%
  group_by(Symbol)%>%
  wilcox_test(expr ~ Group)%>%
  adjust_pvalue(method = 'fdr')
write.table(stat.test,file = 'wilcox.T.stage.xls',sep = '\t',row.names = F,quote = F)
hub_exp2<-lc.tableToNum(hub_exp2)
T.stage_plot <- ggplot(hub_exp2, aes(x=Symbol, 
                                     y=expr,
                                     fill=Group)) +
  #  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#A73030FF", "#0073C2FF","green"), name = "Group")+
  labs(title="Expression between Different T.stage(Validation)", x="", y = "expression level",size=20) +
  stat_compare_means(data = hub_exp2,
                     mapping = aes(group = Group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F,
                     label.y = 4.5) +
  #  geom_signif(comparisons = my_comparisons,
  #              test = t.test,
  #              map_signif_level = T)+
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=8),
        legend.title = element_text(face = "bold", size = 10),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
T.stage_plot
ggsave('03.T.stage(Validation).pdf',T.stage_plot,width = 7,height = 5)
ggsave('03.T.stage(Validation).png',T.stage_plot,width = 7,height = 5)

## N.stage-----
## 样本分组
N.stage.dat <- data.frame(sample=rownames(pd2),group=pd2$N.stage)
N.stage.dat <- na.omit(N.stage.dat)
table(pd2$N.stage)
control.sample <- rownames(pd2)[which(pd2$N.stage=='N0')]
hub_exp2 <- hub_exp2[hub_exp2$sample%in%N.stage.dat$sample,]
hub_exp2$Group<-ifelse(hub_exp2$sample%in%control.sample,'N0','N1')
##分面图形
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(rstatix)
library(lance)
stat.test<-hub_exp2%>%lc.tableToNum()%>%
  group_by(Symbol)%>%
  wilcox_test(expr ~ Group)%>%
  adjust_pvalue(method = 'fdr')
write.table(stat.test,file = 'wilcox.N.stage.xls',sep = '\t',row.names = F,quote = F)
hub_exp2<-lc.tableToNum(hub_exp2)
N.stage_plot <- ggplot(hub_exp2, aes(x=Symbol, 
                                     y=expr,
                                     fill=Group)) +
  #  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#A73030FF", "#0073C2FF","green"), name = "Group")+
  labs(title="Expression between Different N.stage(Validation)", x="", y = "expression level",size=20) +
  stat_compare_means(data = hub_exp2,
                     mapping = aes(group = Group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F,
                     label.y = 4.5) +
  #  geom_signif(comparisons = my_comparisons,
  #              test = N.test,
  #              map_signif_level = T)+
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=8),
        legend.title = element_text(face = "bold", size = 10),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
N.stage_plot
ggsave('04.N.stage(Validation).pdf',N.stage_plot,width = 7,height = 5)
ggsave('04.N.stage(Validation).png',N.stage_plot,width = 7,height = 5)

## M.stage-----
## 样本分组
M.stage.dat <- data.frame(sample=rownames(pd2),group=pd2$M.stage)
M.stage.dat <- na.omit(M.stage.dat)
table(pd2$M.stage)
control.sample <- rownames(pd2)[which(pd2$M.stage=='M0')]
hub_exp2 <- hub_exp2[hub_exp2$sample%in%M.stage.dat$sample,]
hub_exp2$Group<-ifelse(hub_exp2$sample%in%control.sample,'M0','M1')
##分面图形
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(rstatix)
library(lance)
stat.test<-hub_exp2%>%lc.tableToNum()%>%
  group_by(Symbol)%>%
  wilcox_test(expr ~ Group)%>%
  adjust_pvalue(method = 'fdr')
write.table(stat.test,file = 'wilcox.M.stage.xls',sep = '\t',row.names = F,quote = F)
hub_exp2<-lc.tableToNum(hub_exp2)
M.stage_plot <- ggplot(hub_exp2, aes(x=Symbol, 
                                     y=expr,
                                     fill=Group)) +
  #  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#A73030FF", "#0073C2FF","green"), name = "Group")+
  labs(title="Expression between Different M.stage(Validation)", x="", y = "expression level",size=20) +
  stat_compare_means(data = hub_exp2,
                     mapping = aes(group = Group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F,
                     label.y = 4.5) +
  #  geom_signif(comparisons = my_comparisons,
  #              test = M.test,
  #              map_signif_level = T)+
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=8),
        legend.title = element_text(face = "bold", size = 10),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
M.stage_plot
ggsave('05.M.stage(Validation).pdf',M.stage_plot,width = 7,height = 5)
ggsave('05.M.stage(Validation).png',M.stage_plot,width = 7,height = 5)


library(patchwork)
all <- age_plot+Gender_plot+T.stage_plot+N.stage_plot+M.stage_plot+
  plot_layout(ncol = 2)&
  theme_bw()&
  theme(legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=7), )
all
ggsave('06.all_clinical.pdf',all,w=10,h=10)
ggsave('06.all_clinical.png',all,w=11,h=11)
