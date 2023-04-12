rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-300-8/")
if (! dir.exists("./10_modelgene")){
  dir.create("./10_modelgene")
}
setwd("./10_modelgene")

#expreesion------
dat.tcga<-read.delim2("../00_rawdata/dat.fpkm.xls", row.names = 1)%>% lc.tableToNum
colnames(dat.tcga)<-gsub('.','-',colnames(dat.tcga),fixed = T)
dat1 <- read.delim2('../00_rawdata/dat(GSE16088).xls')%>% lc.tableToNum
dat2 <- read.delim2('../00_rawdata/dat(GSE19276).xls')%>% lc.tableToNum
dat3 <- read.delim2('../00_rawdata/dat(GSE99671).xls')%>% lc.tableToNum
group1 <- read.delim2('../00_rawdata/group(GSE16088).xls')
group2 <- read.delim2('../00_rawdata/group(GSE19276).xls')
group3 <- read.delim2('../00_rawdata/group(GSE99671).xls')
diff1 <- read.delim2('../01_DEGs(GSE16088)/DEG_sig(GSE16088).xls')
diff2 <- read.delim2('../03_DEGs(GSE19276)/DEG_sig(GSE19276).xls')
diff3 <- read.delim2('../02_DEGs(GSE99671)/DEG_sig(GSE99671).xls',row.names = 1)
gene <- read.csv('../07_Multivariate_cox/multivariate_cox_result.csv')$X%>%as.data.frame()

##1 GSE16088------
hub_exp<-dat1[gene$.,]
hub_exp2<-hub_exp
hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))

control.sample <- group1$sample[which(group1$group=='control')]
## 样本分组
hub_exp2$Group<-ifelse(hub_exp2$sample%in%control.sample,'control','OS')

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

df<-diff1[stat.test$Symbol,]
stat.test$p.adj<-df$adj.P.Val%>%as.numeric()%>%round(digits = 3)
stat.test$p.adj<-ifelse(stat.test$p.adj<0.001,"***",
                        ifelse(stat.test$p.adj<0.05,"**",
                               ifelse(stat.test$p.adj<0.05,"*",'ns')))
exp_plot <- ggplot(hub_exp2,aes(x = Group, y = expr, color = Group)) +
  #geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#CD3700","#4682B4"), name = "Group")+
  labs(title="Expression(GSE16088)", x="", y = "",size=20) +
  stat_pvalue_manual(stat.test,
                     y.position = c(12,10),
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
  facet_wrap(~Symbol,scales = "free",nrow = 1) 
exp_plot
ggsave('01.expression(GSE16088).pdf',exp_plot,width = 6,height = 4)
ggsave('01.expression(GSE16088).png',exp_plot,width = 6,height = 4)
##2 GSE19276------
hub_exp<-dat2[gene$.,]
hub_exp2<-hub_exp
hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))

control.sample <- group2$sample[which(group2$group=='control')]
## 样本分组
hub_exp2$Group<-ifelse(hub_exp2$sample%in%control.sample,'control','OS')

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

df<-diff2[stat.test$Symbol,]
stat.test$p.adj<-df$adj.P.Val%>%as.numeric()%>%round(digits = 3)
stat.test$p.adj<-ifelse(stat.test$p.adj<0.001,"***",
                        ifelse(stat.test$p.adj<0.05,"**",
                               ifelse(stat.test$p.adj<0.05,"*",'ns')))
exp_plot <- ggplot(hub_exp2,aes(x = Group, y = expr, color = Group)) +
  #geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#CD3700","#4682B4"), name = "Group")+
  labs(title="Expression(GSE19276)", x="", y = "",size=20) +
  stat_pvalue_manual(stat.test,
                     y.position = c(4,8),
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
  facet_wrap(~Symbol,scales = "free",nrow = 1) 
exp_plot
ggsave('02.expression(GSE19276).pdf',exp_plot,width = 6,height = 4)
ggsave('02.expression(GSE99671).png',exp_plot,width = 6,height = 4)
##2 GSE99671------
hub_exp<-dat3[gene$.,]
hub_exp <- edgeR::cpm(hub_exp)%>%as.data.frame()
hub_exp2<-log2(hub_exp+1)
hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))

control.sample <- group3$sample[which(group3$group=='control')]
## 样本分组
hub_exp2$Group<-ifelse(hub_exp2$sample%in%control.sample,'control','OS')

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

df<-diff3[stat.test$Symbol,]
stat.test$p.adj<-df$padj%>%as.numeric()%>%round(digits = 3)
stat.test$p.adj<-ifelse(stat.test$p.adj<0.001,"***",
                        ifelse(stat.test$p.adj<0.05,"**",
                               ifelse(stat.test$p.adj<0.05,"*",'ns')))
exp_plot <- ggplot(hub_exp2,aes(x = Group, y = expr, color = Group)) +
  #geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#CD3700","#4682B4"), name = "Group")+
  labs(title="Expression(GSE99671)", x="", y = "",size=20) +
  stat_pvalue_manual(stat.test,
                     y.position = c(20,20),
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
  facet_wrap(~Symbol,scales = "free",nrow = 1) 
exp_plot
ggsave('03.expression(GSE99671).pdf',exp_plot,width = 6,height = 4)
ggsave('03.expression(GSE99671).png',exp_plot,width = 6,height = 4)

## KM

survival<-read.delim2('../06_univariate_cox/survival.xls') 
dat<-t(dat.tcga)%>%as.data.frame()
dat <- log2(dat+1)
## KM（BNIP3）-----
group<-data.frame(sample=rownames(dat),group=ifelse(dat$BNIP3>median(dat$BNIP3),'High BNIP3','Low BNIP3'))
write.table(group1,file = 'group(BNIP3).xls',row.names = F,sep = '\t',quote = F)
group$group<-as.vector(group$group)
km.dat<-dat%>%as.data.frame()
km.dat$group<-as.vector(group$group,)
km.dat$sample<-rownames(km.dat)
km.dat<-merge(survival,km.dat,by='sample')
rownames(km.dat)<-km.dat$sample
km.dat<-km.dat[,-1]
table(km.dat$group)
kmfit<-survfit(Surv(OS.time, OS) ~ group, data =  km.dat)
cluster_survival_median1 <- ggsurvplot(kmfit,
                                       pval = TRUE, 
                                       conf.int = F,
                                       legend.labs=c("High BNIP3","Low BNIP3" ),
                                       legend.title="group",
                                       title="KM(BNIP3)",
                                       font.main = c(15,"bold"),
                                       risk.table = TRUE, 
                                       risk.table.col = "strata", 
                                       linetype = "strata", 
                                       surv.median.line = "hv", 
                                       ggtheme = theme_bw(), 
                                       palette = c("#A73030FF", "#0073C2FF"))
cluster_survival_median1


## KM（CXCL12）-----
group<-data.frame(sample=rownames(dat),group=ifelse(dat$CXCL12>median(dat$CXCL12),'High CXCL12','Low CXCL12'))
write.table(group1,file = 'group(CXCL12).xls',row.names = F,sep = '\t',quote = F)
group$group<-as.vector(group$group)
km.dat<-dat%>%as.data.frame()
km.dat$group<-as.vector(group$group,)
km.dat$sample<-rownames(km.dat)
km.dat<-merge(survival,km.dat,by='sample')
rownames(km.dat)<-km.dat$sample
km.dat<-km.dat[,-1]
table(km.dat$group)
kmfit<-survfit(Surv(OS.time, OS) ~ group, data =  km.dat)
cluster_survival_median1 <- ggsurvplot(kmfit,
                                       pval = TRUE, 
                                       conf.int = F,
                                       legend.labs=c("High CXCL12","Low CXCL12" ),
                                       legend.title="group",
                                       title="KM(CXCL12)",
                                       font.main = c(15,"bold"),
                                       risk.table = TRUE, 
                                       risk.table.col = "strata", 
                                       linetype = "strata", 
                                       surv.median.line = "hv", 
                                       ggtheme = theme_bw(), 
                                       palette = c("#A73030FF", "#0073C2FF"))
cluster_survival_median1
