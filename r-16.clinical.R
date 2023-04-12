rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-321")
if (! dir.exists("./16_clinical")){
  dir.create("./16_clinical")
}
setwd("./16_clinical")
hubgene<-read.delim2('/data/nas1/luchunlin/project/BJTC-321/08_va/hub_final.xls')
dat = read.delim2("/data/nas1/luchunlin/project/BJTC-321/00_rawdata/dat.final.xls", row.names = 1) %>% lc.tableToNum
hub.exp<-dat[hubgene$hubgene,]
group<-read.delim2("/data/nas1/luchunlin/project/BJTC-321/00_rawdata/group.xls")
clinical<-read.delim2('/data/nas1/luchunlin/project/BJTC-321/00_rawdata/clinical.xls')
table(clinical$age)
clinical$age<-cut(clinical$age,breaks = c(30,50,70),labels = c('30-50','50-70'))
library(ggplot2)
library(ggpubr)
##age-----
hub_exp2<-hub.exp
hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))

## 样本分组
control.sample<-clinical$sample[which(clinical$age=='30-50')]
hub_exp2$Group<-ifelse(hub_exp2$sample%in%control.sample,'30-50','50-70')

##分面图形
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(rstatix)
exp_plot1 <- ggplot(hub_exp2,aes(x = Symbol, y = expr, fill = Group)) +
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
  stat_compare_means(data = hub_exp2,
                     mapping = aes(group = Group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F) +
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
        panel.grid.minor = element_blank())#+facet_wrap(~Symbol,scales = "free",nrow = 3) 
exp_plot1
ggsave('age.pdf',exp_plot,w=8,h=5)
ggsave('age.png',exp_plot,w=8,h=5)

## geneder-----------

hub_exp2<-hub.exp
hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))

## 样本分组
table(clinical$gender)
control.sample<-clinical$sample[which(clinical$gender=='male')]
hub_exp2$Group<-ifelse(hub_exp2$sample%in%control.sample,'male','female')
##分面图形
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(rstatix)
exp_plot2 <- ggplot(hub_exp2,aes(x = Symbol, y = expr, fill = Group)) +
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
  stat_compare_means(data = hub_exp2,
                     mapping = aes(group = Group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F) +
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
        panel.grid.minor = element_blank())#+facet_wrap(~Symbol,scales = "free",nrow = 3) 
exp_plot2
ggsave('gender.pdf',exp_plot,w=8,h=5)
ggsave('gender.png',exp_plot,w=8,h=5)
library(patchwork)
all_clinical_index <- exp_plot1+exp_plot2+
  plot_layout(ncol = 1) &
  theme_bw() &
  theme(legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold",size = 12))
all_clinical_index
ggsave('clinical.pdf',all_clinical_index,w=8,h=8)
ggsave('clinical.png',all_clinical_index,w=8,h=8)
