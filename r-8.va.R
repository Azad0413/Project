## ROC----------
rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-320")
if (! dir.exists("./08_va")){
  dir.create("./08_va")
}
setwd("./08_va")
hubgene<-read.delim2('/data/nas1/luchunlin/project/BJTC-320/06_expression/hubgene.xls')
#hubgene<-hubgene[-6,]%>%as.data.frame()
final.gene<-hubgene
write.table(final.gene,file = 'finalgene.xls',row.names = F,sep = '\t',quote = F)
dat = read.delim2("/data/nas1/luchunlin/project/BJTC-320/00_rawdata/dat_va.xls", row.names = 1) %>% lc.tableToNum
group = read.delim2("/data/nas1/luchunlin/project/BJTC-320/00_rawdata/group_va.xls")
control.sample<-group$sample[which(group$group=='TPH')]

hub_exp<-dat[hubgene$x,]
hub_exp2<-hub_exp
hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))

## 样本分组
hub_exp2$Group<-ifelse(hub_exp2$sample%in%control.sample,'control','HT')

##分面图形
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
exp_plot <- ggplot(hub_exp2,aes(x = Symbol, y = expr, fill = Group)) +
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
                     method = 't.test') +
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
        panel.grid.minor = element_blank())+facet_wrap(~Symbol,scales = "free",nrow = 2) 
exp_plot



group$group = factor(group$group, levels = c("TPH", "HT"))
library(pROC)
library(ggplot2)
hub_exp2<-t(hub_exp)
hub_exp2<-cbind(group,hub_exp2)
## 绘制ROC曲线

roc_FOS<-roc(hub_exp2$group,hub_exp2$FOS,
             levels=c("TPH","HT"))
roc_TNFAIP3<-roc(hub_exp2$group,hub_exp2$TNFAIP3,
                 levels=c("TPH","HT"))
roc_PTK2B<-roc(hub_exp2$group,hub_exp2$PTK2B,
               levels=c("TPH","HT"))
roc_STAT1<-roc(hub_exp2$group,hub_exp2$STAT1,
               levels=c("TPH","HT"))
roc_MMP9<-roc(hub_exp2$group,hub_exp2$MMP9,
              levels=c("TPH","HT"))

plot(roc_FOS,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="FOS ROC curve",
     col="#FF2E63",
     legacy.axes=T)
plot(roc_TNFAIP3,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="TNFAIP3 ROC curve",
     col="#FF2E63",
     legacy.axes=T)
plot(roc_PTK2B,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="PTK2B ROC curve",
     col="#FF2E63",
     legacy.axes=T)
plot(roc_STAT1,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="STAT1 ROC curve",
     col="#FF2E63",
     legacy.axes=T)
plot(roc_MMP9,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="MMP9 ROC curve",
     col="#FF2E63",
     legacy.axes=T)
