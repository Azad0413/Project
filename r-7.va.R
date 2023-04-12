rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/BJTC-345")
if (! dir.exists("./07_va")){
  dir.create("./07_va")
}
setwd("./07_va")

library(GEOquery)
library(Biobase)
library(tidyverse)
library(dplyr)
gset<-getGEO("GSE64634",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
a=gset[[1]]
pd<-pData(a)
gpl<-getGEO("GPL570",destdir = '.')
gpl<-Table(gpl)    
colnames(gpl)
probe2symbol<-gpl %>%
  dplyr::select('ID','Gene Symbol')%>%
  filter('Gene Symbol'!='')%>%
  separate('Gene Symbol',c('symbol','drop'),sep = '///')%>%
  dplyr::select(-drop)
probe2symbol=probe2symbol[probe2symbol$symbol!='',]
dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symbol$ID<-as.character(probe2symbol$ID)
dat<-dat %>%
  inner_join(probe2symbol,by='ID')%>% 
  dplyr::select(-ID)%>%     ## 去除多余信息
  dplyr::select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
group<-data.frame(sample=pd$geo_accession,group=pd$`disease state:ch1`)
table(group$group)
group$group<-ifelse(group$group=='normal','control','NPC')
table(group$group)
write.table(group,file = 'group.va.xls',sep = '\t',row.names = F,quote = F)
write.table(dat,file = 'dat_va.xls',sep = '\t',row.names = T,quote = F)
hubgene<-read.delim2('/data/nas1/luchunlin/project/BJTC-345/05_expression/hubgene.xls')
group$group = factor(group$group, levels = c("NPC", "control"))
hub_exp<-dat[hubgene$symbol,]
library(ROCR)
library(ggplot2)
hub_exp2<-t(hub_exp)

## 绘制ROC曲线
library(pROC)
for (i in c(1:6)) {
  roc<-roc(group$group,hub_exp2[,i],levels=c("NPC", "control"))
  png(paste0(i, ".", colnames(hub_exp2)[i],".png"),width = 300,height = 300)
  plot(roc,
       print.auc=T,
       print.auc.x=0.4,print.auc.y=0.5,
       #auc.polygon=T,
       #auc.polygon.con="#fff7f7",
       grid=c(0.5,0.2),
       grid.col=c("black","black"),
       #print.thres=T,
       main=colnames(hub_exp2)[i],
       col="#FF2E63",
       legacy.axes=T)
  dev.off()
  pdf(paste0(i, ".", colnames(hub_exp2)[i],".pdf"),width = 4,height = 4)
  plot(roc,
       print.auc=T,
       print.auc.x=0.4,print.auc.y=0.5,
       #auc.polygon=T,
       #auc.polygon.con="#fff7f7",
       grid=c(0.5,0.2),
       grid.col=c("black","black"),
       #print.thres=T,
       main=colnames(hub_exp2)[i],
       col="#FF2E63",
       legacy.axes=T)
  dev.off()
  i<-i+1
}
### -exp-----------
hub_exp<-dat[hubgene$symbol,]
hub_exp2<-hub_exp
hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))
control.sample<-group$sample[which(group$group=='control')]
## 样本分组
hub_exp2$Group<-ifelse(hub_exp2$sample%in%control.sample,'control','NPC')

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
  labs(title="Expression of hubgene", x="", y = "",size=20) +
  stat_pvalue_manual(stat.test,
                     y.position = c(13,9,9,10,9,9),
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
        panel.grid.minor = element_blank())+guides(fill='none')+facet_wrap(~Symbol,scales = "free",nrow = 2) 
exp_plot
ggsave('expression.pdf',exp_plot,width = 6,height = 5)
ggsave('expression.png',exp_plot,width = 6,height = 5)
