rm(list = ls())
setwd("/data/nas1/luchunlin/project/YQ444-8/")
if (! dir.exists("./06_DEG(TRP)")){
  dir.create("./06_DEG(TRP)")
}
setwd("./06_DEG(TRP)")

coxgene<-read.delim2('/data/nas1/luchunlin/project/YQ444-8/02_uncox/univariate_cox_result_0.05.xls',row.names = 1)
coxgene<-data.frame(symbol=rownames(coxgene))
dat<-read.delim2('/data/nas1/luchunlin/project/YQ444-8/00_rawdata/dat.fpkm.xls',row.names = 1)%>%lc.tableToNum()
colnames(dat)<-gsub('.','-',colnames(dat),fixed = T) 
dat<-dat[,-c(1:8)]
group<-read.delim2("/data/nas1/luchunlin/project/YQ444-8/03_cluster/cluster.xls")%>%as.data.frame()
group<-group[,c('sample','group')]
table(group$group)
control.sample<-group$sample[which(group$group=='Cluster2')]
## 
exp<-dat[coxgene$symbol,]
exp2<-log2(exp+1)
exp2$Symbol<-rownames(exp2)
exp2<-gather(exp2,
             key = sample,
             value = expr,
             -c('Symbol'))

## 样本分组
exp2$Group<-ifelse(exp2$sample%in%control.sample,'Cluster2','Cluster1')
##分面图形
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(rstatix)
stat.test<-exp2%>%
  group_by(Symbol)%>%
  wilcox_test(expr ~ Group)%>%
  adjust_pvalue(method = 'fdr')
write.table(stat.test,file = 'wilcox.test.xls',sep = '\t',row.names = F,quote = F)
exp_plot <- ggplot(exp2, aes(x=Symbol, 
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
  labs(title="", x="", y = "expression level(log2+1)",size=20) +
  stat_compare_means(data = exp2,
                     mapping = aes(group = Group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F) +
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=15),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=12), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())#facet_wrap(~Symbol,scales = "free",nrow = 4) 
exp_plot
ggsave('expression.pdf',exp_plot,width = 8,height = 5)
ggsave('expression.png',exp_plot,width = 8,height = 5)

DETRPs<-stat.test[which(stat.test$p<0.05),]%>%select('Symbol')
write.table(DETRPs,file = 'DETRPs.xls',sep = '\t',row.names = F,quote = F)
