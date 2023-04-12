## 表达分析----------
rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-334")
if (! dir.exists("./05_expression")){
  dir.create("./05_expression")
}
setwd("./05_expression")

library(lance)
dat = read.delim2("/data/nas1/luchunlin/project/BJTC-334/00_rawdata/dat.final.xls", row.names = 1) %>% lc.tableToNum
dat<-log2(dat+1)
group = read.delim2("/data/nas1/luchunlin/project/BJTC-334/00_rawdata/group.xls")
control.sample<-group$sample[which(group$group=='Control')]
hubgene<-data.frame(hubgene=c('JUN','SRC','MAPK3','MMP9','CASP3','EGFR','IL6','TP53','FOS','IL10'))
write.table(hubgene,file = 'hubgene.xls',sep = '\t',row.names = F,quote = F)
hub_exp<-dat[hubgene$hubgene,]
hub_exp2<-hub_exp
hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))

## 样本分组
hub_exp2$Group<-ifelse(hub_exp2$sample%in%control.sample,'control','CAD')


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
write.table(stat.test,file = 'wilcox_results.xls',sep = '\t',row.names = F,quote = F)
##分面图形
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
exp_plot <- ggplot(hub_exp2,aes(x = Group, y = expr, fill = Group)) +
  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar",
               width=0.1,
               position = position_dodge(0.2)) +
  geom_boxplot(width=0.3,
               position=position_dodge(0.9),
               outlier.shape = NA,
               fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#F08080","#20B2AA"), name = "Group")+
  labs(title="", x="", y = "Expression",size=20) +
  stat_compare_means(data = hub_exp2,
                     mapping = aes(group = Group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     label.y = 0.23,
                     label.x = 1.5) +
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
        panel.grid.minor = element_blank())+facet_wrap(~Symbol,scales = "free",nrow = 3) 
exp_plot
ggsave(filename = 'expression.pdf',exp_plot,w=10,h=9)
ggsave(filename = 'expression.png',exp_plot,w=10,h=9)
