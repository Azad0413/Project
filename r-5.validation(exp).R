rm(list = ls())
# 01 获取数据集--------------
setwd("/data/nas1/luchunlin/project/GY0315/")
if (! dir.exists("./05_validation")){
  dir.create("./05_validation")
}
setwd("./05_validation")
library(GEOquery)
library(Biobase)
library(limma)
library(tidyverse)
gset<-getGEO("GSE62402",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
a=gset[[1]]
pd<-pData(a)
gpl<-getGEO("GPL5175",destdir = '.')
gpl<-Table(gpl)    
colnames(gpl)
probe2symbol<-gpl %>%
  dplyr::select('ID','gene_assignment')%>%
  filter('gene_assignment'!='')%>%
  separate('gene_assignment',c('drop','symbol'),sep = '//')%>%
  dplyr::select(-drop)
probe2symbol=probe2symbol[probe2symbol$symbol!='',]
probe2symbol<-na.omit(probe2symbol)
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
## 
rownames(dat)<-gsub(' ','',rownames(dat),fixed = T)
write.table(dat,file = 'dat(GSE62402).xls',sep = '\t',row.names = T,quote = F)
group<-data.frame(sample=pd$geo_accession,group=pd$title)
table(group$group)
group$group<-c(rep('High BMD',5),rep('Low BMD',5))
#group$group<-ifelse(group$group=='bone mineral density: high BMD','High BMD','Low BMD')
write.table(group,file = 'group(GSE62402).xls',sep = '\t',row.names = F,quote = F)
table(group$group)
group$group<-gsub(' ','_',group$group,fixed = T)
control.sample<-group$sample[which(group$group=='Low_BMD')]
hubgene<-read.delim2('../03_DEmodERS/intersect.xls')
hub_exp<-dat[hubgene$.,]
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
# stat.test$p.adj<-df$adj.P.Val%>%as.numeric()%>%round(digits = 3)
# stat.test$p.adj<-ifelse(stat.test$p.adj<0.001,"***",
#                         ifelse(stat.test$p.adj<0.05,"**",
#                                ifelse(stat.test$p.adj<0.05,"*",'ns')))
exp_plot <- ggplot(hub_exp2,aes(x = Group, y = expr, fill = Group)) +
  geom_violin(trim=F) + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA,fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#CD3700","#4682B4"), name = "Group")+
  labs(title="Expression", x="", y = "",size=20) +
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
  guides(fill='none')+
  facet_wrap(~Symbol,scales = "free",nrow = 3) 
exp_plot
ggsave('expression.pdf',exp_plot,width = 8,height = 7)
ggsave('expression.png',exp_plot,width = 8,height = 7)


##----------

# gset<-getGEO("GSE7429",
#              destdir = '.',
#              GSEMatrix = T,
#              getGPL = F)
# expr<-as.data.frame(exprs(gset[[1]]))
# a=gset[[1]]
# pd<-pData(a)
# gpl<-getGEO("GPL96",destdir = '.')
# gpl<-Table(gpl)    
# colnames(gpl)
# probe2symbol<-gpl %>%
#   dplyr::select('ID','Gene Symbol')%>%
#   filter('Gene Symbol'!='')%>%
#   separate('Gene Symbol',c('symbol','drop'),sep = '///')%>%
#   dplyr::select(-drop)
# probe2symbol=probe2symbol[probe2symbol$symbol!='',]
# dat<-expr
# dat$ID<-rownames(dat)
# dat$ID<-as.character(dat$ID)
# probe2symbol$ID<-as.character(probe2symbol$ID)
# dat<-dat %>%
#   inner_join(probe2symbol,by='ID')%>% 
#   dplyr::select(-ID)%>%     ## 去除多余信息
#   dplyr::select(symbol,everything())%>%     ## 重新排列
#   mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
#   arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
#   distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
#   dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
#   tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
# ## 
# write.table(dat,file = 'dat(GSE7429).xls',sep = '\t',row.names = T,quote = F)
# group<-data.frame(sample=pd$geo_accession,group=pd$title)
# group$group<-c(rep('High BMD',10),rep('Low BMD',10))
# write.table(group,file = 'group(GSE7429).xls',sep = '\t',row.names = F,quote = F)
# 
# table(group$group)
# group$group<-gsub(' ','_',group$group,fixed = T)
# control.sample<-group$sample[which(group$group=='Low_BMD')]
# hubgene<-read.delim2('../03_DEmodERS/intersect.xls')
# hub_exp<-dat[hubgene$.,]
# hub_exp2<-hub_exp
# hub_exp2$Symbol<-rownames(hub_exp2)
# hub_exp2<-gather(hub_exp2,
#                  key = sample,
#                  value = expr,
#                  -c("Symbol"))
# 
# ## 样本分组
# hub_exp2$Group<-ifelse(hub_exp2$sample%in%control.sample,'Low BMD','High BMD')
# 
# ##分面图形
# library(ggplot2)
# library(tidyverse)
# library(ggpubr)
# library(ggsci)
# library(rstatix)
# stat.test<-hub_exp2%>%
#   group_by(Symbol)%>%
#   wilcox_test(expr ~ Group)%>%
#   adjust_pvalue(method = 'fdr')
# # stat.test$p.adj<-df$adj.P.Val%>%as.numeric()%>%round(digits = 3)
# # stat.test$p.adj<-ifelse(stat.test$p.adj<0.001,"***",
# #                         ifelse(stat.test$p.adj<0.05,"**",
# #                                ifelse(stat.test$p.adj<0.05,"*",'ns')))
# exp_plot <- ggplot(hub_exp2,aes(x = Group, y = expr, color = Group)) +
#   geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
#   #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
#   stat_boxplot(geom="errorbar", 
#                width=0.1,
#                position = position_dodge(0.9)) +
#   geom_boxplot(width=0.7,
#                position=position_dodge(0.9),
#                outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
#   scale_fill_manual(values= c("#4682B4","#CD3700"), name = "Group")+
#   labs(title="Expression", x="", y = "",size=20) +
#   # stat_pvalue_manual(stat.test,
#   #                    y.position = c(9,11,9,11,11,9,9,13,13,13),
#   #                    size = 3.2,
#   #                    family = "Times",
#   #                    label = "p.adj",
#   #                    #parse = T,
#   #                    face = "bold")+
#   theme_bw()+
#   theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=15),
#         axis.text.x=element_text(angle=0,hjust=0.5,colour="black",face="bold",size=8), 
#         axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
#         axis.title.x=element_text(size=16,face="bold"),
#         axis.title.y=element_text(size=16,face="bold"),
#         legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
#         legend.title = element_text(face = "bold", size = 12),
#         legend.position = "top",
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())+
#   guides(fill='none')+
#   facet_wrap(~Symbol,scales = "free",nrow = 3) 
# exp_plot
# ggsave('expression.pdf',exp_plot,width = 8,height = 7)
# ggsave('expression.png',exp_plot,width = 8,height = 7)
# 
