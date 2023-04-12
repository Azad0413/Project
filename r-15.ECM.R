rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-321")
if (! dir.exists("./15_ECM")){
  dir.create("./15_ECM")
}
setwd("./15_ECM")
library(tidyverse)
library(lance)
dat<-read.delim2('/data/nas1/luchunlin/project/BJTC-321/00_rawdata/dat.final.xls',row.names = 1)%>%lc.tableToNum()
group<-read.delim2("/data/nas1/luchunlin/project/BJTC-321/00_rawdata/group.xls")
ecm.gene<-data.frame(symbol=c('COL1A1','COL3A1','COL5A1','FN1','ELN','VCAN','TNC','DCN','VTN','POSTN','VIM'))
table(group$group)
control.sample<-group$sample[which(group$group=='control')]
## 
hub_exp<-dat[ecm.gene$symbol,]
hub_exp2<-hub_exp
hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))

## 样本分组
hub_exp2$Group<-ifelse(hub_exp2$sample%in%control.sample,'control','asthma')

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

violin_plot <- ggplot(hub_exp2,aes(x = Symbol, y = expr, color = Group)) +
  #  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#A73030FF", "#0073C2FF"), name = "Group")+
  labs(title="Expression of ECM genes", x="", y = "log2(expr+1)",size=20) +
  stat_compare_means(data = hub_exp2,
                     mapping = aes(group = Group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F) +
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
violin_plot
ggsave('expression.pdf',violin_plot,width = 8,height = 6)
ggsave('expression.png',violin_plot,width = 8,height = 6)

### 相关性--------
hubgene<-read.delim2('/data/nas1/luchunlin/project/BJTC-321/08_va/hub_final.xls')
hub.dat<-dat[hubgene$hubgene,]

library(ggcorrplot)
library(corrplot)
## 提取诊断基因表达矩阵
exp<-rbind(hub.dat,hub_exp)
hub_corr<-round(cor(t(exp)),3)
## 检验基因之间的相关性p值
## 计算相关性系数并显示基因之间的相关性。相关性系数大于0为正相关，小于0为负相关。
## p小于0.05认为相关性显著
hub_p.mat<-round(cor_pmat(t(exp)),3)
col1 <- colorRampPalette(c("#4169E1","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF",
                           "#FDDBC7","#F4A582","#D6604D", "#B2182B","#CC0000","#990000"))
pdf('correlatin.pdf',w=8,h=8)
hub_corr_plot<-corrplot(hub_corr,
                        method = "circle",
                        is.corr = T,
                        type = "lower",
                        p.mat = hub_p.mat,
                        insig = "blank",
                        outline = "white",
                        addCoef.col ="black",
                        tl.col = 'black',
                        col = col1(200),
                        tl.offset = 0.5,
                        number.font = 2,
                        number.cex = 0.8)
dev.off()
png('correlatin.png',w=700,h=700)
hub_corr_plot<-corrplot(hub_corr,
                        method = "circle",
                        is.corr = T,
                        type = "lower",
                        p.mat = hub_p.mat,
                        insig = "blank",
                        outline = "white",
                        addCoef.col ="black",
                        tl.col = 'black',
                        col = col1(200),
                        tl.offset = 0.5,
                        number.font = 2,
                        number.cex = 0.8)
dev.off()
write.table(hub_corr,file = 'hub_corr.xls',
            sep = '\t',
            row.names = T)
