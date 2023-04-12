rm(list = ls())
setwd("/data/nas1/luchunlin/project/SJZZK-431-10/")
if (! dir.exists("./15_checkpoint")){
  dir.create("./15_checkpoint")
}
setwd("./15_checkpoint")

dat<-read.delim2('../00_rawdata/dat.xls',row.names = 1)%>%lc.tableToNum()

risk<-read.delim2('../08_risk/risk.xls')
high.sample<-risk$id[which(risk$risk==0)]
low.sample<-risk$id[which(risk$risk==1)]
risk2 <- risk
risk2$risk_label <- ifelse(risk$risk == 0, "High", "Low")
dat<-dat[,risk$id]
dim(dat)
##   25172    75
## 把检验点基因提取出来
checkpoint <- data.frame(V1=c('PDCD1','CD274','CTLA4','LGALS9','HAVCR2','PDCD1LG2','LAG3','TIGIT'))
checkpoint <- checkpoint$V1
length(checkpoint)
checkpoint_dat <- dat[checkpoint,]
# checkpoint_dat <- dat[which(rownames(dat)%in%checkpoint),]
dim(checkpoint_dat)
#checkpoint_dat<-log2(checkpoint_dat+1)
checkpoint_dat$gene<-rownames(checkpoint_dat)

library(tidyr)
library(ggplot2)
library(ggpubr)
library(Ipaper)
head(checkpoint_dat[,1:3])
violin_dat <- gather(checkpoint_dat, key=sample, value='expr', -c("gene"))
head(violin_dat)
violin_dat$group <- ifelse(violin_dat$sample %in% high.sample,
                           "High", "Low") 
head(violin_dat)
colnames(violin_dat)
library(rstatix)
stat.test<-violin_dat%>%
  group_by(gene)%>%
  wilcox_test(expr ~ group)%>%
  adjust_pvalue(method = 'fdr')
write.table(stat.test,file = 'checkpoint_result.xls',
            sep = '\t',
            row.names = F)

DE.checkpoint <- subset(stat.test,stat.test$p < 0.05)
##2
# violin_dat<-violin_dat[violin_dat$gene%in%DE.checkpoint$gene,]

violin_plot <- ggplot(violin_dat, aes(x=gene, 
                                      y=expr,
                                      fill=group)) +
  #  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#A73030FF", "#0073C2FF"), name = "Group")+
  labs(title="Immune Checkpoint", x="", y = "expression level",size=20) +
  stat_compare_means(data = violin_dat,
                     mapping = aes(group = group),
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
        panel.grid.minor = element_blank())+
  facet_wrap(~gene,scales = "free",nrow = 2) 
violin_plot
ggsave('01.checkpoint.pdf',violin_plot,w=8,h=7)
ggsave('01.checkpoint.png',violin_plot,w=8,h=7)

##相关性分析------
cor_plot_dat <- data.frame(sample=risk$id,riskScore=risk$riskScore,risk_group=risk$risk)
DE.checkpoint_dat <- t(dat[which(rownames(dat)%in%DE.checkpoint$gene),])%>%as.data.frame()
dim(checkpoint_dat)
#DE.checkpoint_dat<-log2(DE.checkpoint_dat+1)
DE.checkpoint_dat$sample=rownames(DE.checkpoint_dat)
cor_plot_dat <- merge(DE.checkpoint_dat,cor_plot_dat,by='sample')
cor_plot_dat$risk_group <- factor(cor_plot_dat$risk_group,
                                  levels = c(0,1),
                                  labels = c('High Risk','Low Risk'))

dim(cor_plot_dat)
library(ggplot2)
library(ggpubr)

cor_plot_dat$riskScore<-as.numeric(cor_plot_dat$riskScore)
#install.packages('ggside')
#install.packages('ggstatsplot')
colnames(cor_plot_dat)
library(ggstatsplot)
checkpoint_cor1 <- ggscatterstats(data = cor_plot_dat,
                                  x = riskScore,
                                  y = LGALS9,
                                  centrality.para = "mean",
                                  margins = "both",
                                  xfill = "#A73030FF",
                                  yfill = "#0073C2FF",
                                  type = "pearson",
                                  ylab = "LGALS9 Expression Level",
                                  marginal.type = "histogram",
                                  title = "Relationship between Checkpoint and riskScore"
)+theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(size = 12))
checkpoint_cor1
ggsave(filename = "02.LGALS9_riskscore_cor.png", height = 7, width = 8,checkpoint_cor1)
ggsave(filename = "02.LGALS9_riskscore_cor.pdf", height = 7, width = 8,checkpoint_cor1)
library(ggstatsplot)
checkpoint_cor2 <- ggscatterstats(data = cor_plot_dat,
                                  x = riskScore,
                                  y = PDCD1LG2,
                                  centrality.para = "mean",
                                  margins = "both",
                                  xfill = "#A73030FF",
                                  yfill = "#0073C2FF",
                                  type = "pearson",
                                  ylab = "PDCD1LG2 Expression Level",
                                  marginal.type = "histogram",
                                  title = "Relationship between Checkpoint and riskScore"
)+theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(size = 12))
checkpoint_cor2
ggsave(filename = "03.PDCD1LG2_riskscore_cor.png", height = 7, width = 8,checkpoint_cor2)
ggsave(filename = "03.PDCD1LG2_riskscore_cor.pdf", height = 7, width = 8,checkpoint_cor2)

l
library(patchwork)
all_correlation <-checkpoint_cor1+checkpoint_cor2+
  plot_layout(ncol = 2) &
  theme_bw() &
  theme(legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold",size = 5))
all_correlation
ggsave('04.correlation_all.pdf',all_correlation,w=12,h=6)
ggsave('04.correlation_all.png',all_correlation,w=12,h=6)
