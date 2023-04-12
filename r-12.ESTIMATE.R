rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-399-11/")
if (! dir.exists("./12_ESTIMATE")){
  dir.create("./12_ESTIMATE")
}
setwd("./12_ESTIMATE")
dat<-read.delim2("../00_rawdata/dat.fpkm.xls", row.names = 1)%>% lc.tableToNum
dat.pcg <- read.delim2('../00_rawdata/dat.pcg.xls',row.names = 1)%>%lc.tableToNum()
dat <- dat[rownames(dat.pcg),]
dat <- na.omit(dat)
colnames(dat)<-gsub('.','-',colnames(dat),fixed = T)
group<-read.delim2('../08_risk/risk.xls')%>%dplyr::select(c('id','riskScore','risk'))
colnames(group)<-c('sample','riskScore','label')
group$label <- ifelse(group$label=='1','Low risk','High risk')
Low.sample<-group$sample[which(group$label=='Low risk')]
High.sample<-group$sample[which(group$label=='High risk')]
dat<-dat[,group$sample]
group_estimate<-group$label%>%as.factor()
design<-model.matrix(~0 + group_estimate)
rownames(design)<-group$sample
colnames(design)<-levels(group_estimate)
design<-as.data.frame(design)
Low<-rownames(design)[which(design$`Low risk`==1)]
High<-rownames(design)[which(design$`High risk`==1)]
length(Low)
length(High)
#install.packages("estimate", repos="http://R-Forge.R-project.org")
library(estimate)
expr_train <- log2(dat+1)
#expr_train <- dat
write.table(expr_train, 
            'expr.txt', 
            col.names = T, 
            row.names = T, 
            quote = F, sep="\t")
# 生成expr_train.gct
filterCommonGenes(input.f = './expr.txt', 
                  output.f = 'expr_train.gct', 
                  id = 'GeneSymbol')
# [1] "Merged dataset includes 10221 genes (191 mismatched)."

# 生成train_purity.gct
estimateScore('expr_train.gct', 'train_purity.gct', platform="affymetrix")
es_score <- read.table('train_purity.gct', skip = 2, header = T, check.names = F)
immu_score <- es_score[,3:length(es_score)]
rownames(immu_score) <- es_score$NAME
write.table(es_score,
            file = "es_score.xls",
            sep = "\t",
            quote = F,
            row.names = F)
## 小提琴图--------
violin_dat <- data.frame(t(immu_score))
rownames(violin_dat)<-group$sample
violin_dat$sample <- rownames(violin_dat)
violin_dat$group <- ifelse(violin_dat$sample %in% Low.sample,
                           "Low Risk", "High Risk")
library(tidyr)
library(ggplot2)
library(ggpubr)
library(Ipaper)
{
  p1 <- ggplot(violin_dat, aes(x=group,y=StromalScore, fill=group))+
    geom_violin() + #绘制小提琴图
    stat_boxplot(geom="errorbar",
                 width=0.1,
                 position = position_dodge(0.9)) +
    geom_boxplot(width=0.4,
                 position=position_dodge(0.9),
                 outlier.shape = NA,
                 fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
    geom_point(aes(fill = group),
               size = 0.5,
               position = position_dodge(0.9))+
    scale_fill_manual(values = c("#FF6347","#6495ED"), name = "Group") + 
    stat_compare_means(aes(group = group),method = 'wilcox')+ 
    theme_bw()+ #背景变为白色
    ylim(-2000,2000) +
    theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
          panel.grid = element_blank()) + 
    theme(axis.text.x=element_text(family="Times",size=12,face="bold"), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
          axis.text.y=element_text(family="Times",size=12,face="bold"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
          axis.title.y=element_text(family="Times",size = 15,face="bold"), #设置y轴标题的字体属性
          axis.title.x=element_text(family="Times",size = 10,face="bold"),
          plot.title = element_text(hjust = 0.5,family="Times",size = 15,face="bold"),
          legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold"))+
    labs(title="Stromal Score", x="", y="Stromal Score")+
    guides(fill='none')
  p1
  p2 <- ggplot(violin_dat, aes(x=group,y=ImmuneScore, fill=group))+
    geom_violin() + #绘制小提琴图
    stat_boxplot(geom="errorbar",
                 width=0.1,
                 position = position_dodge(0.9)) +
    geom_boxplot(width=0.4,
                 position=position_dodge(0.9),
                 outlier.shape = NA,
                 fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
    geom_point(aes(fill = group),
               size = 0.5,
               position = position_dodge(0.9))+
    scale_fill_manual(values = c("#FF6347","#6495ED"), name = "Group") +
    stat_compare_means(aes(group = group),method = 'wilcox')+ 
    theme_bw()+ #背景变为白色
    ylim(-1500,3000) +
    theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
          panel.grid = element_blank()) + 
    theme(axis.text.x=element_text(family="Times",size=12,face="bold"), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
          axis.text.y=element_text(family="Times",size=12,face="bold"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
          axis.title.y=element_text(family="Times",size = 15,face="bold"), #设置y轴标题的字体属性
          axis.title.x=element_text(family="Times",size = 10,face="bold"),
          plot.title = element_text(hjust = 0.5,family="Times",size = 15,face="bold"),     
          legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold"))+
    labs(title="Immune Score", x="", y="Immune Score")+
    guides(fill='none')
  p2
  p3 <- ggplot(violin_dat, aes(x=group, y=ESTIMATEScore, fill=group))+
    geom_violin() + #绘制小提琴图
    stat_boxplot(geom="errorbar",
                 width=0.1,
                 position = position_dodge(0.9)) +
    geom_boxplot(width=0.4,
                 position=position_dodge(0.9),
                 outlier.shape = NA,
                 fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
    geom_point(aes(fill = group),
               size = 0.5,
               position = position_dodge(0.9))+
    scale_fill_manual(values = c("#FF6347","#6495ED"), name = "Group") +
    stat_compare_means(aes(group = group),method = 'wilcox')+ 
    theme_bw()+ #背景变为白色
    theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
          panel.grid = element_blank()) + 
    ylim(-2500, 4500) +
    theme(axis.text.x=element_text(family="Times",size=12,face="bold"), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
          axis.text.y=element_text(family="Times",size=12,face="bold"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
          axis.title.y=element_text(family="Times",size = 15,face="bold"), #设置y轴标题的字体属性
          axis.title.x=element_text(family="Times",size = 10,face="bold"),
          plot.title = element_text(hjust = 0.5,family="Times",size = 15,face="bold"),
          legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold"))+
    labs(title="ESTIMATE Score", x="", y="ESTIMATE Score")+
    guides(fill='none')
  p3
  p4 <- ggplot(violin_dat, aes(x=group, y=TumorPurity, fill=group))+
    geom_violin() + #绘制小提琴图
    stat_boxplot(geom="errorbar",
                 width=0.1,
                 position = position_dodge(0.9)) +
    geom_boxplot(width=0.4,
                 position=position_dodge(0.9),
                 outlier.shape = NA,
                 fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
    geom_point(aes(fill = group),
               size = 0.5,
               position = position_dodge(0.9))+
    scale_fill_manual(values = c("#FF6347","#6495ED"), name = "Group") +
    stat_compare_means(aes(group = group),method = 'wilcox')+ 
    theme_bw()+ #背景变为白色
    theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
          panel.grid = element_blank()) + 
    ylim(0.3, 1.2) +
    theme(axis.text.x=element_text(family="Times",size=12,face="bold"), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
          axis.text.y=element_text(family="Times",size=12,face="bold"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
          axis.title.y=element_text(family="Times",size = 15,face="bold"), #设置y轴标题的字体属性
          axis.title.x=element_text(family="Times",size = 10,face="bold"),
          plot.title = element_text(hjust = 0.5,family="Times",size = 15,face="bold"),
          legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold"))+
    labs(title="Tumor Purity", x="", y="Tumor Purity")+
    guides(fill='none')
  p4
  p5 <- cowplot::plot_grid(p1,p2,p3,p4,
                           nrow = 2, 
                           align = 'h', 
                           vjust = -0.3)
  p5
}
ggsave(filename = '01.estimate.all.pdf',p5,w=10,h=8)
ggsave(filename = '01.estimate.all.png',p5,w=10,h=8)


cor_plot_dat <- merge(violin_dat,group,by='sample')
cor_plot_dat$riskScore <- as.numeric(cor_plot_dat$riskScore)
library(ggstatsplot)
colnames(cor_plot_dat)
cor1 <- ggscatterstats(data = cor_plot_dat,
                                  x = riskScore,
                                  y = StromalScore,
                                  centrality.para = "mean",
                                  margins = "both",
                                  xfill = "#A73030FF",
                                  yfill = "#0073C2FF",
                                  type = "pearson",
                                  ylab = "StromalScore",
                                  marginal.type = "histogram",
                                  title = "Relationship between StromalScore and riskScore"
)+theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(size = 12))
cor1

colnames(cor_plot_dat)
cor2 <- ggscatterstats(data = cor_plot_dat,
                       x = riskScore,
                       y = ImmuneScore,
                       centrality.para = "mean",
                       margins = "both",
                       xfill = "#A73030FF",
                       yfill = "#0073C2FF",
                       type = "pearson",
                       ylab = "ImmuneScore",
                       marginal.type = "histogram",
                       title = "Relationship between ImmuneScore and riskScore"
)+theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(size = 12))
cor2

colnames(cor_plot_dat)
cor3 <- ggscatterstats(data = cor_plot_dat,
                       x = riskScore,
                       y = ESTIMATEScore,
                       centrality.para = "mean",
                       margins = "both",
                       xfill = "#A73030FF",
                       yfill = "#0073C2FF",
                       type = "pearson",
                       ylab = "ESTIMATEScore",
                       marginal.type = "histogram",
                       title = "Relationship between ESTIMATEScore and riskScore"
)+theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(size = 12))
cor3

colnames(cor_plot_dat)
cor4 <- ggscatterstats(data = cor_plot_dat,
                       x = riskScore,
                       y = TumorPurity,
                       centrality.para = "mean",
                       margins = "both",
                       xfill = "#A73030FF",
                       yfill = "#0073C2FF",
                       type = "pearson",
                       ylab = "TumorPurity",
                       marginal.type = "histogram",
                       title = "Relationship between TumorPurity and riskScore"
)+theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(size = 12))
cor4

cor5 <- cowplot::plot_grid(cor1,cor2,cor3,cor4,
                         nrow = 2, 
                         align = 'h', 
                         vjust = -0.3)
cor5

ggsave(filename = '02.cor.all.pdf',cor5,w=12,h=8)
ggsave(filename = '02.cor.all.png',cor5,w=13,h=9)
