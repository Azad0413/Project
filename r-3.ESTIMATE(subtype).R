rm(list = ls())
setwd("/data/nas1/luchunlin/project/HF-0103-1/")
if (! dir.exists("./03_ESTIMATE(subtype)")){
  dir.create("./03_ESTIMATE(subtype)")
}
setwd("./03_ESTIMATE(subtype)")

dat<-read.delim2("../00_rawdata/dat(GSE10846).xls", row.names = 1)%>% lc.tableToNum
group<-read.delim2('../01_subtype/cluster.xls')
table(group$cluster)
colnames(group)<-c('sample','group')
cluster1_sample <- group$sample[which(group$group=='cluster 1')]
cluster2_sample <- group$sample[which(group$group=='cluster 2')]

dat<-dat[,group$sample]
group_estimate<-group$group%>%as.factor()
design<-model.matrix(~0 + group_estimate)
rownames(design)<-group$sample
colnames(design)<-levels(group_estimate)
design<-as.data.frame(design)
cluster1<-rownames(design)[which(design$`cluster 1`==1)]
cluster2<-rownames(design)[which(design$`cluster 2`==1)]
# cluster3<-rownames(design)[which(design$`cluster 3`==1)]

length(cluster1)
length(cluster2)
# length(cluster3)
#install.packages("estimate", repos="http://R-Forge.R-project.org")
library(estimate)
# expr_train <- log2(dat+1)
expr_train <- dat
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
violin_dat$group <- ifelse(violin_dat$sample %in% cluster1_sample,"cluster 1",
                           ifelse(violin_dat$sample%in%cluster2_sample,"cluster 2","cluster 3"))
library(tidyr)
library(ggplot2)
library(ggpubr)
library(Ipaper)
# Kruskal-Wallis
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
    scale_fill_manual(values = c("#FF6347","#6495ED","purple"), name = "Group") + 
    stat_compare_means(aes(group = group),method = 'kruskal.test')+ 
    theme_bw()+ #背景变为白色
    ylim(-1500,2500) +
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
    scale_fill_manual(values = c("#FF6347","#6495ED","purple"), name = "Group") +
    stat_compare_means(aes(group = group),method = 'kruskal.test')+ 
    theme_bw()+ #背景变为白色
    ylim(0,5000) +
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
    scale_fill_manual(values = c("#FF6347","#6495ED","purple"), name = "Group") +
    stat_compare_means(aes(group = group),method = 'kruskal.test')+ 
    theme_bw()+ #背景变为白色
    theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
          panel.grid = element_blank()) + 
    ylim(500, 6500) +
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
    scale_fill_manual(values = c("#FF6347","#6495ED","purple"), name = "Group") +
    stat_compare_means(aes(group = group),method = 'kruskal.test')+ 
    theme_bw()+ #背景变为白色
    theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
          panel.grid = element_blank()) + 
    ylim(0, 1) +
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
