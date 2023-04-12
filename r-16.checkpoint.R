rm(list = ls())
setwd("/data/nas1/luchunlin/project/SJZZK-428-10/")
if (! dir.exists("./16_checkpoint")){
  dir.create("./16_checkpoint")
}
setwd("./16_checkpoint")
dat.tcga<-read.delim2("../00_rawdata/dat.fpkm.xls", row.names = 1)%>% lc.tableToNum
colnames(dat.tcga)<-gsub('.','-',colnames(dat.tcga),fixed = T)
risk<-read.delim2('../09_risk/risk.xls')
high.sample<-risk$id[which(risk$risk==0)]
low.sample<-risk$id[which(risk$risk==1)]
risk2 <- risk
risk2$risk_label <- ifelse(risk$risk == 0, "High", "Low")
dat.tcga<-dat.tcga[,risk$id]
dim(dat.tcga)
## 把检验点基因提取出来
checkpoint <- read.table("/data/nas1/luchunlin/pipeline/Checkpoint/checkpoint.txt",
                         header = F)
checkpoint <- checkpoint$V1
length(checkpoint)
checkpoint_dat <- dat.tcga[which(rownames(dat.tcga)%in%checkpoint),]
dim(checkpoint_dat)
checkpoint_dat<-log2(checkpoint_dat+1)
checkpoint_dat$gene<-rownames(checkpoint_dat)

##PART A expression-----------
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
write.table(stat.test,file = 'wilcox_result.xls',
            sep = '\t',
            row.names = F)
DE.checkpoint <- subset(stat.test,stat.test$p < 0.05)
##14
violin_dat<-violin_dat[violin_dat$gene%in%DE.checkpoint$gene,]

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
  labs(title="Immune Checkpoint", x="", y = "log2(expr+1)",size=20) +
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
        panel.grid.minor = element_blank())#+facet_wrap(~gene,scales = "free",nrow = 3) 
violin_plot
ggsave('01.checkpoint.pdf',violin_plot,w=7,h=5)
ggsave('01.checkpoint.png',violin_plot,w=7,h=5)

## PART B correlation--------
## 添加相关性
#表达值进行log(1+)转化，使数据更服从正态分布，减少离散度极大值影响
#hubgene<-read.delim2('/data/nas1/luchunlin/project/YQ444-8/12_Lasso/lasso_genes.csv',header = F)
#dat.final<-dat.tcga[hubgene$V1,]
#dat.final<-log2(dat.final+1)
riskscore<-risk[,"riskScore"]%>%as.data.frame()
riskscore<-as.numeric(riskscore$.)
riskscore<-t(riskscore)
colnames(riskscore)<-risk$id
rownames(riskscore)<-'riskscore'
DE.checkpoint_dat <- checkpoint_dat[DE.checkpoint$gene,]
corr.dat<-t(rbind(DE.checkpoint_dat[,-379],riskscore))
#基因表达值的相关性分析，以Pearson相关系数为例
gene_cor <- cor(corr.dat, method = 'spearman')

#去除基因的自相关，也就是对角线的值
diag(gene_cor) <- 0
gene_cor<-gene_cor[,-c(1:14)]
#gene_cor<-gene_cor[c(1:46),]
gene_cor  #最终的基因间表达值Pearson相关性矩阵
#将获得的相关性矩阵转换为两两对应的数据框结构
gene_cor <- reshape2::melt(gene_cor)
gene_cor <- subset(gene_cor, value != 0)  #去除0值的相关性
head(gene_cor)  #前两列是两个基因名称，第三列为两个基因的相关性
gene_cor$Var1<-rownames(gene_cor)
gene_cor$Var2<-rep('riskscore',14)
gene_cor<-gene_cor[,c(2,3,1)]
class(gene_cor$value)
## abs >3的
#gene_cor<-gene_cor[which(abs(gene_cor$value)>0.3),]

library(circlize)
library(ComplexHeatmap)
pdf(file = '02.correlation.pdf',w=6,h=6)
circos.par("track.height" = 0.2)
chordDiagram(gene_cor, 
             annotationTrack = c('grid', 'name', 'axis'), #绘制外周圆弧区，显示名称和刻度轴
             col = colorRamp2(c(-1, 0, 1), c('green', 'white', 'red'), transparency = 0.5), #根据相关性大小展示连线的颜色范围
             annotationTrackHeight = c(0.05, 0.05), #名称离圆弧的距离，以及圆弧的宽度
)
# next the grid graphics are added directly to the plot
col_fun = colorRamp2(c(-1,0,1), c("green", "white", "red"))
# where circlize has created.
lgd_links = Legend(at = c(-1,0,1), col_fun = col_fun, 
                   title_position = "topleft", title = "correlation")
draw(lgd_links, x = unit(5, "mm"), y = unit(5, "mm"), just= c("left", "bottom"))
dev.off()

png(file = '02.correlation.png',w=500,h=500)
circos.par("track.height" = 0.2)
chordDiagram(gene_cor, 
             annotationTrack = c('grid', 'name', 'axis'), #绘制外周圆弧区，显示名称和刻度轴
             col = colorRamp2(c(-1, 0, 1), c('green', 'white', 'red'), transparency = 0.5), #根据相关性大小展示连线的颜色范围
             annotationTrackHeight = c(0.05, 0.05), #名称离圆弧的距离，以及圆弧的宽度
)
# next the grid graphics are added directly to the plot
col_fun = colorRamp2(c(-1,0,1), c("green", "white", "red"))
# where circlize has created.
lgd_links = Legend(at = c(-1,0,1), col_fun = col_fun, 
                   title_position = "topleft", title = "correlation")
draw(lgd_links, x = unit(5, "mm"), y = unit(5, "mm"), just= c("left", "bottom"))
dev.off()
