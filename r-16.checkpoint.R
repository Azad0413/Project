rm(list = ls())
setwd("/data/nas1/luchunlin/project/HZ0301-3/")
if (! dir.exists("./16_checkpoint")){
  dir.create("./16_checkpoint")
}
setwd("./16_checkpoint")
dat.tcga<-read.csv("../00_rawdata/dat.fpkm.xls", row.names = 1,sep = '\t')
colnames(dat.tcga)<-gsub('.','-',colnames(dat.tcga),fixed = T)
risk<-read.delim2('../07_risk/risk.xls')
high.sample<-risk$id[which(risk$risk==0)]
low.sample<-risk$id[which(risk$risk==1)]
risk2 <- risk
risk2$risk_label <- ifelse(risk$risk == 0, "High", "Low")
dat.tcga<-dat.tcga[,risk$id]
dim(dat.tcga)
## 把检验点基因提取出来
## 把检验点基因提取出来
checkpoint <- read.table("/data/nas1/luchunlin/pipeline/Checkpoint/checkpoint48.txt",
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
##19
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
ggsave('01.checkpoint.pdf',violin_plot,w=9,h=5)
ggsave('01.checkpoint.png',violin_plot,w=9,h=5)
modelgene <- read.delim2('../06_Lasso/lasso_genes.csv',header = F)
gene.exp <- dat.tcga[modelgene$V1,]
gene.exp <- log2(gene.exp+1)
table(hub_exp2$type)

## 相关性----------
cor.dat <- checkpoint_dat[DE.checkpoint$gene,]%>%select(-'gene')
lasso.gene <- read.delim2('../06_Lasso/lasso_genes.csv',header = F)
gene.exp <- dat.tcga[lasso.gene$V1,]
gene.exp <- log2(gene.exp+1)
gene.exp <- gene.exp[,colnames(cor.dat)]
nc<-t(rbind(cor.dat,gene.exp))
m=rcorr(nc,type = 'spearman')$r[1:nrow(cor.dat),(ncol(nc)-length(lasso.gene$V1)+1):ncol(nc)]
m<-t(m)
p=rcorr(nc,type = 'spearman')$P[1:nrow(cor.dat),(ncol(nc)-length(lasso.gene$V1)+1):ncol(nc)]
p<-t(p)
library(dplyr)

tmp <- p
tmp[tmp<0.0001] <- '****'
tmp[tmp<0.001] <- '***'
tmp[tmp<0.01] <- '**'
tmp[tmp<0.05] <- '*'
tmp[tmp>0.05] <- 'ns'

cor <- m
# cor <- signif(cor,3)
cor <- round(cor,digits = 3)
#cor[abs(cor)<0.1] <- ''

textMatrix = paste(cor,"\n",
                   tmp, sep = "")
dim(textMatrix) = dim(m)
colnames(textMatrix)<-colnames(m)
rownames(textMatrix)<-rownames(m)
library(WGCNA)
pdf(file = '02.correlation.pdf',w=12,h=6)
par(mar = c(6, 6, 3, 1));
labeledHeatmap(Matrix = m, 
               xLabels = colnames(m), 
               yLabels = rownames(m), 
               cex.lab = 1, 
               # ySymbols = colnames(MEs_col)[-20], 
               colorLabels = FALSE, 
               colors = blueWhiteRed(50, 0.8),
               textMatrix = textMatrix, 
               setStdMargins = FALSE, 
               cex.text = 0.8, 
               zlim = c(-1,1),
               main = paste("modelgene-checkpoint correlation"))
dev.off()
png(file = '02.correlation.png',w=800,h=400)
par(mar = c(6, 6, 3, 1));
labeledHeatmap(Matrix = m, 
               xLabels = colnames(m), 
               yLabels = rownames(m), 
               cex.lab = 1, 
               # ySymbols = colnames(MEs_col)[-20], 
               colorLabels = FALSE, 
               colors = blueWhiteRed(50, 0.8),
               textMatrix = textMatrix, 
               setStdMargins = FALSE, 
               cex.text = 0.8, 
               zlim = c(-1,1),
               main = paste("modelgene-checkpoint correlation"))
dev.off()
