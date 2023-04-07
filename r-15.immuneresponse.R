rm(list = ls())
setwd("/data/nas1/luchunlin/project/HZ0301-3/")
if (! dir.exists("./15_immueresponse")){
  dir.create("./15_immueresponse")
}
setwd("./15_immueresponse")
dat.tcga<-read.delim2('../00_rawdata/dat.fpkm.xls',row.names = 1)%>%lc.tableToNum()
colnames(dat.tcga)<-gsub('.','-',colnames(dat.tcga),fixed = T)
risk<-read.delim2('../07_risk/risk.xls')
dat.tcga <- dat.tcga[,risk$id]
high.sample<-risk$id[which(risk$risk==0)]
low.sample<-risk$id[which(risk$risk==1)]
dat <- log2(dat.tcga+1)
genelist <- read.table(file = 'GeneList.txt',header = T,sep = '\t')
geneset <- genelist[,c(1,6)]

gene_list <- split(as.matrix(geneset)[,1],
                   geneset[,2])
dat2 <- as.matrix(dat)
score = gsva(dat2, gene_list,
             method = "ssgsea",
             ssgsea.norm = TRUE,
             verbose = TRUE)
write.table(score,
            file = "immueresponse_score.xls",
            sep = "\t",
            quote = F)


hub_exp2<-score%>%as.data.frame()
hub_exp2$type<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = score,
                 -c("type"))
hub_exp2 <- na.omit(hub_exp2)
## 样本分组
hub_exp2$Group<-ifelse(hub_exp2$sample%in%low.sample,'Low_risk','High_risk')
##分面图形
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(rstatix)
library(lance)
stat.test<-hub_exp2%>%lc.tableToNum()%>%
  group_by(type)%>%
  wilcox_test(score ~ Group)%>%
  adjust_pvalue(method = 'fdr')
write.table(stat.test,file = 'wilcox.test.xls',sep = '\t',row.names = F,quote = F)
hub_exp2<-lc.tableToNum(hub_exp2)
violin_plot <- ggplot(hub_exp2, aes(x=type, 
                                    y=score,
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
  labs(title="ImmunGenset Score between Different Risk group", x="", y = "ImmunGenset Score",size=20) +
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
        axis.text.x=element_text(angle=75,hjust=1,colour="black",face="bold",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=8),
        legend.title = element_text(face = "bold", size = 10),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
violin_plot

ggsave('01.score.pdf',violin_plot,width = 11,height = 7)
ggsave('01.score.png',violin_plot,width = 11,height = 7)


## 生物标志物与基因集的相关性
# devtools::install_git("https://gitee.com/dr_yingli/ggcor")
# library(ggcor)
modelgene <- read.delim2('../06_Lasso/lasso_genes.csv',header = F)
gene.exp <- dat.tcga[modelgene$V1,]
gene.exp <- log2(gene.exp+1)
table(hub_exp2$type)

## 相关性----------
cor.dat <- score[,colnames(gene.exp)]
nc<-t(rbind(cor.dat,gene.exp))
m=rcorr(nc,type = 'spearman')$r[1:nrow(cor.dat),(ncol(nc)-length(modelgene$V1)+1):ncol(nc)]
m<-t(m)
p=rcorr(nc,type = 'spearman')$P[1:nrow(cor.dat),(ncol(nc)-length(modelgene$V1)+1):ncol(nc)]
p<-t(p)
library(dplyr)

tmp <- p
tmp[tmp<0.0001] <- '****'
tmp[tmp<0.001] <- '***'
tmp[tmp<0.01] <- '**'
tmp[tmp<0.05] <- '*'
tmp[tmp>0.05] <- 'ns'

cor <- m
cor <- signif(cor,3)
#cor[abs(cor)<0.1] <- ''

textMatrix = paste(cor,"\n",
                   tmp, sep = "")
dim(textMatrix) = dim(m)
colnames(textMatrix)<-colnames(m)
rownames(textMatrix)<-rownames(m)
library(WGCNA)
pdf(file = '02.correlation(modelgene).pdf',w=12,h=6)
par(mar = c(12, 10, 3, 1));
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
               main = paste("modelgene-DEcells correlation"))
dev.off()
png(file = '02.correlation(modelgene).png',w=800,h=400)
par(mar = c(12, 10, 3, 1));
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
               main = paste("modelgene-DEcells correlation"))
dev.off()
##最正相关和最负相关的
cor_plot_dat <- data.frame(sample=colnames(gene.exp),t(cor.dat[c('Interferons','TGFb_Family_Member_Receptor'),]),
                           t(gene.exp))
cor_plot_dat <- cor_plot_dat[risk$id,]
class(risk$riskScore)
risk$riskScore <- as.numeric(risk$riskScore)
cor_plot_dat$risk_group <- ifelse(risk$riskScore>median(risk$riskScore),'High risk','Low risk')
dim(cor_plot_dat)
library(ggplot2)
library(ggpubr)
## CXCL12-Antimicrobials--------
#install.packages('ggside')
#install.packages('ggstatsplot')
colnames(cor_plot_dat)
library(ggstatsplot)
cor1 <- ggscatterstats(data = cor_plot_dat,
                       x = MYB,
                       y = Interferons,
                       centrality.para = "mean",
                       margins = "both",
                       xfill = "#A73030FF",
                       yfill = "#0073C2FF",
                       type = "spearman",
                       ylab = "Interferons",
                       marginal.type = "histogram",
                       title = "Relationship between Interferons and MYB"
)+theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(size = 12))
cor1
ggsave(filename = "03.MYB_Interferons.png", height = 7, width = 8,cor1)
ggsave(filename = "03.MYB_Interferons.pdf", height = 7, width = 8,cor1)

##2-------
colnames(cor_plot_dat)
library(ggstatsplot)
colnames(cor_plot_dat)
cor2 <- ggscatterstats(data = cor_plot_dat,
                       x = LRFN4,
                       y = TGFb_Family_Member_Receptor,
                       centrality.para = "mean",
                       margins = "both",
                       xfill = "#A73030FF",
                       yfill = "#0073C2FF",
                       type = "spearman",
                       ylab = "TGFb_Family_Member_Receptor",
                       marginal.type = "histogram",
                       title = "Relationship between TGFb_Family_Member_Receptor and LRFN4"
)+theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(size = 12))
cor2
ggsave(filename = "04.LRFN4_TGFb_Family_Member_Receptor.png", height = 7, width = 8,cor2)
ggsave(filename = "04.LRFN4_TGFb_Family_Member_Receptor.pdf", height = 7, width = 8,cor2)

##差异免疫通路与风险评分相关性--------
risk<-read.delim2('../07_risk/risk.xls')
high.sample<-risk$id[which(risk$risk==0)]
low.sample<-risk$id[which(risk$risk==1)]
dat<-dat[,risk$id]
group<-data.frame(sample=risk$id,group=ifelse(risk$risk=='0','High risk','Low risk'))

DEpath <- score[stat.test$type[which(stat.test$p<0.05)],]
DEpath <- t(DEpath)%>%as.data.frame()

DEpath <- DEpath[risk$id,]

library(psych)
x <-as.numeric(risk$riskScore)
y <-DEpath
library(psych)
d <- corr.test(y,x,use="complete",method = 'spearman')
r <- data.frame(d$r)
p <- data.frame(d$p)
correlation<-data.frame(rownames(p),r$d.r,p$d.p)
colnames(correlation) <- c("cell","Correlation","p.value")
write.table(correlation,file = 'correlation.xls',sep = '\t',row.names = F,quote = F)
correlation<-correlation[correlation$p.value<0.05,]
correlation$sig[correlation$p.value>0.05] <- "ns"   
correlation$sig[correlation$p.value<0.05&correlation$p.value>0.01] <- "*"   
correlation$sig[correlation$p.value<0.01&correlation$p.value>0.001] <- "**"  
correlation$sig[correlation$p.value<0.001&correlation$p.value>0.0001] <- "***"   
correlation$sig[correlation$p.value<0.0001] <- "****"
correlation$'correlation'<-abs(correlation$Correlation)
#correlation_sig<-correlation[c(3,9,10,11,23,29,7,12,18,31,1,14,15,25,28),]
library("viridis")
#trace(ggdotchart,edit=T)
p0<-ggdotchart(correlation, x = "cell", y = "Correlation",
               
               dot.size ='correlation',
               # filled.contour(p.value,col=Lab.palette(20)),
               color ='p.value',
               # label='sig',
               
               #  font.label = list( size = 9),
               
               #ylab="mgp = c(3, 2, 0)",
               # font.label = list(size=10,face="bold",color='black',position="dodge",
               #                   vjust=0.5),
               #  y.label=0.7,
               # palette = palette,
               # 按照cyl填充颜色
               # palette = palette, # 修改颜色
               sorting = "descending",
               add = "segments",                             # 添加棒子
               rotate = TRUE,
               ggtheme = theme_pubr(),                        # 改变主题
               ylab="Correlation Coefficient (r)",
               xlab='',
               # dot.size = 6 ,
               title='riskScore'
               
)+scale_colour_gradient( high = "#4682B4",low = "#66CDAA")

p10<-p0+theme(legend.position = "right",
              panel.background = element_blank())+geom_hline(aes(yintercept=0),linetype="dashed",lwd = 0.2)+
  theme(axis.title.x =element_text(size=15,family = "Times", face = "bold"),
        axis.text.x =element_text(size=12,family = "Times", face = "bold"),
        axis.title.y =element_text(size=20,family = "Times", face = "bold"),
        axis.text.y=element_text(size=12,family = "Times", face = "bold"),
        plot.title=element_text(size=20,family = "Times", face = "bold",hjust=0.5),
        legend.text = element_text(size = 16, family = "Times"),
        legend.title = element_text(size = 18, family = "Times",face = "bold"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

p10
ggsave(filename = '05.cor.pdf',p10,w=8,h=5)  
ggsave(filename = '05.cor.png',p10,w=8,h=5)

