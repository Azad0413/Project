rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-300-8/")
if (! dir.exists("./19_HLA")){
  dir.create("./19_HLA")
}
setwd("./19_HLA")
dat.tcga<-read.delim2('../00_rawdata/dat.fpkm.xls',row.names = 1)%>%lc.tableToNum()
colnames(dat.tcga)<-gsub('.','-',colnames(dat.tcga),fixed = T)
risk<-read.delim2('../08_risk/risk.xls')
dat.tcga <- dat.tcga[,risk$id]
high.sample<-risk$id[which(risk$risk==0)]
low.sample<-risk$id[which(risk$risk==1)]
risk2 <- risk
risk2$risk_label <- ifelse(risk$risk == 0, "High", "Low")
dim(dat.tcga)
## 把HLA基因表达矩阵提取出来
library(readxl)
HLA <- read_xlsx('HLA_family.xlsx')
HLA.dat<-dat.tcga[HLA$Symbol,]
HLA.dat<-log2(HLA.dat+1)
HLA.dat <- na.omit(HLA.dat)
HLA.dat2 <- HLA.dat
HLA.dat2$gene<-rownames(HLA.dat2)
violin_dat <- gather(HLA.dat2, key=sample, value='expr', -c("gene"))
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
write.table(stat.test,file = 'HLA_exp_result.xls',
            sep = '\t',
            row.names = F)
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
  labs(title="Expression of HLA family genes", x="", y = "log2(expr+1)",size=20) +
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
        panel.grid.minor = element_blank())
violin_plot
ggsave('01.HLA_expression.pdf',violin_plot,w=9,h=5)
ggsave('01.HLA_expression.png',violin_plot,w=9,h=5)


modelgene <- read.delim2('../06_univariate_cox/univariate_cox_result_0.05.xls')
gene.exp <- dat.tcga[rownames(modelgene),]
gene.exp <- log2(gene.exp+1)%>%t
hubgene <- data.frame(hubgene=rownames(modelgene))
gene.exp <- t(gene.exp)

nc<-t(rbind(HLA.dat,gene.exp))
m=rcorr(nc)$r[1:nrow(HLA.dat),(ncol(nc)-length(hubgene$hubgene)+1):ncol(nc)]
m<-t(m)
p=rcorr(nc)$P[1:nrow(HLA.dat),(ncol(nc)-length(hubgene$hubgene)+1):ncol(nc)]
p<-t(p)
library(dplyr)
library(dplyr)


tmp = matrix(case_when(p<0.0001~"****",
                       p<0.0001~"***",
                       p<0.01~"**",
                       p<0.05~"*",
                       T~""),nrow = nrow(p))
cor <- m
cor <- signif(cor,3)
cor[abs(cor)<0.2] <- ''

textMatrix = paste(cor,"\n",
                   tmp, sep = "")
dim(textMatrix) = dim(m)
colnames(textMatrix)<-colnames(m)
rownames(textMatrix)<-rownames(m)
#source("modified_pheatmap.R")
#library(pheatmap)
#trace(pheatmap,edit = T)
pdf(file = '02.heatmap.pdf',w=15,h=5)
pheatmap(m,
         display_numbers =textMatrix,
         angle_col =45,
         color = colorRampPalette(c("#92b7d1", "white", "#d71e22"))(100),
         border_color = "white",
         treeheight_col = 0,
         treeheight_row = 0,
         fontsize_row = 12,
         fontsize_col = 12,
         fontsize = 15)
dev.off()
png(file = '02.heatmap.png',w=1400,h=400)
pheatmap(m,
         display_numbers =textMatrix,
         angle_col =45,
         color = colorRampPalette(c("#92b7d1", "white", "#d71e22"))(100),
         border_color = "white",
         treeheight_col = 0,
         treeheight_row = 0,
         fontsize_row = 12,
         fontsize_col = 12,
         fontsize = 15)
dev.off()

##最正相关和最负相关的
cor_plot_dat <- data.frame(sample=colnames(gene.exp),t(HLA.dat[c('HLA-DOA','HLA-DMA'),]),
                           t(gene.exp))
cor_plot_dat <- cor_plot_dat[risk$id,]
cor_plot_dat$risk_group <- ifelse(risk$riskScore>median(risk$riskScore),'High risk','Low risk')
dim(cor_plot_dat)
library(ggplot2)
library(ggpubr)
## CXCL12-HLA-DOA--------
#install.packages('ggside')
#install.packages('ggstatsplot')
colnames(cor_plot_dat) <- gsub('.','-',colnames(cor_plot_dat),fixed = T)
library(ggstatsplot)
cor1 <- ggscatterstats(data = cor_plot_dat,
                       x = CXCL12,
                       y = 'HLA-DOA',
                       centrality.para = "mean",
                       margins = "both",
                       xfill = "#A73030FF",
                       yfill = "#0073C2FF",
                       type = "pearson",
                       ylab = "HLA-DOA",
                       marginal.type = "histogram",
                       title = "Relationship between HLA-DOA and CXCL12"
)+theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(size = 12))
cor1
ggsave(filename = "03.CXCL12_HLA-DOA.png", height = 7, width = 8,cor1)
ggsave(filename = "03.CXCL12_HLA-DOA.pdf", height = 7, width = 8,cor1)
##高低风险组CXCL12-HLA-DOA的表达差异
exp1 <- t(cor_plot_dat[,c(2,5)])%>%as.data.frame()
exp1$gene <- rownames(exp1)
violin_dat1 <- gather(exp1, key=sample, value='expr', -c("gene"))
head(violin_dat1)
violin_dat1$group <- ifelse(violin_dat1$sample %in% high.sample,
                            "High", "Low") 
head(violin_dat1)
colnames(violin_dat1)
library(rstatix)
stat.test<-violin_dat1%>%
  group_by(gene)%>%
  wilcox_test(expr ~ group)%>%
  adjust_pvalue(method = 'fdr')
write.table(stat.test,file = 'wilcox_result(CXCL12-HLA-DOA).xls',
            sep = '\t',
            row.names = F)
# violin_dat<-violin_dat[violin_dat$gene%in%DE.checkpoint$gene,]

violin_plot1 <- ggplot(violin_dat1, aes(x=group, 
                                        y=expr,
                                        fill=group)) +
  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.3,
               position=position_dodge(0.9),
               outlier.shape = NA,fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#A73030FF", "#0073C2FF"), name = "Group")+
  labs(title="", x="", y = "log2(expr+1)",size=20) +
  stat_compare_means(data = violin_dat1,
                     mapping = aes(group = group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F,label.x = 1.4) +
  # geom_signif(comparisons = list(c('High','Low')),
  #             test = wilcox.test,
  #             map_signif_level = T,
  #             y_position = c(4))+
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
  facet_wrap(~gene,scales = "free",nrow = 1) 
violin_plot1
ggsave('04.expression(CXCL12_HLA-DOA).pdf',violin_plot1,w=6,h=4.5)
ggsave('04.expression(CXCL12_HLA-DOA).png',violin_plot1,w=6,h=4.5)

## BNIP3- HLA-DMA-------
colnames(cor_plot_dat)
library(ggstatsplot)
cor2 <- ggscatterstats(data = cor_plot_dat,
                       x = BNIP3,
                       y = 'HLA-DMA',
                       centrality.para = "mean",
                       margins = "both",
                       xfill = "#A73030FF",
                       yfill = "#0073C2FF",
                       type = "pearson",
                       ylab = "HLA-DMA",
                       marginal.type = "histogram",
                       title = "Relationship between HLA-DMA and BNIP3"
)+theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(size = 12))
cor2
ggsave(filename = "05.BNIP3_HLA-DMA.png", height = 7, width = 8,cor2)
ggsave(filename = "05.BNIP3_HLA-DMA.pdf", height = 7, width = 8,cor2)
##高低风险组BNIP3-HLA-DMA的表达差异
exp2 <- t(cor_plot_dat[,c(3,4)])%>%as.data.frame()
exp2$gene <- rownames(exp2)
violin_dat2 <- gather(exp2, key=sample, value='expr', -c("gene"))
head(violin_dat2)
violin_dat2$group <- ifelse(violin_dat2$sample %in% high.sample,
                            "High", "Low") 
head(violin_dat2)
colnames(violin_dat2)
library(rstatix)
stat.test<-violin_dat2%>%
  group_by(gene)%>%
  wilcox_test(expr ~ group)%>%
  adjust_pvalue(method = 'fdr')
write.table(stat.test,file = 'wilcox_result(BNIP3-HLA-DMA).xls',
            sep = '\t',
            row.names = F)
# violin_dat<-violin_dat[violin_dat$gene%in%DE.checkpoint$gene,]

violin_plot2 <- ggplot(violin_dat2, aes(x=group, 
                                        y=expr,
                                        fill=group)) +
  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.3,
               position=position_dodge(0.9),
               outlier.shape = NA,fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#A73030FF", "#0073C2FF"), name = "Group")+
  labs(title="", x="", y = "log2(expr+1)",size=20) +
  stat_compare_means(data = violin_dat2,
                     mapping = aes(group = group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F,label.x = 1.4) +
  # geom_signif(comparisons = list(c('High','Low')),
  #             test = wilcox.test,
  #             map_signif_level = T,
  #             y_position = c(4))+
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
  facet_wrap(~gene,scales = "free",nrow = 1) 
violin_plot2
ggsave('06.expression(BNIP3_HLA-DMA).pdf',violin_plot2,w=6,h=4.5)
ggsave('06.expression(BNIP3_HLA-DMA).png',violin_plot2,w=6,h=4.5)
