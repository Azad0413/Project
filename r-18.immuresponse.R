rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-300-8/")
if (! dir.exists("./18_immueresponse")){
  dir.create("./18_immueresponse")
}
setwd("./18_immueresponse")
dat.tcga<-read.delim2('../00_rawdata/dat.fpkm.xls',row.names = 1)%>%lc.tableToNum()
colnames(dat.tcga)<-gsub('.','-',colnames(dat.tcga),fixed = T)
risk<-read.delim2('../08_risk/risk.xls')
dat.tcga <- dat.tcga[,risk$id]
high.sample<-risk$id[which(risk$risk==0)]
low.sample<-risk$id[which(risk$risk==1)]

genelist <- read.table(file = 'GeneList.txt',header = T,sep = '\t')
hub.exp <- dat.tcga[genelist$Symbol,]
hub_exp2<-log2(hub.exp+1)
hub_exp2$type<-genelist$Category
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
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
  wilcox_test(expr ~ Group)%>%
  adjust_pvalue(method = 'fdr')
write.table(stat.test,file = '01.wilcox.test.xls',sep = '\t',row.names = F,quote = F)
hub_exp2<-lc.tableToNum(hub_exp2)
violin_plot <- ggplot(hub_exp2, aes(x=type, 
                                    y=expr,
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

ggsave('01.expression.pdf',violin_plot,width = 11,height = 7)
ggsave('01.expression.png',violin_plot,width = 11,height = 7)


## 生物标志物与基因集的相关性
#devtools::install_git("https://gitee.com/dr_yingli/ggcor")
library(ggcor)
modelgene <- read.delim2('../06_univariate_cox/univariate_cox_result_0.05.xls')
gene.exp <- dat.tcga[rownames(modelgene),]
gene.exp <- log2(gene.exp+1)%>%t
table(hub_exp2$type)
type <- data.frame(table(hub_exp2$type))
hub.exp3 <- log2(hub.exp+1)
hub.exp3$type <- genelist$Category
hub.exp3 <- na.omit(hub.exp3)
Antigen_Processing_and_Presentation <- hub.exp3[which(hub.exp3$type=='Antigen_Processing_and_Presentation'),]%>%select(-type)
Antimicrobials <- hub.exp3[which(hub.exp3$type=='Antimicrobials'),]%>%select(-type)
BCRSignalingPathway <- hub.exp3[which(hub.exp3$type=='BCRSignalingPathway'),]%>%select(-type)
Chemokine_Receptors <- hub.exp3[which(hub.exp3$type=='Chemokine_Receptors'),]%>%select(-type)
Chemokines <- hub.exp3[which(hub.exp3$type=='Chemokines'),]%>%select(-type)
Cytokine_Receptors <- hub.exp3[which(hub.exp3$type=='Cytokine_Receptors'),]%>%select(-type)
Cytokines <- hub.exp3[which(hub.exp3$type=='Cytokines'),]%>%select(-type)
Interferon_Receptor <- hub.exp3[which(hub.exp3$type=='Interferon_Receptor'),]%>%select(-type)
Interferons <- hub.exp3[which(hub.exp3$type=='Interferons'),]%>%select(-type)
Interleukins <- hub.exp3[which(hub.exp3$type=='Interleukins'),]%>%select(-type)
Interleukins_Receptor <- hub.exp3[which(hub.exp3$type=='Interleukins_Receptor'),]%>%select(-type)
NaturalKiller_Cell_Cytotoxicity <- hub.exp3[which(hub.exp3$type=='NaturalKiller_Cell_Cytotoxicity'),]%>%select(-type)
TCRsignalingPathway <- hub.exp3[which(hub.exp3$type=='TCRsignalingPathway'),]%>%select(-type)
TGFb_Family_Member <- hub.exp3[which(hub.exp3$type=='TGFb_Family_Member'),]%>%select(-type)
TGFb_Family_Member_Receptor <- hub.exp3[which(hub.exp3$type=='TGFb_Family_Member_Receptor'),]%>%select(-type)
TNF_Family_Members <- hub.exp3[which(hub.exp3$type=='TNF_Family_Members'),]%>%select(-type)
TNF_Family_Members_Receptors <- hub.exp3[which(hub.exp3$type=='TNF_Family_Members_Receptors'),]%>%select(-type)
colMeans(Antigen_Processing_and_Presentation)

cor.dat <- data.frame(row.names = type$Var1,
                      rbind(colMeans(Antigen_Processing_and_Presentation),
                            colMeans(Antimicrobials),
                            colMeans(BCRSignalingPathway),
                            colMeans(Chemokine_Receptors),
                            colMeans(Chemokines),
                            colMeans(Cytokine_Receptors),
                            colMeans(Cytokines),
                            colMeans(Interferon_Receptor),
                            colMeans(Interferons),
                            colMeans(Interleukins),
                            colMeans(Interleukins_Receptor),
                            colMeans(NaturalKiller_Cell_Cytotoxicity),
                            colMeans(TCRsignalingPathway),
                            colMeans(TGFb_Family_Member),
                            colMeans(TGFb_Family_Member_Receptor),
                            colMeans(TNF_Family_Members),
                            colMeans(TNF_Family_Members_Receptors)))
colnames(cor.dat) <- gsub('.','-',colnames(cor.dat),fixed = T)
## 相关性----------

colnames(gene.exp)
hubgene <- data.frame(hubgene=rownames(modelgene))
gene.exp <- t(gene.exp)
nc<-t(rbind(cor.dat,gene.exp))
m=rcorr(nc)$r[1:nrow(cor.dat),(ncol(nc)-length(hubgene$hubgene)+1):ncol(nc)]
m<-t(m)
p=rcorr(nc)$P[1:nrow(cor.dat),(ncol(nc)-length(hubgene$hubgene)+1):ncol(nc)]
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
cor[abs(cor)<0.21] <- ''

textMatrix = paste(cor,"\n",
                   tmp, sep = "")
dim(textMatrix) = dim(m)
colnames(textMatrix)<-colnames(m)
rownames(textMatrix)<-rownames(m)
source("modified_pheatmap.R")
library(pheatmap)
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
cor_plot_dat <- data.frame(sample=colnames(gene.exp),t(cor.dat[c('Antimicrobials','Antigen_Processing_and_Presentation'),]),
                           t(gene.exp))
cor_plot_dat <- cor_plot_dat[risk$id,]
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
                                  x = CXCL12,
                                  y = Antimicrobials,
                                  centrality.para = "mean",
                                  margins = "both",
                                  xfill = "#A73030FF",
                                  yfill = "#0073C2FF",
                                  type = "pearson",
                                  ylab = "Antimicrobials",
                                  marginal.type = "histogram",
                                  title = "Relationship between Antimicrobials and CXCL12"
)+theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(size = 12))
cor1
ggsave(filename = "03.CXCL12_Antimicrobials.png", height = 7, width = 8,cor1)
ggsave(filename = "03.CXCL12_Antimicrobials.pdf", height = 7, width = 8,cor1)
##高低风险组CXCL12-Antimicrobials的表达差异
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
write.table(stat.test,file = 'wilcox_result(CXCL12-Antimicrobials).xls',
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
ggsave('04.expression(CXCL12_Antimicrobials).pdf',violin_plot1,w=6,h=4.5)
ggsave('04.expression(CXCL12_Antimicrobials).png',violin_plot1,w=6,h=4.5)

## BNIP3- Antigen_Processing_and_Presentation-------
colnames(cor_plot_dat)
library(ggstatsplot)
cor2 <- ggscatterstats(data = cor_plot_dat,
                       x = BNIP3,
                       y = Antigen_Processing_and_Presentation,
                       centrality.para = "mean",
                       margins = "both",
                       xfill = "#A73030FF",
                       yfill = "#0073C2FF",
                       type = "pearson",
                       ylab = "Antigen_Processing_and_Presentation",
                       marginal.type = "histogram",
                       title = "Relationship between Antigen_Processing_and_Presentation and BNIP3"
)+theme(axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        plot.title = element_text(size = 12))
cor2
ggsave(filename = "05.BNIP3_Antigen_Processing_and_Presentation.png", height = 7, width = 8,cor2)
ggsave(filename = "05.BNIP3_Antigen_Processing_and_Presentation.pdf", height = 7, width = 8,cor2)
##高低风险组BNIP3-Antigen_Processing_and_Presentation的表达差异
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
write.table(stat.test,file = 'wilcox_result(BNIP3-Antigen_Processing_and_Presentation).xls',
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
ggsave('06.expression(BNIP3_Antigen_Processing_and_Presentation).pdf',violin_plot2,w=6,h=4.5)
ggsave('06.expression(BNIP3_Antigen_Processing_and_Presentation).png',violin_plot2,w=6,h=4.5)
