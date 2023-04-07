rm(list = ls())
setwd("/data/nas1/luchunlin/project/HZ0301-3/")
if (! dir.exists("./12_cancerpathway")){
  dir.create("./12_cancerpathway")
}
setwd("./12_cancerpathway")
library(lance)
library(tidyverse)
dat.tcga<-read.delim2('../00_rawdata/dat.fpkm.xls',row.names = 1)%>%lc.tableToNum()
colnames(dat.tcga)<-gsub('.','-',colnames(dat.tcga),fixed = T)
risk<-read.delim2('../07_risk/risk.xls')
dat.tcga <- dat.tcga[,risk$id]
high.sample<-risk$id[which(risk$risk==0)]
low.sample<-risk$id[which(risk$risk==1)]
dat <- log2(dat.tcga+1)
# dat <- dat.tcga
geneset <- readxl::read_xlsx('13cancerpathway.xlsx')
colnames(geneset)

gene_list <- split(as.matrix(geneset)[,2],
                   geneset[,1])
dat2 <- as.matrix(dat)
library(GSVA)
score = gsva(dat2, gene_list,
             method = "ssgsea",
             ssgsea.norm = TRUE,
             verbose = TRUE)
write.table(score,
            file = "cancerpathway_score.xls",
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
  labs(title="Cancer-related Pathway Score between Different Risk group", x="", y = "Cancer-related Pathway Score",size=20) +
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

ggsave('01.score.pdf',violin_plot,width = 10,height = 6)
ggsave('01.score.png',violin_plot,width = 10,height = 6)


# 绘图 ----------------------------------------------------------------------
#devtools::install_git("https://gitee.com/dr_yingli/ggcor") 
library(ggcor)
library(vegan)
score1<- t(score[,risk$id])%>%as.data.frame()%>%lc.tableToNum()
riskScore <- as.numeric(risk$riskScore)
names(riskScore) <- risk$id
riskScore <- as.data.frame(riskScore)

class(riskScore$riskScore)
colnames(riskScore)

library(Hmisc)
cor.dat <- t(score1)
riskScore <- t(riskScore)
nc<-t(rbind(cor.dat,riskScore))

class(nc)
m=rcorr(nc,type = 'spearman')$r[1:nrow(cor.dat),ncol(nc)]

p=rcorr(nc,type = 'spearman')$P[1:nrow(cor.dat),ncol(nc)]


correlation <- data.frame(sepc='riskScore',env=rownames(cor.dat),r=m,p=p)%>%
  mutate(r_value = cut(r, breaks = c(-Inf,-0.5,-0.25,0.25, 0.5, Inf), 
                       labels = c('<=-0.5','-0.5- -0.25','-0.25-0.25','0.25-0.5', '>=0.5'), right = FALSE),
         p_value = cut(p, breaks = c(-Inf, 0.001, 0.01, 0.05, 1),   ##Inf最大
                       labels = c('<0.001', '<0.01', '<0.05', 'ns'), right = FALSE))

library(ggcor)

library(psych)
corr <- corr.test(score1,method = 'spearman')

pdf('02.heatmap.pdf',w=10,h=10)

quickcor(corr,type = "lower") +
  geom_circle2()+
  scale_fill_gradient2(low = "#9483E1",mid = "white",high = "#E11953",midpoint=0) +
  geom_mark(sig.thres = 0.05, size = 2.5, colour = "black")+
  # geom_colour()+
  #geom_square() +
  anno_link(aes(colour = p_value, size = r_value), data = correlation) +  
  scale_size_manual(values = c(2,1.5,1,1.5)) +
  scale_colour_manual(values = c("#6B6EAE","#D6398A",'orange',"#A2A2A288")) +  ##调整线条颜色scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) 
  guides(size = guide_legend(title = "Spearman's correlation",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Spearman's pvalue", 
                               override.aes = list(size = 3), 
                               order = 2),
         fill = guide_colorbar(title = "Spearman's correlation", order = 3))
dev.off()

png('02.heatmap.png',w=700,h=700)

quickcor(corr,type = "lower") +
  geom_circle2()+
  scale_fill_gradient2(low = "#9483E1",mid = "white",high = "#E11953",midpoint=0) +
  geom_mark(sig.thres = 0.05, size = 2.5, colour = "black")+
  # geom_colour()+
  #geom_square() +
  anno_link(aes(colour = p_value, size = r_value), data = correlation) +  
  scale_size_manual(values = c(2,1.5,1,1.5)) +
  scale_colour_manual(values = c("#6B6EAE","#D6398A",'orange',"#A2A2A288")) +  ##调整线条颜色scale_colour_manual(values = c("#D95F02", "#1B9E77", "#A2A2A288")) 
  guides(size = guide_legend(title = "Spearman's correlation",
                             override.aes = list(colour = "grey35"), 
                             order = 2),
         colour = guide_legend(title = "Spearman's pvalue", 
                               override.aes = list(size = 3), 
                               order = 2),
         fill = guide_colorbar(title = "Spearman's correlation", order = 3))
dev.off()
