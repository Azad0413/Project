rm(list = ls())
#08 风险评分与临床指标相关性分析-------------
setwd("/data/nas1/luchunlin/project/YCZK-127/")
if (! dir.exists("./07_clinical_index")){
  dir.create("./07_clinical_index")
}
setwd("./07_clinical_index")
library(readr)
library(readxl)
survival<-read.delim2("/data/nas1/luchunlin/project/YCZK-127/03_univariate_cox/survival.xls")
phenotype<-read_tsv(file = 'TCGA-BRCA.GDC_phenotype.tsv.gz')
train_phenotype<-data.frame(sample=phenotype$submitter_id.samples,
                            TNM.stage=phenotype$tumor_stage.diagnoses,
                            T.stage=phenotype$pathologic_T,
                            N.stage=phenotype$pathologic_N,
                            M.stage=phenotype$pathologic_M)

train_phenotype<-merge(train_phenotype,survival,by='sample')
write.table(train_phenotype,file = 'phenotype.xls',sep = '\t',row.names = F,quote = F)
train_phenotype2<-train_phenotype
train_phenotype2$OS<-as.numeric(train_phenotype2$OS)
train_phenotype2$OS.time<-as.numeric(train_phenotype2$OS.time)
table(train_phenotype2$TNM.stage)
train_phenotype2$TNM.stage<-gsub('not reported',NA,train_phenotype2$TNM.stage,fixed = T)
train_phenotype2$TNM.stage<-gsub('ia','',train_phenotype2$TNM.stage,fixed = T)
train_phenotype2$TNM.stage<-gsub('ib','',train_phenotype2$TNM.stage,fixed = T)
train_phenotype2$TNM.stage<-gsub('ic','',train_phenotype2$TNM.stage,fixed = T)
train_phenotype2$TNM.stage<-gsub('stage x',NA,train_phenotype2$TNM.stage,fixed = T)
train_phenotype2$TNM.stage<-gsub('stage iv','Stage 3/4',train_phenotype2$TNM.stage,fixed = T)
train_phenotype2$TNM.stage<-gsub('stage iii','Stage 3/4',train_phenotype2$TNM.stage,fixed = T)
train_phenotype2$TNM.stage<-gsub('stage ii','Stage 2',train_phenotype2$TNM.stage,fixed = T)
train_phenotype2$TNM.stage<-gsub('stage i','Stage 1',train_phenotype2$TNM.stage,fixed = T)
train_phenotype2$TNM.stage<-gsub('stage ','Stage 1',train_phenotype2$TNM.stage,fixed = T)

table(train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('a','',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('b','',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('c','',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('d','',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('TX',NA,train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('T3','T3/4',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('T4','T3/4',train_phenotype2$T.stage)

table(train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('a','',train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('b','',train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('c','',train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub(' (i-)','',train_phenotype2$N.stage,fixed = T)
train_phenotype2$N.stage<-gsub(' (i+)','',train_phenotype2$N.stage,fixed = T)
train_phenotype2$N.stage<-gsub(' (mol+)','',train_phenotype2$N.stage,fixed = T)
train_phenotype2$N.stage<-gsub('mi','',train_phenotype2$N.stage,fixed = T)
train_phenotype2$N.stage<-gsub('NX',NA,train_phenotype2$N.stage,fixed = T)
train_phenotype2$N.stage<-gsub('N2','N2/3',train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('N3','N2/3',train_phenotype2$N.stage)


table(train_phenotype2$M.stage)
train_phenotype2$M.stage<-gsub('cM0 (i+)','M0',train_phenotype2$M.stage,fixed = T)
train_phenotype2$M.stage<-gsub('MX',NA,train_phenotype2$M.stage,fixed = T)
library(tidyverse)
library(lance)
colnames(train_phenotype2)<-c('id','Stage','T.stage','N.stage','M.stage','OS','OS.time')
risk<-read.delim2('/data/nas1/luchunlin/project/YCZK-127/05_risk/risk.xls')%>%lc.tableToNum()
sub_risk <- subset(risk, select = c(id, riskScore))
train_phenotype3 <- merge(sub_risk,
                          train_phenotype2,
                          by = "id")
train_phenotype3$Group<-ifelse(train_phenotype3$riskScore>median(train_phenotype3$riskScore),'High risk','Low risk')

write.table(train_phenotype3,file = "clinical_risk.xls",row.names = F,sep = "\t",quote = F)

library(ggpubr)
library(Ipaper)
library(ggthemes)
# table(train_phenotype3$Stage)
# ## 08-1 TNM stage-----
# my_comparisons <- list(c("Stage 1","Stage 2"),c("Stage 2","Stage 3/4"),c("Stage 1","Stage 3/4"))
# stage_data <- data.frame(riskScore = train_phenotype3$riskScore,
#                          stage = factor(train_phenotype3$Stage,
#                                         levels = c("Stage 1","Stage 2", "Stage 3/4")))
# stage_data <- na.omit(stage_data)
# stage<-ggplot(stage_data,aes(x = stage, y = riskScore, fill = stage)) +
#   geom_boxplot(alpha=0.7) +
#   scale_y_continuous(name = "riskScore",expand = c(0.1,0.1))+
#   scale_x_discrete(name = "") +
#   ggtitle("Stage") +
#   theme_bw() +
#   geom_signif(comparisons = my_comparisons,
#               test = wilcox.test,
#               map_signif_level = T,
#               y_position = c(2,3,4))+
#   theme(plot.title = element_text(size = 16, face =  "bold"),
#         text = element_text(size = 20),
#         axis.title = element_text(face="bold"),
#         axis.text.x=element_text(size = 20))+
#   guides(fill='none')
# stage

## 08-2 T stage-----
colnames(train_phenotype3)
train_phenotype3<-na.omit(train_phenotype3)
exp_plot <- ggplot(train_phenotype3,aes(x = `T.stage`, y = riskScore, fill = Group)) +
  #geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#CD3700","#4682B4"), name = "Group")+
  labs(title="", x="", y = "",size=20) +
  stat_compare_means(data = train_phenotype3,
                     mapping = aes(group = Group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F) +
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18),
        axis.text.x=element_text(angle=0,hjust=0.5,colour="black",face="bold",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+facet_wrap(~`T.stage`,scales = "free",nrow = 1) 
exp_plot

exp_plot2 <- ggplot(train_phenotype3,aes(x = `N.stage`, y = riskScore, fill = Group)) +
  #geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#CD3700","#4682B4"), name = "Group")+
  labs(title="", x="", y = "",size=20) +
  stat_compare_means(data = train_phenotype3,
                     mapping = aes(group = Group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F) +
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18),
        axis.text.x=element_text(angle=0,hjust=0.5,colour="black",face="bold",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill='none')+
  facet_wrap(~`N.stage`,scales = "free",nrow = 1) 

exp_plot2

exp_plot3 <- ggplot(train_phenotype3,aes(x = `M.stage`, y = riskScore, fill = Group)) +
  #geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#CD3700","#4682B4"), name = "Group")+
  labs(title="", x="", y = "",size=20) +
  stat_compare_means(data = train_phenotype3,
                     mapping = aes(group = Group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F) +
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18),
        axis.text.x=element_text(angle=0,hjust=0.5,colour="black",face="bold",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
#  guides(fill='none')+
  facet_wrap(~`M.stage`,scales = "free",nrow = 1) 

exp_plot3

# table(train_phenotype3$T.stage)
# my_comparisons <- list(c("T1","T2"),c("T2","T3/4"),c("T1","T3/4"))
# T.stage_data <- data.frame(riskScore = train_phenotype3$riskScore,
#                          stage = factor(train_phenotype3$T.stage,
#                                         levels = c("T1","T2", "T3/4")))
# T.stage_data <- na.omit(T.stage_data)
# T.stage<-ggplot(T.stage_data,aes(x = stage, y = riskScore, fill = stage)) +
#   geom_boxplot(alpha=0.7) +
#   scale_y_continuous(name = "riskScore",expand = c(0.1,0.1))+
#   scale_x_discrete(name = "") +
#   ggtitle("T Stage") +
#   theme_bw() +
#   geom_signif(comparisons = my_comparisons,
#               test = wilcox.test,
#               map_signif_level = T,
#               y_position = c(2,3,4))+
#   theme(plot.title = element_text(size = 16, face =  "bold"),
#         text = element_text(size = 20),
#         axis.title = element_text(face="bold"),
#         axis.text.x=element_text(size = 20))+
#   guides(fill='none')
# T.stage
# 
# ## 08-3 N stage-----
# table(train_phenotype3$N.stage)
# my_comparisons <- list(c("N0","N1"),c("N1","N2/3"),c("N0","N2/3"))
# N.stage_data <- data.frame(riskScore = train_phenotype3$riskScore,
#                            stage = factor(train_phenotype3$N.stage,
#                                           levels = c("N0","N1", "N2/3")))
# N.stage_data <- na.omit(N.stage_data)
# N.stage<-ggplot(N.stage_data,aes(x = stage, y = riskScore, fill = stage)) +
#   geom_boxplot(alpha=0.7) +
#   scale_y_continuous(name = "riskScore",expand = c(0.1,0.1))+
#   scale_x_discrete(name = "") +
#   ggtitle("N Stage") +
#   theme_bw() +
#   geom_signif(comparisons = my_comparisons,
#               test = wilcox.test,
#               map_signif_level = T,
#               y_position = c(2,3,4))+
#   theme(plot.title = element_text(size = 16, face =  "bold"),
#         text = element_text(size = 20),
#         axis.title = element_text(face="bold"),
#         axis.text.x=element_text(size = 20))+
#   guides(fill='none')
# N.stage
# 
# ## 08-3 N stage-----
# table(train_phenotype3$M.stage)
# my_comparisons <- list(c("M0","M1"))
# M.stage_data <- data.frame(riskScore = train_phenotype3$riskScore,
#                            stage = factor(train_phenotype3$M.stage,
#                                           levels = c("M0","M1")))
# M.stage_data <- na.omit(M.stage_data)
# M.stage<-ggplot(M.stage_data,aes(x = stage, y = riskScore, fill = stage)) +
#   geom_boxplot(alpha=0.7) +
#   scale_y_continuous(name = "riskScore",expand = c(0.1,0.1))+
#   scale_x_discrete(name = "") +
#   ggtitle("M Stage") +
#   theme_bw() +
#   geom_signif(comparisons = my_comparisons,
#               test = wilcox.test,
#               map_signif_level = T,
#               y_position = c(2))+
#   theme(plot.title = element_text(size = 16, face =  "bold"),
#         text = element_text(size = 20),
#         axis.title = element_text(face="bold"),
#         axis.text.x=element_text(size = 20))+
#   guides(fill='none')
# M.stage

# library(patchwork)
# all_clinical_index <- stage + T.stage + N.stage + M.stage+ 
#   plot_layout(ncol = 2) & 
#   theme_bw() & 
#   theme(legend.title = element_blank(),
#         legend.position = "top",
#         legend.text = element_text(face = "bold"),
#         plot.title = element_text(hjust = 0.5, face = "bold"),
#         panel.grid = element_blank(),
#         axis.title = element_text(face = "bold"),
#         axis.text = element_text(face = "bold",size = 12))
# all_clinical_index
