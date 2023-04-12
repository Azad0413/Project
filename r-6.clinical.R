rm(list = ls())
#08 风险评分与临床指标相关性分析-------------
setwd("/data/nas1/luchunlin/project/JNZK-204(modify)//")
if (! dir.exists("./06_clinical")){
  dir.create("./06_clinical")
}
setwd("./06_clinical")
library(readr)
library(readxl)
survival<-read.delim2("/data/nas1/yangly/Project/JNZK204/clean.data/survival.tsv")
survival <- survival[,-3]
survival$sample <- gsub('.','-',survival$sample,fixed = T)
phenotype<-read.delim2('../pheno.tsv')
train_phenotype<-data.frame(sample=phenotype$submitter_id.samples,
                            T.stage=phenotype$pathologic_T,
                            Grade=phenotype$neoplasm_histologic_grade,
                            Gender=phenotype$gender.demographic,
                            ATI=phenotype$adjacent_hepatic_tissue_inflammation_extent_type)

train_phenotype<-merge(train_phenotype,survival,by='sample')
write.table(train_phenotype,file = 'phenotype.xls',sep = '\t',row.names = F,quote = F)
train_phenotype2<-train_phenotype
train_phenotype2$OS<-as.numeric(train_phenotype2$OS)
train_phenotype2$OS.time<-as.numeric(train_phenotype2$OS.time)

table(train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('a','',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('b','',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('TX',NA,train_phenotype2$T.stage)
train_phenotype2[train_phenotype2==''] <- NA

table(train_phenotype2$Grade)
table(train_phenotype2$ATI)

table(train_phenotype2$Gender)

library(tidyverse)
library(lance)
colnames(train_phenotype2)
colnames(train_phenotype2)<-c('id','T.stage','Grade','Gender','ATI','OS','OS.time')
risk<-read.delim2('../06_clinical/risk.xls')%>%lc.tableToNum()
sub_risk <- subset(risk, select = c(sample, risk))
colnames(sub_risk) <- c('id','riskScore')
sub_risk$id <- gsub('.','-',sub_risk$id,fixed = T)
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
exp_plot <- ggplot(train_phenotype3,aes(x = Group, y = riskScore, fill = Group)) +
  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA,fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#CD3700","#4682B4"), name = "Group")+
  labs(title="", x="", y = "",size=20) +
  stat_compare_means(data = train_phenotype3,
                     mapping = aes(group = Group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F,label.y = 1.5,label.x = 1.3) +
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
        panel.grid.minor = element_blank())+facet_wrap(~`T.stage`,scales = "free",nrow = 1) +
guides(fill='none')
exp_plot
ggsave(exp_plot,filename = '01.Tstage.pdf',w=6,h=4)
ggsave(exp_plot,filename = '01.Tstage.png',w=6,h=4)

exp_plot2 <- ggplot(train_phenotype3,aes(x = Group, y = riskScore, fill = Group)) +
  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA,fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#CD3700","#4682B4"), name = "Group")+
  labs(title="", x="", y = "",size=20) +
  stat_compare_means(data = train_phenotype3,
                     mapping = aes(group = Group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F,label.y = 1.5,label.x = 1.4) +
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
  facet_wrap(~Grade,scales = "free",nrow = 1) 

exp_plot2

ggsave(exp_plot2,filename = '02.Grade.pdf',w=6,h=4)
ggsave(exp_plot2,filename = '02.Grade.png',w=6,h=4)

exp_plot3 <- ggplot(train_phenotype3,aes(x = Group, y = riskScore, fill = Group)) +
  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA,fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#CD3700","#4682B4"), name = "Group")+
  labs(title="", x="", y = "",size=20) +
  stat_compare_means(data = train_phenotype3,
                     mapping = aes(group = Group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F,label.y = 1,label.x = 1.45) +
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
  facet_wrap(~Gender,scales = "free",nrow = 1) 

exp_plot3
ggsave(exp_plot3,filename = '03.Gender.pdf',w=6,h=4)
ggsave(exp_plot3,filename = '03.Gender.png',w=6,h=4)



ati <- na.omit(train_phenotype3)
exp_plot4 <- ggplot(ati,aes(x = Group, y = riskScore, fill = Group)) +
  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA,fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#CD3700","#4682B4"), name = "Group")+
  labs(title="", x="", y = "",size=20) +
  stat_compare_means(data = ati,
                     mapping = aes(group = Group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F,label.y = 1,label.x = 1.45) +
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
  facet_wrap(~ATI,scales = "free",nrow = 1) 

exp_plot4
ggsave(exp_plot4,filename = '04.Adjacent-Tissue-Inflammation.pdf',w=6,h=4)
ggsave(exp_plot4,filename = '04.Adjacent-Tissue-Inflammation.png',w=6,h=4)

library(patchwork)
all <- exp_plot+exp_plot2+exp_plot3+exp_plot4&
  theme(axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=10))

all
ggsave(all,filename = '05.all.clinical.pdf',w=9,h=8)
ggsave(all,filename = '05.all.clinical.png',w=9,h=8)
