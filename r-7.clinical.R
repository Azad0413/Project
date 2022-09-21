rm(list = ls())
#07 风险评分与临床指标相关性分析-------------
setwd("/data/nas1/luchunlin/project/BJTC-357/")
if (! dir.exists("./07_clinical")){
  dir.create("./07_clinical")
}
setwd("./07_clinical")
library(readr)
library(readxl)
survival<-read.delim2('/data/nas1/luchunlin/project/BJTC-357/00_rawdata/survival.xls')
phenotype<-read.delim2('/data/nas1/luchunlin/project/BJTC-357/00_rawdata/phenotype.xls')
train.data<-read.delim2('../03_univariate_cox/train.data.xls',row.names = 1)
phenotype<-phenotype[phenotype$sample%in%rownames(train.data),]
train_phenotype<-merge(phenotype,survival,by='sample')
table(train_phenotype$gender)
train_phenotype$gender<-factor(ifelse(train_phenotype$gender=='Sex: female','female','male'))
write.table(train_phenotype,file = 'phenotype.xls',sep = '\t',row.names = F,quote = F)
train_phenotype2<-train_phenotype
train_phenotype2$OS<-as.numeric(train_phenotype2$OS)
train_phenotype2$OS.time<-as.numeric(train_phenotype2$OS.time)
table(train_phenotype2$TNM.stage)
train_phenotype2$TNM.stage<-gsub('tnm stage: III','Stage 3',train_phenotype2$TNM.stage,fixed = T)
train_phenotype2$TNM.stage<-gsub('tnm stage: II','Stage 2',train_phenotype2$TNM.stage,fixed = T)
train_phenotype2$TNM.stage<-gsub('tnm stage: I','Stage 1',train_phenotype2$TNM.stage,fixed = T)

table(train_phenotype2$T.stage)

train_phenotype2$T.stage<-gsub('t stage: T1','T1/2',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('t stage: T2','T1/2',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('t stage: T3','T3',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('t stage: T4','T4',train_phenotype2$T.stage)

table(train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('n stage: N0','N0',train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('n stage: N1','N1',train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('n stage: N2','N2/N3',train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('n stage: N3','N2/N3',train_phenotype2$N.stage)

library(tidyverse)
library(lance)
colnames(train_phenotype2)
colnames(train_phenotype2)<-c('id','Gender','T.stage','N.stage','TNM.stage','OS','OS.time')
risk<-read.delim2('/data/nas1/luchunlin/project/BJTC-357/05_risk/risk.xls')%>%lc.tableToNum()
sub_risk <- subset(risk, select = c(id, riskScore))
train_phenotype3 <- merge(sub_risk,
                          train_phenotype2,
                          by = "id")
#train_phenotype3$Group<-ifelse(train_phenotype3$riskScore>median(train_phenotype3$riskScore),'High risk','Low risk')

write.table(train_phenotype3,file = "clinical_risk.xls",row.names = F,sep = "\t",quote = F)

library(ggpubr)
library(Ipaper)
library(ggthemes)

## 07-1 TNM stage-----
table(train_phenotype3$TNM.stage)
my_comparisons <- list(c("Stage 1","Stage 2"),c("Stage 2","Stage 3"),c("Stage 1","Stage 3"))
stage_data <- data.frame(riskScore = train_phenotype3$riskScore,
                         stage = factor(train_phenotype3$TNM.stage,
                                        levels = c("Stage 1","Stage 2", "Stage 3")))
stage_data <- na.omit(stage_data)
stage<-ggplot(stage_data,aes(x = stage, y = riskScore, fill = stage)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "riskScore",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("Stage") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              test = wilcox.test,
              map_signif_level = T,
              y_position = c(17,18,19))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 20))+
  guides(fill='none')
stage

## 07-2 T stage-----

table(train_phenotype3$T.stage)
my_comparisons <- list(c("T1/2","T3"),c("T1/2","T4"),c("T3","T4"))
T.stage_data <- data.frame(riskScore = train_phenotype3$riskScore,
                           stage = factor(train_phenotype3$T.stage,
                                          levels = c("T1/2","T3","T4")))
T.stage_data <- na.omit(T.stage_data)
T.stage<-ggplot(T.stage_data,aes(x = stage, y = riskScore, fill = stage)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "riskScore",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("T Stage") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              test = wilcox.test,
              map_signif_level = T,
              y_position = c(17,19,18))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 20))+
  guides(fill='none')
T.stage

## 07-3 N stage-----

table(train_phenotype3$N.stage)

my_comparisons <- list(c("N0","N1"),c('N0','N2/N3'),c("N1",'N2/N3'))
N.stage_data <- data.frame(riskScore = train_phenotype3$riskScore,
                           stage = factor(train_phenotype3$N.stage,
                                          levels = c("N0","N1",'N2/N3')))
N.stage_data <- na.omit(N.stage_data)
N.stage<-ggplot(N.stage_data,aes(x = stage, y = riskScore, fill = stage)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "riskScore",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("N Stage") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              test = wilcox.test,
              map_signif_level = T,
              y_position = c(17,19,18))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 20))+
  guides(fill='none')
N.stage
# ## gender----
table(train_phenotype3$Gender)
my_comparisons <- list(c("female","male"))
gender_data <- data.frame(riskScore = train_phenotype3$riskScore,
                          stage = factor(train_phenotype3$Gender,
                                         levels = c("female","male")))
gender_data <- na.omit(gender_data)
gender<-ggplot(gender_data,aes(x = stage, y = riskScore, fill = stage)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "riskScore",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("Gender") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              test = wilcox.test,
              map_signif_level = T,
              y_position = c(17))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 20))+
  guides(fill='none')
gender

library(patchwork)
all_clinical_index <- stage + T.stage + N.stage +gender+
  plot_layout(ncol = 2) &
  theme_bw() &
  theme(legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold",size = 12))
all_clinical_index
ggsave(filename = 'clinical.pdf',all_clinical_index,w=8,h=7.5)
ggsave(filename = 'clinical.png',all_clinical_index,w=8,h=7.5)
