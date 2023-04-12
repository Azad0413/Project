rm(list = ls())
#08 风险评分与临床指标相关性分析-------------
setwd("/data/nas1/luchunlin/project/LZZK-512/")
if (! dir.exists("./08_clinical_index")){
  dir.create("./08_clinical_index")
}
setwd("./08_clinical_index")
library(readr)
library(readxl)
survival<-read.delim2("/data/nas1/luchunlin/project/LZZK-512/03_univariate_cox/survival.xls")
phenotype<-read_tsv(file = 'TCGA-BLCA.GDC_phenotype.tsv.gz')

train_phenotype<-data.frame(sample=phenotype$submitter_id.samples,
                            age=phenotype$age_at_initial_pathologic_diagnosis,
                            gender=phenotype$gender.demographic,
                            stage=phenotype$tumor_stage.diagnoses,
                            T.stage=phenotype$clinical_T
                            )
train_phenotype<-merge(train_phenotype,survival,by='sample')
write.table(train_phenotype,file = 'phenotype.xls',sep = '\t',row.names = F,quote = F)
train_phenotype2<-train_phenotype
train_phenotype2$OS<-as.numeric(train_phenotype2$OS)
train_phenotype2$OS.time<-as.numeric(train_phenotype2$OS.time)
table(train_phenotype2$stage)
train_phenotype2$stage<-gsub('not reported',NA,train_phenotype2$stage,fixed = T)
train_phenotype2$stage<-gsub('stage iv','Stage 3/4',train_phenotype2$stage,fixed = T)
train_phenotype2$stage<-gsub('stage iii','Stage 3/4',train_phenotype2$stage,fixed = T)
train_phenotype2$stage<-gsub('stage ii','Stage 1/2',train_phenotype2$stage,fixed = T)
train_phenotype2$stage<-gsub('stage i','Stage 1/2',train_phenotype2$stage,fixed = T)

table(train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('a','',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('b','',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('TX',NA,train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('T3','T3/4',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('T4','T3/4',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('T1','T1/2',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('T2','T1/2',train_phenotype2$T.stage)

table(train_phenotype2$age)
train_phenotype2$age<-ifelse(train_phenotype2$age>60,'>60','≤60')

library(tidyverse)
library(lance)
colnames(train_phenotype2)
colnames(train_phenotype2)<-c('id','Age','Gender','Stage','T.stage','OS','OS.time')
risk<-read.delim2('/data/nas1/luchunlin/project/LZZK-512/06_risk/risk.xls')%>%lc.tableToNum()
sub_risk <- subset(risk, select = c(id, riskScore))
train_phenotype3 <- merge(sub_risk,
                          train_phenotype2,
                          by = "id")
train_phenotype3$Group<-ifelse(train_phenotype3$riskScore>median(train_phenotype3$riskScore),'High risk','Low risk')

write.table(train_phenotype3,file = "clinical_risk.xls",row.names = F,sep = "\t",quote = F)

library(ggpubr)
library(Ipaper)
library(ggthemes)
table(train_phenotype3$Stage)
## 08-1 TNM stage-----
my_comparisons <- list(c("Stage 1/2","Stage 3/4"))
stage_data <- data.frame(riskScore = train_phenotype3$riskScore,
                         stage = factor(train_phenotype3$Stage,
                                        levels = c("Stage 1/2","Stage 3/4")))
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
              y_position = c(3))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 20))+
  guides(fill='none')
stage

## 08-2 T stage-----

table(train_phenotype3$T.stage)
my_comparisons <- list(c("T1/2","T3/4"))
T.stage_data <- data.frame(riskScore = train_phenotype3$riskScore,
                         stage = factor(train_phenotype3$T.stage,
                                        levels = c("T1/2","T3/4")))
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
              y_position = c(2,3,4))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 20))+
  guides(fill='none')
T.stage
# 
## 08-3 age-----
table(train_phenotype3$Age)
my_comparisons <- list(c(">60","≤60"))
age_data <- data.frame(riskScore = train_phenotype3$riskScore,
                           age = factor(train_phenotype3$Age,
                                          levels = c(">60","≤60")))
age_data <- na.omit(age_data)
age<-ggplot(age_data,aes(x = age, y = riskScore, fill = age)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "riskScore",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("Age") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              test = wilcox.test,
              map_signif_level = T,
              y_position = c(3))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 20))+
  guides(fill='none')
age
# 
# ## 08-3gender-----
table(train_phenotype3$Gender)
my_comparisons <- list(c("female","male"))
gender_data <- data.frame(riskScore = train_phenotype3$riskScore,
                           gender = factor(train_phenotype3$Gender,
                                          levels = c("female","male")))
gender_data <- na.omit(gender_data)
gender<-ggplot(gender_data,aes(x = gender, y = riskScore, fill = gender)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "riskScore",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("Gender") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              test = wilcox.test,
              map_signif_level = T,
              y_position = c(3))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 20))+
  guides(fill='none')
gender

library(patchwork)
all_clinical_index <- stage + T.stage + age + gender+
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
ggsave('clinical.pdf',all_clinical_index,w=6,h=6)
ggsave('clinical.png',all_clinical_index,w=6,h=6)
