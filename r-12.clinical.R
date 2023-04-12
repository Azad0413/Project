rm(list = ls())
#08 风险评分与临床指标相关性分析-------------
setwd("/data/nas1/luchunlin/project/BJTC-406-12/")
if (! dir.exists("./12_clinical")){
  dir.create("./12_clinical")
}
setwd("./12_clinical")
library(readr)
library(readxl)
survival<-read.delim2("../07_univariate_cox/survival.xls")
phenotype<-read_tsv(file = '/data/nas1/luchunlin/TCGA_phenotype/TCGA-LUAD.GDC_phenotype.tsv.gz')
# age、gender、stage、smoking、T stage、N stage、M stage
train_phenotype<-data.frame(sample=phenotype$submitter_id.samples,
                            age=phenotype$age_at_initial_pathologic_diagnosis,
                            gender=phenotype$gender.demographic,
                            #     smoking=phenotype$tobacco_smoking_history,
                            # stage=phenotype$tumor_stage.diagnoses,
                            T.stage=phenotype$pathologic_T,
                            N.stage=phenotype$pathologic_N,
                            M.stage=phenotype$pathologic_M)

train_phenotype<-merge(train_phenotype,survival,by='sample')
write.table(train_phenotype,file = 'phenotype.xls',sep = '\t',row.names = F,quote = F)
train_phenotype2<-train_phenotype
train_phenotype2$OS<-as.numeric(train_phenotype2$OS)
train_phenotype2$OS.time<-as.numeric(train_phenotype2$OS.time)

table(train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('a','',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('b','',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('TX',NA,train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('T3','T3/4',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('T4','T3/4',train_phenotype2$T.stage)

table(train_phenotype2$N.stage)

train_phenotype2$N.stage<-gsub('NX',NA,train_phenotype2$N.stage,fixed = T)
train_phenotype2$N.stage<-gsub('N2','N2/3',train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('N3','N2/3',train_phenotype2$N.stage)


table(train_phenotype2$M.stage)
train_phenotype2$M.stage<-gsub('MX',NA,train_phenotype2$M.stage,fixed = T)
train_phenotype2$M.stage<-gsub('a','',train_phenotype2$M.stage,fixed = T)
train_phenotype2$M.stage<-gsub('b','',train_phenotype2$M.stage,fixed = T)

table(train_phenotype2$age)

train_phenotype2$age <- ifelse(train_phenotype2$age>65,">65","<=65")

table(train_phenotype2$gender)

library(tidyverse)
library(lance)
colnames(train_phenotype2)
colnames(train_phenotype2)<-c('id','Age','Gender','T.stage','N.stage','M.stage','OS','OS.time')
risk<-read.delim2('../09_risk/risk.xls')%>%lc.tableToNum()
sub_risk <- subset(risk, select = c(id, riskScore))
train_phenotype3 <- merge(sub_risk,
                          train_phenotype2,
                          by = "id")
train_phenotype3$Group<-ifelse(train_phenotype3$riskScore>median(train_phenotype3$riskScore),'High risk','Low risk')

write.table(train_phenotype3,file = "clinical_risk.xls",row.names = F,sep = "\t",quote = F)

library(ggpubr)
library(Ipaper)
library(ggthemes)
##Age-----
table(train_phenotype3$Age)
my_comparisons <- list(c("<=65",">65"))
age_data <- data.frame(riskScore = train_phenotype3$riskScore,
                       stage = factor(train_phenotype3$Age,
                                      levels = c("<=65",">65")))
age_data <- na.omit(age_data)
age<-ggplot(age_data,aes(x = stage, y = riskScore, fill = stage)) +
  geom_violin() + #绘制小提琴图
  stat_boxplot(geom="errorbar",
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.4,
               position=position_dodge(0.9),
               outlier.shape = NA,
               fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  geom_point(aes(fill = stage),
             size = 0.5,
             position = position_dodge(0.9))+
  scale_y_continuous(name = "riskScore",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("Age") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              test = wilcox.test,
              map_signif_level = T,
              y_position = c(2.5))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 20))+
  guides(fill='none')
age

##gender-----
table(train_phenotype3$Gender)
my_comparisons <- list(c("female","male"))
gender_data <- data.frame(riskScore = train_phenotype3$riskScore,
                       gender = factor(train_phenotype3$Gender,
                                      levels = c("female","male")))
gender_data <- na.omit(gender_data)
gender<-ggplot(gender_data,aes(x = gender, y = riskScore, fill = gender)) +
  geom_violin() + #绘制小提琴图
  stat_boxplot(geom="errorbar",
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.4,
               position=position_dodge(0.9),
               outlier.shape = NA,
               fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  geom_point(aes(fill = gender),
             size = 0.5,
             position = position_dodge(0.9))+
  scale_y_continuous(name = "riskScore",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("Gender") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              test = wilcox.test,
              map_signif_level = T,
              y_position = c(2.5))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 20))+
  guides(fill='none')
gender

## T.stage-----
table(train_phenotype3$T.stage)
my_comparisons <- list(c("T1","T2"),c("T2","T3/4"),c("T1","T3/4"))
T.stage_data <- data.frame(riskScore = train_phenotype3$riskScore,
                           stage = factor(train_phenotype3$T.stage,
                                          levels = c("T1","T2", "T3/4")))
T.stage_data <- na.omit(T.stage_data)
T.stage<-ggplot(T.stage_data,aes(x = stage, y = riskScore, fill = stage)) +
  geom_violin() + #绘制小提琴图
  stat_boxplot(geom="errorbar",
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.4,
               position=position_dodge(0.9),
               outlier.shape = NA,
               fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  geom_point(aes(fill = stage),
             size = 0.5,
             position = position_dodge(0.9))+
  scale_y_continuous(name = "riskScore",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("T Stage") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              test = wilcox.test,
              map_signif_level = T,
              y_position = c(2.5,3,3.5))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 20))+
  guides(fill='none')
T.stage

## 08-3 N stage-----
table(train_phenotype3$N.stage)
my_comparisons <- list(c("N0","N1"),c("N1","N2/3"),c("N0","N2/3"))
N.stage_data <- data.frame(riskScore = train_phenotype3$riskScore,
                           stage = factor(train_phenotype3$N.stage,
                                          levels = c("N0","N1", "N2/3")))
N.stage_data <- na.omit(N.stage_data)
N.stage<-ggplot(N.stage_data,aes(x = stage, y = riskScore, fill = stage)) +
  geom_violin() + #绘制小提琴图
  stat_boxplot(geom="errorbar",
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.4,
               position=position_dodge(0.9),
               outlier.shape = NA,
               fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  geom_point(aes(fill = stage),
             size = 0.5,
             position = position_dodge(0.9))+
  scale_y_continuous(name = "riskScore",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("N Stage") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              test = wilcox.test,
              map_signif_level = T,
              y_position = c(2.3,3,3.5))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 20))+
  guides(fill='none')
N.stage
#
# ## 08-3 M stage-----
table(train_phenotype3$M.stage)
my_comparisons <- list(c("M0","M1"))
M.stage_data <- data.frame(riskScore = train_phenotype3$riskScore,
                           stage = factor(train_phenotype3$M.stage,
                                          levels = c("M0","M1")))
M.stage_data <- na.omit(M.stage_data)
M.stage<-ggplot(M.stage_data,aes(x = stage, y = riskScore, fill = stage)) +
  geom_violin() + #绘制小提琴图
  stat_boxplot(geom="errorbar",
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.4,
               position=position_dodge(0.9),
               outlier.shape = NA,
               fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  geom_point(aes(fill = stage),
             size = 0.5,
             position = position_dodge(0.9))+
  scale_y_continuous(name = "riskScore",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("M Stage") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              test = wilcox.test,
              map_signif_level = T,
              y_position = c(2.5))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 20))+
  guides(fill='none')
M.stage

library(patchwork)
all_clinical_index <- age + gender + T.stage + N.stage + M.stage+
  plot_layout(ncol = 3) &
  theme_bw() &
  theme(legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold",size = 12))
all_clinical_index
ggsave(all_clinical_index,filename = '01.clinical.pdf',w=7,h=6)
ggsave(all_clinical_index,filename = '01.clinical.png',w=7,h=6)


