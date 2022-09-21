rm(list = ls())
### 验证集----------
setwd("/data/nas1/luchunlin/project/BJTC-258")
if (! dir.exists("./08_clinical")){
  dir.create("./08_clinical")
}
setwd("./08_clinical")

clinical<-read_tsv(file ='/data/nas1/luchunlin/project/BJTC-258/07_progmodel/TCGA-OV.GDC_phenotype.tsv')
phenotype<-data.frame(sample=clinical$submitter_id.samples,
                      stage=clinical$clinical_stage,
                      age=clinical$age_at_index.demographic)
score<-read.delim2('/data/nas1/luchunlin/project/BJTC-258/05_survival/clinical.xls')
train_phenotype<-merge(phenotype,score,by='sample')#%>%select(-'group')
train_phenotype$OS<-as.numeric(train_phenotype$OS)
train_phenotype$OS.time<-as.numeric(train_phenotype$OS.time)
train_phenotype2<-train_phenotype
table(train_phenotype2$stage)
train_phenotype2$stage<-gsub('A','',train_phenotype2$stage)
train_phenotype2$stage<-gsub('B','',train_phenotype2$stage)
train_phenotype2$stage<-gsub('C','',train_phenotype2$stage)
train_phenotype2$stage<-gsub('Stage IV','Stage4',train_phenotype2$stage)
train_phenotype2$stage<-gsub('Stage III','Stage3',train_phenotype2$stage)
train_phenotype2$stage<-gsub('Stage II','Stage1/2',train_phenotype2$stage)
train_phenotype2$stage<-gsub('Stage I','Stage1/2',train_phenotype2$stage)
table(train_phenotype2$age)
train_phenotype2$age<-ifelse(train_phenotype2$age>60,'>60','<=60')

library(ggpubr)
library(Ipaper)
library(ggthemes)
table(train_phenotype2$stage)
train_phenotype2$scores<-as.numeric(train_phenotype2$scores)
## 08-1 TNM stage-----
my_comparisons <- list(c("Stage1/2","Stage3"),c("Stage3","Stage4"),c("Stage1/2","Stage4"))
stage_data <- data.frame(AgingScore = train_phenotype2$scores,
                         stage = factor(train_phenotype2$stage,
                                        levels = c("Stage1/2","Stage3","Stage4")))
stage_data <- na.omit(stage_data)
stage<-ggplot(stage_data,aes(x = stage, y = AgingScore, fill = stage)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "AgingScore",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("Stage") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              test = wilcox.test,
              map_signif_level = T,
              y_position = c(3,4,5))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 20))+
  guides(fill='none')
stage
## AGE-------
table(train_phenotype2$age)
my_comparisons <- list(c(">60","<=60"))
age_data <- data.frame(AgingScore = train_phenotype2$scores,
                       age = factor(train_phenotype2$age,
                                    levels = c(">60","<=60")))
age_data <- na.omit(age_data)
age<-ggplot(age_data,aes(x = age, y = AgingScore, fill = age)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "AgingScore",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("Age") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              test = wilcox.test,
              map_signif_level = T,
              y_position = c(3,4,5))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 20))+
  guides(fill='none')
age

library(patchwork)
all_clinical_index <- age+stage+
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
ggsave('01.clinical.pdf',all_clinical_index,w=7,h=4)
ggsave('01.clinical.png',all_clinical_index,w=7,h=4)
