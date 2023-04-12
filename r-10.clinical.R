rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-399-11/")
if (! dir.exists("./10_clinical")){
  dir.create("./10_clinical")
}
setwd("./10_clinical")
library(readr)
library(readxl)
survival<-read.delim2('../06_univariate_cox/survival.xls')
phenotype<-read.delim2('/data/nas1/luchunlin/TCGA_phenotype/TCGA-GBM.GDC_phenotype.tsv.gz')
#年龄、性别
train_phenotype<-data.frame(sample=phenotype$submitter_id.samples,
                            age=phenotype$age_at_initial_pathologic_diagnosis,
                            gender=phenotype$gender.demographic)
train_phenotype<-merge(train_phenotype,survival,by='sample')
write.table(train_phenotype,file = 'phenotype.xls',sep = '\t',row.names = F,quote = F)
train_phenotype2<-train_phenotype
train_phenotype2$OS<-as.numeric(train_phenotype2$OS)
train_phenotype2$OS.time<-as.numeric(train_phenotype2$OS.time)

table(train_phenotype2$age)
median(train_phenotype2$age)
train_phenotype2$age<-ifelse(train_phenotype2$age>50,'>50','<=50')

table(train_phenotype2$gender)


library(tidyverse)
library(lance)
colnames(train_phenotype2)
colnames(train_phenotype2)<-c('id','Age','Gender','OS','OS.time')
risk<-read.delim2('../08_risk/risk.xls')%>%lc.tableToNum()
sub_risk <- subset(risk, select = c(id, riskScore))
train_phenotype3 <- merge(sub_risk,
                          train_phenotype2,
                          by = "id")
train_phenotype3$Group<-ifelse(train_phenotype3$riskScore>median(train_phenotype3$riskScore),'High risk','Low risk')

write.table(train_phenotype3,file = "clinical_risk.xls",row.names = F,sep = "\t",quote = F)

library(ggpubr)
library(Ipaper)
library(ggthemes)
## 08-3 age-----
table(train_phenotype3$Age)
my_comparisons <- list(c(">50","<=50"))
age_data <- data.frame(riskScore = train_phenotype3$riskScore,
                       age = factor(train_phenotype3$Age,
                                    levels = c(">50","<=50")))
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
              y_position = c(5))+
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
              y_position = c(5.5))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 20))+
  guides(fill='none')
gender



library(patchwork)
all_clinical_index <- age + gender+
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
ggsave('01.clinical.pdf',all_clinical_index,w=5,h=3)
ggsave('01.clinical.png',all_clinical_index,w=5,h=3)

