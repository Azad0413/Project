rm(list = ls())
#08 风险评分与临床指标相关性分析-------------
setwd("/data/nas1/luchunlin/project/JNZK-207/")
if (! dir.exists("./07_clinical_index")){
  dir.create("./07_clinical_index")
}
setwd("./07_clinical_index")
survival<-read.delim2("/data/nas1/luchunlin/project/JNZK-207/03_univariate_cox/survival.xls")
phenotype<-read_tsv(file = 'MMRF-COMMPASS.Xena_phenotype.tsv')
phenotype <- phenotype[phenotype$samples.submitter_id%in%survival$sample,]
write.table(phenotype,file = 'clinical(MMRF-COMMPASS).xls',sep = '\t',row.names = F,quote = F)

phenotype<-data.frame(sample=phenotype$samples.submitter_id,
                      gender=phenotype$demographic.gender,
                      race=phenotype$demographic.race,
                      stage=phenotype$diagnoses.tumor_stage)

phenotype<-merge(survival,phenotype,by='sample')


write.table(phenotype,file = 'phenotype.xls',sep = '\t',row.names = F,quote = F)
train_phenotype2<-phenotype
train_phenotype2$OS<-as.numeric(train_phenotype2$OS)
train_phenotype2$OS.time<-as.numeric(train_phenotype2$OS.time)
train_phenotype2$stage<-gsub('III','Stage3',train_phenotype2$stage)
train_phenotype2$stage<-gsub('II','Stage2',train_phenotype2$stage)
train_phenotype2$stage<-gsub('I','Stage1',train_phenotype2$stage)
train_phenotype2$stage<-gsub('Unknown',NA,train_phenotype2$stage)
train_phenotype2$stage<-gsub('unknown',NA,train_phenotype2$stage)
train_phenotype2$race<-gsub('not reported',NA,train_phenotype2$race)
train_phenotype2$race<-gsub('NA',NA,train_phenotype2$race)
train_phenotype2$race<-gsub('not allowed to collect',NA,train_phenotype2$race)
train_phenotype2$race<-gsub('other',NA,train_phenotype2$race)

colnames(train_phenotype2)<-c('id','OS','OS.time','gender','race','stage')
risk<-read.delim2('/data/nas1/luchunlin/project/JNZK-207/05_risk/risk.xls')%>%lc.tableToNum()
sub_risk <- subset(risk, select = c(id, riskScore))
train_phenotype3 <- merge(sub_risk,
                          train_phenotype2,
                          by = "id")
write.table(train_phenotype3,file = "clinical_risk.xls",row.names = F,sep = "\t",quote = F)

library(ggpubr)
library(Ipaper)
library(ggthemes)

## 08-1 stage-----
my_comparisons <- list(c("Stage1","Stage2"),c("Stage2","Stage3"),c("Stage1","Stage3"))
stage_data <- data.frame(riskScore = train_phenotype3$riskScore,
                         stage = factor(train_phenotype3$stage,
                                        levels = c("Stage1","Stage2", "Stage3")))
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
              y_position = c(6,7,8))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 20))+
  guides(fill='none')
stage
## 08-2 gender-----
my_comparisons <- list(c("male", "female"))
gender<-ggplot(train_phenotype3,aes(x = gender, y = riskScore, fill = gender)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "riskScore",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("Gender") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              test = t.test,
              map_signif_level = T,
              y_position = c(8))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 20))+
  guides(fill='none')
gender
## 08-2 race-----
my_comparisons <- list(c("white","asian"),c("asian","black or african american"),c("white","black or african american"))
stage_race <- data.frame(riskScore = train_phenotype3$riskScore,
                         race = factor(train_phenotype3$race,
                                       levels = c("white","asian", "black or african american")))
stage_race <- na.omit(stage_race)

race<-ggplot(stage_race,aes(x = race, y = riskScore, fill = race)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "riskScore",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("Race") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              test = t.test,
              map_signif_level = T,
              y_position = c(6,7,8))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 18),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13))+
  guides(fill='none')
race
library(patchwork)
all_clinical_index <- stage + gender + 
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
