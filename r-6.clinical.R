rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-441-3/")
if (! dir.exists("./06_clinical")){
  dir.create("./06_clinical")
}
setwd("./06_clinical")
dat<-read.csv('../00_rawdata/dat.fpkm.xls',sep = '\t',row.names = 1)
colnames(dat) <- gsub('.','-',colnames(dat),fixed = T)
hub_exp<-dat['UGCG',]
group <- read.delim2('../05_survival/group(UGCG).xls')

survival<-read.delim2("/data/nas1/luchunlin/TCGA_survival/TCGA-LAML.survival.tsv")
survival <- survival[,-3]
survival$sample <- gsub('.','-',survival$sample,fixed = T)
phenotype<-read.delim2('/data/nas1/luchunlin/TCGA_phenotype/TCGA-LAML.GDC_phenotype.tsv.gz')
train_phenotype<-data.frame(sample=phenotype$submitter_id.samples,
                            Age=phenotype$age_at_initial_pathologic_diagnosis,
                            Gender=phenotype$gender.demographic,
                            FAB_classification=phenotype$leukemia_french_american_british_morphology_code)

train_phenotype<-merge(train_phenotype,survival,by='sample')
table(train_phenotype$FAB_classification)
train_phenotype$FAB_classification<-gsub(' Undifferentiated','',train_phenotype$FAB_classification,fixed = T)
train_phenotype$FAB_classification<-gsub('Not Classified',NA,train_phenotype$FAB_classification)

table(train_phenotype$Age)
median(train_phenotype$Age)
train_phenotype$Age <- ifelse(train_phenotype$Age>60,'>60','<=60')

table(train_phenotype$Gender)

colnames(train_phenotype)

hub.dat <- log2(t(hub_exp)+1)%>%as.data.frame()%>%rownames_to_column(var = 'sample')


train_phenotype2 <- merge(train_phenotype,
                          hub.dat,
                          by = "sample")
write.table(train_phenotype2,
            file = "clinical.csv",
            row.names = T,
            sep = "\t",
            quote = F)
library(ggpubr)
library(Ipaper)
library(ggthemes)
## 10-1 gender-----
my_comparisons <- list(c("male", "female"))
gender<-ggplot(train_phenotype2,aes(x = Gender, y = UGCG, fill = Gender)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "UGCG expression level",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("Gender") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              map_signif_level = T,
              y_position = c(7))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 20))+
  guides(fill='none')
gender

## 10-2 age-----
my_comparisons <- list(c(">60", "<=60"))
age<-ggplot(train_phenotype2,aes(x = Age, y = UGCG, fill = Age)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "UGCG expression level",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("Age") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              map_signif_level = T,
              y_position = c(7))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 20))+
  guides(fill='none')
age

## 10-2 FAB_classification-----
my_comparisons <- list(c("M0","M1"),c("M1","M2"),c("M2","M3"),c("M3","M4"),c("M4","M5"),c("M5","M6"),c("M6","M7"))
stage_data <- data.frame(expression = train_phenotype2$UGCG,
                         stage = factor(train_phenotype2$FAB_classification,
                                        levels = c("M0","M1","M2","M3","M4","M5","M6","M7")))
stage_data <- na.omit(stage_data)
stage<-ggplot(stage_data,aes(x = stage, y = expression, fill = stage)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "UGCG expression level",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("FAB classification") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              map_signif_level = T,
              y_position = c(7,8,9,7,8,9,7))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 20))+
  guides(fill='none')
stage
library(patchwork)
all_clinical_index <- age + gender + stage +
  plot_layout(ncol = 3,widths = c(1,1,2)) & 
  theme_bw() & 
  theme(legend.title = element_blank(),
        legend.position = "top",  
        legend.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold",size = 12))
all_clinical_index

ggsave(filename = '01.clinical.pdf',all_clinical_index,w=10,h=4)
ggsave(filename = '01.clinical.png',all_clinical_index,w=10,h=4)
