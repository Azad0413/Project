## 03 临床指标-----
rm(list = ls())
setwd("/data/nas1/luchunlin/project/LZZK-504")
if (! dir.exists("./02_clinical")){
  dir.create("./02_clinical")
}
setwd("./02_clinical")
clinical<-read_tsv(file ='TCGA-OV.GDC_phenotype.tsv')
phenotype<-data.frame(sample=clinical$submitter_id.samples,
                      stage=clinical$clinical_stage,
                      age=clinical$age_at_index.demographic)

phenotype$stage<-gsub('A','',phenotype$stage)
phenotype$stage<-gsub('B','',phenotype$stage)
phenotype$stage<-gsub('C','',phenotype$stage)
phenotype$age<-cut(phenotype$age,breaks = c(30,40,50,60,70,80,90),labels = c('30-40','40-50','50-60','60-70','70-80','80-90'))

dat<-read.delim2('/data/nas1/luchunlin/project/LLZK-505/00_rawdata/dat.fpkm.xls',row.names = 1)%>%lc.tableToNum()
##列名改回去
colnames<-data.frame(sample=colnames(dat))
colnames$sample<-gsub('.','-',colnames$sample,fixed = T)
colnames(dat)<-colnames$sample
dat<-t(dat)%>%as.data.frame()
## CIRBP-------
group<-data.frame(sample=rownames(dat),group=ifelse(dat$CIRBP>median(dat$CIRBP),'High CIRBP','Low CIRBP'),expression=log2(dat$CIRBP+1))
write.table(group,file = 'group.xls',sep = '\t',quote = F,row.names = F)
phenotype$stage<-gsub('IV','4',phenotype$stage)
phenotype$stage<-gsub('III','3',phenotype$stage)
phenotype$stage<-gsub('II','1-2',phenotype$stage)
phenotype$stage<-gsub('I','1-2',phenotype$stage)
table(phenotype$age)
### 30-50  50-70  70-90  
phenotype$age<-gsub('30-40','30-50',phenotype$age)
phenotype$age<-gsub('40-50','30-50',phenotype$age)
phenotype$age<-gsub('50-60','50-70',phenotype$age)
phenotype$age<-gsub('60-70','50-70',phenotype$age)
phenotype$age<-gsub('70-80','70-90',phenotype$age)
phenotype$age<-gsub('80-90','70-90',phenotype$age)
colnames(phenotype)
colnames(phenotype)<-c('sample','TNM stage','Age')
phenotype2 <- merge(phenotype,group1,by = "sample")
library(ggpubr)
library(Ipaper)
library(ggthemes)
table(phenotype$`TNM stage`)
my_comparisons <- list(c("Stage 1-2","Stage 3"),c("Stage 1-2","Stage 4"),c("Stage 3","Stage 4"))
stage_data <- data.frame(expression = phenotype2$expression,
                         stage = factor(phenotype2$`TNM stage`,
                                        levels = c("Stage 1-2","Stage 3","Stage 4")))
stage_data <- na.omit(stage_data)
stage<-ggplot(stage_data,aes(x = stage, y = expression, fill = stage)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "log2(expression+1)",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("stage(CIRBP)") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              map_signif_level = T,
              y_position = c(15,17,16))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 17),
        axis.title.y = element_text(size = 17))+
  guides(fill='none')
stage
table(phenotype$Age)
my_comparisons <- list(c("30-50","50-70"),c("50-70","70-90"),c("30-50","70-90"))
age_data <- data.frame(expression = phenotype2$expression,
                       age = factor(phenotype2$Age,
                                    levels = c("30-50","50-70","70-90")))
age_data <- na.omit(age_data)
age<-ggplot(age_data,aes(x = age, y = expression, fill = age)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "log2(expression+1)",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("age(CIRBP)") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              map_signif_level = T,
              y_position = c(15,16,17))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 17),
        axis.title.y = element_text(size = 17))+
  guides(fill='none')
age
library(patchwork)
all_clinical_index <- stage + age +
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
