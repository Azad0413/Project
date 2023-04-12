rm(list = ls())
setwd("/data/nas1/luchunlin/project/JNZK-255-2/")
if (! dir.exists("./08_clinical")){
  dir.create("./08_clinical")
}
setwd("./08_clinical")

library(readr)
library(readxl)
dat.tcga<-read.csv('../00_rawdata/dat.fpkm.xls',sep = '\t',row.names = 1)
colnames(dat.tcga)<-gsub('.','-',colnames(dat.tcga),fixed = T)
group <- read.delim2('../06_risk/risk.xls')

survival<-read.delim2("/data/nas1/luchunlin/TCGA_survival/TCGA-BLCA.survival.tsv")
survival <- survival[,-3]
survival$sample <- gsub('.','-',survival$sample,fixed = T)
phenotype<-read.delim2('/data/nas1/luchunlin/TCGA_phenotype/TCGA-BLCA.GDC_phenotype.tsv.gz')
train_phenotype<-data.frame(sample=phenotype$submitter_id.samples,
                            Age=phenotype$age_at_initial_pathologic_diagnosis,
                            Gender=phenotype$gender.demographic,
                            Grade=phenotype$neoplasm_histologic_grade,
                            Stage=phenotype$tumor_stage.diagnoses,
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
train_phenotype2[train_phenotype2==''] <- NA
train_phenotype2$T.stage<-gsub('T0','T0/1/2',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('T1','T0/1/2',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('T2','T0/1/2',train_phenotype2$T.stage)

table(train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('NX',NA,train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('N2','N2/3',train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('N3','N2/3',train_phenotype2$N.stage)
# train_phenotype2$N.stage[train_phenotype2$N.stage==''] <- NA


table(train_phenotype2$M.stage)
train_phenotype2$M.stage<-gsub('MX',NA,train_phenotype2$M.stage)

table(train_phenotype2$Stage)
train_phenotype2$Stage<-gsub('not reported',NA,train_phenotype2$Stage)
# train_phenotype2$Stage<-gsub('a','',train_phenotype2$Stage)
# train_phenotype2$Stage<-gsub('b','',train_phenotype2$Stage)
# train_phenotype2$Stage<-gsub('c','',train_phenotype2$Stage)
train_phenotype2$Stage<-gsub('stage iv','Stage4',train_phenotype2$Stage)
train_phenotype2$Stage<-gsub('stage iii','Stage3',train_phenotype2$Stage)
train_phenotype2$Stage<-gsub('stage ii','Stage1/2',train_phenotype2$Stage)
train_phenotype2$Stage<-gsub('stage i','Stage1/2',train_phenotype2$Stage)

table(train_phenotype2$Grade)
# train_phenotype2$Grade<-gsub('GX',NA,train_phenotype2$Grade)
# train_phenotype2$Grade<-gsub('G1','G1/2',train_phenotype2$Grade)
# train_phenotype2$Grade<-gsub('G2','G1/2',train_phenotype2$Grade)

table(train_phenotype2$Age)
median(train_phenotype2$Age)
train_phenotype2$Age <- ifelse(train_phenotype2$Age>70,'>70','<=70')

table(train_phenotype2$Gender)

library(tidyverse)
colnames(train_phenotype2)
colnames(train_phenotype2)[1]<-'id'
risk<-read.csv('../06_risk/risk.xls',sep = '\t')
sub_risk <- subset(risk, select = c(id, riskScore))
colnames(sub_risk) <- c('id','riskScore')
sub_risk$id <- gsub('.','-',sub_risk$id,fixed = T)
train_phenotype3 <- merge(sub_risk,
                          train_phenotype2,
                          by = "id")
train_phenotype3$Group<-ifelse(train_phenotype3$riskScore>median(train_phenotype3$riskScore),'High risk','Low risk')

write.table(train_phenotype3,file = "clinical_risk.xls",row.names = F,sep = "\t",quote = F)

library(ggpubr)
library(ggthemes)
table(train_phenotype3$Stage)

# ## 08-1 TNM stage-----
my_comparisons <- list(c("Stage1/2","Stage3"),c("Stage1/2","Stage4"),c("Stage3","Stage4"))
stage_data <- data.frame(riskScore = train_phenotype3$riskScore,
                         stage = factor(train_phenotype3$Stage,
                                        levels = c("Stage1/2","Stage3","Stage4")))
stage_data <- na.omit(stage_data)

stage <- ggplot(stage_data,aes(x = stage, y = riskScore, fill = stage)) +
  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar",
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.4,
               position=position_dodge(0.9),
               outlier.shape = NA,fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#CD3700","#4682B4","darkgreen","orange"), name = "Group")+
  labs(title="", x="", y = "",size=20) +
  geom_signif(comparisons = my_comparisons,
              test = wilcox.test,
              map_signif_level = T,
              y_position = c(3,3.5,4))+
  # stat_compare_means(data = train_phenotype3,
  #                    mapping = aes(group = stage),
  #                    label ="p",
  #                    method = 'kruskal.test',
  #                    paired = F,label.y = 4.5,label.x = 1.3) +
  ggtitle('Stage')+
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
  guides(fill='none')
stage

## 08-2 T stage-----
table(train_phenotype3$T.stage)
my_comparisons <- list(c("T0/1/2","T3"),c("T0/1/2","T4"),c("T3","T4"))
Tstage_data <- data.frame(riskScore = train_phenotype3$riskScore,
                          stage = factor(train_phenotype3$T.stage,
                                         levels = c("T0/1/2","T3", "T4")))
Tstage_data <- na.omit(Tstage_data)

Tstage <- ggplot(Tstage_data,aes(x = stage, y = riskScore, fill = stage)) +
  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.4,
               position=position_dodge(0.9),
               outlier.shape = NA,fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#CD3700","#4682B4","darkgreen"), name = "Group")+
  labs(title="", x="", y = "",size=20) +
  geom_signif(comparisons = my_comparisons,
              test = wilcox.test,
              map_signif_level = T,
              y_position = c(3,3.5,4))+
  # stat_compare_means(data = train_phenotype3,
  #                    mapping = aes(group = stage),
  #                    label ="p",
  #                    method = 'kruskal.test',
  #                    paired = F,label.y = 4.5,label.x = 1.3) +
  ggtitle('T.stage')+
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
  guides(fill='none')
Tstage

## 08-2 M stage-----
table(train_phenotype3$M.stage)
my_comparisons <- list(c("M0","M1"))
Mstage_data <- data.frame(riskScore = train_phenotype3$riskScore,
                          stage = factor(train_phenotype3$M.stage,
                                         levels = c("M0","M1")))
Mstage_data <- na.omit(Mstage_data)

Mstage <- ggplot(Mstage_data,aes(x = stage, y = riskScore, fill = stage)) +
  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.4,
               position=position_dodge(0.9),
               outlier.shape = NA,fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#CD3700","#4682B4","darkgreen"), name = "Group")+
  labs(title="", x="", y = "",size=20) +
  geom_signif(comparisons = my_comparisons,
              test = wilcox.test,
              map_signif_level = T,
              y_position = c(3))+
  # stat_compare_means(data = train_phenotype3,
  #                    mapping = aes(group = stage),
  #                    label ="p",
  #                    method = 'kruskal.test',
  #                    paired = F,label.y = 4.5,label.x = 1.3) +
  ggtitle('M.stage')+
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
  guides(fill='none')
Mstage

## 08-2 N stage-----
table(train_phenotype3$N.stage)
my_comparisons <- list(c("N0","N1"),c('N1','N2/3'),c('N0','N2/3'))
Nstage_data <- data.frame(riskScore = train_phenotype3$riskScore,
                          stage = factor(train_phenotype3$N.stage,
                                         levels = c("N0","N1","N2/3")))
Nstage_data <- na.omit(Nstage_data)

Nstage <- ggplot(Nstage_data,aes(x = stage, y = riskScore, fill = stage)) +
  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.4,
               position=position_dodge(0.9),
               outlier.shape = NA,fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#CD3700","#4682B4","darkgreen",'orange'), name = "Group")+
  labs(title="", x="", y = "",size=20) +
  geom_signif(comparisons = my_comparisons,
              test = wilcox.test,
              map_signif_level = T,
              y_position = c(3,3.5,4))+
  # stat_compare_means(data = train_phenotype3,
  #                    mapping = aes(group = stage),
  #                    label ="p",
  #                    method = 'kruskal.test',
  #                    paired = F,label.y = 4.5,label.x = 1.3) +
  ggtitle('N.stage')+
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
  guides(fill='none')
Nstage

## 08-2 N stage-----
table(train_phenotype3$Grade)
my_comparisons <- list(c("High Grade","Low Grade"))
grade_data <- data.frame(riskScore = train_phenotype3$riskScore,
                         grade = factor(train_phenotype3$Grade,
                                        levels = c("High Grade","Low Grade")))
grade_data <- na.omit(grade_data)

grade <- ggplot(grade_data,aes(x = grade, y = riskScore, fill = grade)) +
  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.4,
               position=position_dodge(0.9),
               outlier.shape = NA,fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#CD3700","#4682B4","darkgreen"), name = "Group")+
  labs(title="", x="", y = "",size=20) +
  geom_signif(comparisons = my_comparisons,
              test = wilcox.test,
              map_signif_level = T,
              y_position = c(3.5))+
  # stat_compare_means(data = train_phenotype3,
  #                    mapping = aes(group = stage),
  #                    label ="p",
  #                    method = 'kruskal.test',
  #                    paired = F,label.y = 4.5,label.x = 1.3) +
  ggtitle('Grade')+
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
  guides(fill='none')
grade

## 08-2 N stage-----
table(train_phenotype3$Age)
my_comparisons <- list(c("<=70",">70"))
Age_data <- data.frame(riskScore = train_phenotype3$riskScore,
                       age = factor(train_phenotype3$Age,
                                    levels = c("<=70",">70")))
Age_data <- na.omit(Age_data)

Age <- ggplot(Age_data,aes(x = age, y = riskScore, fill = age)) +
  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.4,
               position=position_dodge(0.9),
               outlier.shape = NA,fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#CD3700","#4682B4","darkgreen"), name = "Group")+
  labs(title="", x="", y = "",size=20) +
  geom_signif(comparisons = my_comparisons,
              test = wilcox.test,
              map_signif_level = T,
              y_position = c(3.5))+
  # stat_compare_means(data = train_phenotype3,
  #                    mapping = aes(group = stage),
  #                    label ="p",
  #                    method = 'kruskal.test',
  #                    paired = F,label.y = 4.5,label.x = 1.3) +
  ggtitle('Age')+
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
  guides(fill='none')
Age

## 08-2 Gender-----
table(train_phenotype3$Gender)
my_comparisons <- list(c("female","male"))
Gender_data <- data.frame(riskScore = train_phenotype3$riskScore,
                          Gender = factor(train_phenotype3$Gender,
                                          levels = c("female","male")))
Gender_data <- na.omit(Gender_data)

Gender <- ggplot(Gender_data,aes(x = Gender, y = riskScore, fill = Gender)) +
  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.4,
               position=position_dodge(0.9),
               outlier.shape = NA,fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#CD3700","#4682B4","darkgreen"), name = "Group")+
  labs(title="", x="", y = "",size=20) +
  geom_signif(comparisons = my_comparisons,
              test = wilcox.test,
              map_signif_level = T,
              y_position = c(3.5))+
  # stat_compare_means(data = train_phenotype3,
  #                    mapping = aes(group = stGender),
  #                    label ="p",
  #                    method = 'kruskal.test',
  #                    paired = F,label.y = 4.5,label.x = 1.3) +
  ggtitle('Gender')+
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
  guides(fill='none')
Gender


library(patchwork)
all_plot <- Age+Gender+grade+stage+Tstage+Nstage+Mstage+
  plot_layout(ncol = 4) &
  theme_bw() &
  theme(legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text.x = element_text(angle = 45,hjust = 1,face = "bold",size = 12),
        axis.text.y = element_text(face = 'bold',size = 12))
all_plot

ggsave(filename = '01.clinical.pdf',all_plot,w=10,h=8)
ggsave(filename = '01.clinical.png',all_plot,w=10,h=8)

