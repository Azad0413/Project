rm(list = ls())
setwd("/data/nas1/luchunlin/project/SJZZK-428-10/")
if (! dir.exists("./01_mRNAsi")){
  dir.create("./01_mRNAsi")
}
setwd("./01_mRNAsi")
library(lance)
library(tidyverse)
dat_tcga <- read.delim2('../00_rawdata/dat.fpkm.xls',row.names = 1)%>%lc.tableToNum()
colnames(dat_tcga) <- gsub('.','-',colnames(dat_tcga),fixed = T)
###样本fpkm表达量
exp.ss<-log2(dat_tcga+0.01)
##1 RNA干性评分-----
##导入模型
load('Stemness_index.rda')
Weight<-mRNAsi$Weight
names(Weight)<-mRNAsi$HUGO
class(Weight)
common <- intersect(names(Weight), rownames(exp.ss))
X <- exp.ss[common, ]
w <- Weight[common]
w <-as.numeric(w )
score <- apply(X, 2, function(z) {cor(z, w, method="sp", use="complete.obs")})
score <- score - min(score)
score <- score / max(score)
score <-data.frame(score )
colnames(score)<-'mRNAsi'
write.table(score,file = 'mRNAsi.txt',sep = '\t',quote = F)

## 与临床特征的相关性------
group <- data.frame(sample=rownames(score),group=score$mRNAsi)
group$group <- ifelse(group$group>median(group$group),'High','Low')

write.table(group,file = 'group(mRNAsi).xls',sep = '\t',row.names = F,quote = F)
##clinical--------

clinical<-read_tsv(file ='/data/nas1/luchunlin/project/BJTC-258/07_progmodel/TCGA-OV.GDC_phenotype.tsv')
phenotype<-data.frame(sample=clinical$submitter_id.samples,
                      stage=clinical$clinical_stage,
                      age=clinical$age_at_index.demographic)
group <- data.frame(id=rownames(score),mRNAsi=score$mRNAsi)
survival <- read.delim2('../07_univariate_cox/survival.xls')
colnames(survival) <- c('id','OS','OS.time')
colnames(phenotype)<-c('id','stage','age')
train_phenotype<-merge(phenotype,group,by='id')#%>%select(-'group')
train_phenotype <- merge(train_phenotype,survival,by='id')
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
train_phenotype2$mRNAsi<-as.numeric(train_phenotype2$mRNAsi)
train_phenotype2$Survival <- ifelse(train_phenotype2$OS=='1','Dead','Alive')
library(tidyverse)
library(lance)
colnames(train_phenotype2)
colnames(train_phenotype2)<-c('id','Stage','Age','mRNAsi','OS','OS.time','Survival')
train_phenotype3<-train_phenotype2
train_phenotype3$Group<-ifelse(train_phenotype3$mRNAsi>median(train_phenotype3$mRNAsi),'High','Low')
write.table(train_phenotype3,file = "clinical_score.xls",row.names = F,sep = "\t",quote = F)


## 堆叠图--------
table(train_phenotype3$Group)
Low<-subset(train_phenotype3,Group=='Low')
High<-subset(train_phenotype3,Group=='High')
## AGE-----
## Low_risk
age.data1 <- data.frame(group=Low$Group,age=Low$Age)
age.data1<-as.data.frame(sort(table(age.data1$age)))
colnames(age.data1)<-c('Age','Low')
## High_risk
age.data2<-data.frame(group=High$Group,age=High$Age)
age.data2<-as.data.frame(sort(table(age.data2$age)))
colnames(age.data2)<-c('Age','High')
age.data<-merge(age.data1,age.data2,by='Age')

## 堆叠图
library(reshape2)
library(plyr)
library(RColorBrewer)
mydata <- melt(age.data,id.vars="Age",variable.name="Group",value.name="Quantity")
age<-ggplot(mydata,aes(x=Group,y=Quantity,fill=Age))+
  geom_bar(position = "fill",stat="identity",alpha=0.7)+
  theme_bw()+
  scale_fill_manual(values=brewer.pal(3, "Set2"))+ 
  theme(axis.title.x =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.x =element_text(size=18, face = "bold",family='Times'),
        axis.title.y =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.y=element_text(size=18,  face = "bold",family='Times'),
        legend.title=element_text(size=20, color='black', face = "bold",family='Times'),
        legend.text=element_text(size=18, color='black', face = "bold",family='Times'),
        title=element_text(size=20, color='black', face = "bold",family='Times'),
        strip.text = element_text(size = 14,family = "Times", face = "bold"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(x="",y="Percentage",fill="")
age
# ggsave(filename = '01.age.pdf',age,w=7,h=6)
# ggsave(filename = '01.age.png',age,w=7,h=6)
## stage
## Low_risk
stage.data1<-data.frame(group=Low$Group,stage=Low$Stage)
stage.data1<-as.data.frame(sort(table(stage.data1$stage)))
colnames(stage.data1)<-c('stage','Low')
## High_risk
stage.data2<-data.frame(group=High$Group,stage=High$Stage)
stage.data2<-as.data.frame(sort(table(stage.data2$stage)))
colnames(stage.data2)<-c('stage','High')
stage.data<-merge(stage.data1,stage.data2,by='stage')

# i<-1
# while(i<4){
#   stage.data$All[i]<-sum(stage.data$Low_risk[i],stage.data$High_risk[i])
#   i<-i+1
# }

## 堆叠图
library(reshape2)
library(plyr)
library(RColorBrewer)
mydata <- melt(stage.data,id.vars="stage",variable.name="Group",value.name="Quantity")
stage<-ggplot(mydata,aes(x=Group,y=Quantity,fill=stage))+
  geom_bar(position = "fill",stat="identity",alpha=0.7)+
  theme_bw()+
  scale_fill_manual(values=brewer.pal(3, "Set2"))+ 
  theme(axis.title.x =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.x =element_text(size=18, face = "bold",family='Times'),
        axis.title.y =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.y=element_text(size=18,  face = "bold",family='Times'),
        legend.title=element_text(size=20, color='black', face = "bold",family='Times'),
        legend.text=element_text(size=18, color='black', face = "bold",family='Times'),
        title=element_text(size=20, color='black', face = "bold",family='Times'),
        strip.text = element_text(size = 14,family = "Times", face = "bold"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(x="",y="Percentage",fill="")
stage

## survival-----
## Low_risk
survival.data1<-data.frame(group=Low$Group,survival=Low$Survival)
survival.data1<-as.data.frame(sort(table(survival.data1$survival)))
colnames(survival.data1)<-c('survival','Low')
## High_risk
survival.data2<-data.frame(group=High$Group,survival=High$Survival)
survival.data2<-as.data.frame(sort(table(survival.data2$survival)))
colnames(survival.data2)<-c('survival','High')
survival.data<-merge(survival.data1,survival.data2,by='survival')
# i<-1
# while(i<4){
#   survival.data$All[i]<-sum(survival.data$Low_risk[i],survival.data$High_risk[i])
#   i<-i+1
# }

## 堆叠图
library(reshape2)
library(plyr)
library(RColorBrewer)
mydata <- melt(survival.data,id.vars="survival",variable.name="Group",value.name="Quantity")
survival<-ggplot(mydata,aes(x=Group,y=Quantity,fill=survival))+
  geom_bar(position = "fill",stat="identity",alpha=0.7)+
  theme_bw()+
  scale_fill_manual(values=brewer.pal(3, "Set2"))+ 
  theme(axis.title.x =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.x =element_text(size=18, face = "bold",family='Times'),
        axis.title.y =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.y=element_text(size=18,  face = "bold",family='Times'),
        legend.title=element_text(size=20, color='black', face = "bold",family='Times'),
        legend.text=element_text(size=18, color='black', face = "bold",family='Times'),
        title=element_text(size=20, color='black', face = "bold",family='Times'),
        strip.text = element_text(size = 14,family = "Times", face = "bold"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(x="",y="Percentge",fill="")
survival
library(patchwork)
all_clinical_index <-age+stage+ survival+
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
ggsave('01.clinical.pdf',all_clinical_index,w=8,h=4)
ggsave('01.clinical.png',all_clinical_index,w=8,h=4)
### 卡方检验------------
mydata<-train_phenotype3
##age
age.ka<-xtabs(~mydata$Age+mydata$Group,data = mydata)
chisq.test(age.ka)
# X-squared = 0.68345, df = 1, p-value = 0.4084

stage.ka<-xtabs(~mydata$Stage+mydata$Group,data = mydata)
stage.ka
chisq.test(stage.ka)
# X-squared = 3.7426, df = 2, p-value = 0.1539
survival.ka<-xtabs(~mydata$Survival+mydata$Group,data = mydata)
survival.ka
chisq.test(survival.ka)
# X-squared = 1.886, df = 1, p-value = 0.1697

## KM-------
library("survival")
library("survminer")
group <- data.frame(sample=rownames(score),group=score$mRNAsi)
group$group <- ifelse(group$group>median(group$group),'High','Low')
survival<-read.delim2('/data/nas1/luchunlin/TCGA_survival/TCGA-OV.survival.tsv')
survival<-survival[,-3]
dat <- merge(group,survival,by = "sample")
# res.cut<-surv_cutpoint(dat,time = 'OS.time',event = 'OS',variables = 'group')
# summary(res.cut)
# cutpoint<-'0.5793966'
# dat$group<-ifelse(dat$group>cutpoint,'High','Low')
library(survival)
##计算最佳截断值
dat$OS.time <- as.numeric(dat$OS.time)
dat$OS <- as.numeric(dat$OS)
# 去除掉"score"列
dat$group <- factor(dat$group,levels = c('High','Low'))
# 进行KM生存分析
kmfit <- survfit(Surv(OS.time,OS) ~ group, data = dat)

# 绘制KM曲线
KM <- ggsurvplot(kmfit,
                 pval = TRUE, 
                 conf.int = F,
                 legend.labs=c("High","Low" ),
                 legend.title="group",
                 title="KM",
                 font.main = c(15,"bold"),
                 risk.table = TRUE, 
                 risk.table.col = "strata", 
                 linetype = "strata", 
                 surv.median.line = "hv", 
                 ggtheme = theme_bw(), 
                 palette = c("#A73030FF", "#0073C2FF"))
KM
# surv_pvalue(kmfit,method = 'Gehan-Breslow')



# library(ggpubr)
# library(Ipaper)
# library(ggthemes)
# table(train_phenotype2$Stage)
# ## 08-1 TNM stage-----
# my_comparisons <- list(c("Stage1/2","Stage3"),c("Stage3","Stage4"),c("Stage1/2","Stage4"))
# stage_data <- data.frame(mRNAsi = train_phenotype2$mRNAsi,
#                          stage = factor(train_phenotype2$Stage,
#                                         levels = c("Stage1/2","Stage3","Stage4")))
# stage_data <- na.omit(stage_data)
# stage<-ggplot(stage_data,aes(x = stage, y = mRNAsi, fill = stage)) +
#   geom_boxplot(alpha=0.7) +
#   scale_y_continuous(name = "mRNAsi",expand = c(0.1,0.1))+
#   scale_x_discrete(name = "") +
#   ggtitle("Stage") +
#   theme_bw() +
#   geom_signif(comparisons = my_comparisons,
#               test = wilcox.test,
#               map_signif_level = T,
#               y_position = c(1.5,2,2.5))+
#   theme(plot.title = element_text(size = 16, face =  "bold"),
#         text = element_text(size = 20),
#         axis.title = element_text(face="bold"),
#         axis.text.x=element_text(size = 20))+
#   guides(fill='none')
# stage
# ## AGE-------
# table(train_phenotype2$Age)
# my_comparisons <- list(c(">60","<=60"))
# 
# age_data <- data.frame(mRNAsi = train_phenotype2$mRNAsi,
#                        age = factor(train_phenotype2$Age,
#                                     levels = c(">60","<=60")))
# age_data <- na.omit(age_data)
# 
# age<-ggplot(age_data,aes(x = age, y = mRNAsi, fill = age)) +
#   geom_boxplot(alpha=0.7) +
#   scale_y_continuous(name = "mRNAsi",expand = c(0.1,0.1))+
#   scale_x_discrete(name = "") +
#   ggtitle("Age") +
#   theme_bw() +
#   geom_signif(comparisons = my_comparisons,
#               test = wilcox.test,
#               map_signif_level = T,
#               y_position = c(1.5))+
#   theme(plot.title = element_text(size = 16, face =  "bold"),
#         text = element_text(size = 20),
#         axis.title = element_text(face="bold"),
#         axis.text.x=element_text(size = 20))+
#   guides(fill='none')
# age
# 
# library(patchwork)
# all_clinical_index <- stage+age+
#   plot_layout(ncol = 2) &
#   theme_bw() &
#   theme(legend.title = element_blank(),
#         legend.position = "top",
#         legend.text = element_text(face = "bold"),
#         plot.title = element_text(hjust = 0.5, face = "bold"),
#         panel.grid = element_blank(),
#         axis.title = element_text(face = "bold"),
#         axis.text = element_text(face = "bold",size = 12))
# all_clinical_index
# ggsave('01.clinical.pdf',all_clinical_index,w=6,h=4)
# ggsave('01.clinical.png',all_clinical_index,w=6,h=4)