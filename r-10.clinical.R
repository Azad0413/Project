rm(list = ls())
setwd("/data/nas1/luchunlin/project/JNZK-214-8/")
if (! dir.exists("./10_clinical")){
  dir.create("./10_clinical")
}
setwd("./10_clinical")
library(readxl)                                               
library(readr)
survival<-read.delim2('../07_univariate_cox/survival.xls')
phenotype<-read.delim2('/data/nas1/luchunlin/TCGA_phenotype/TCGA-COAD.GDC_phenotype.tsv.gz')

train_phenotype<-data.frame(sample=phenotype$submitter_id.samples,
                            age=phenotype$age_at_initial_pathologic_diagnosis,
                            gender=phenotype$gender.demographic,
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
train_phenotype2$T.stage<-gsub('Tis',NA,train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('T1','T1/2',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('T2','T1/2',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('T3','T3/4',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('T4','T3/4',train_phenotype2$T.stage)
table(train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('a','',train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('b','',train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('c','',train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('N1','N1/2',train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('N2','N1/2',train_phenotype2$N.stage)
table(train_phenotype2$M.stage)
train_phenotype2$M.stage<-gsub('a','',train_phenotype2$M.stage)
train_phenotype2$M.stage<-gsub('b','',train_phenotype2$M.stage)
train_phenotype2$M.stage<-gsub('MX',NA,train_phenotype2$M.stage)
train_phenotype2[train_phenotype2=='']<-NA 
table(train_phenotype2$age)
train_phenotype2$age<-ifelse(train_phenotype2$age>60,'>60','<=60')
train_phenotype2$status <-ifelse(train_phenotype2$OS=='1','Dead','Alive') 
table(train_phenotype2$gender)

library(tidyverse)
library(lance)
colnames(train_phenotype2)
colnames(train_phenotype2)<-c('id','Age','Gender','T.stage','N.stage','M.stage','OS','OS.time','Survival')
risk<-read.delim2('../09_risk/risk.xls')%>%lc.tableToNum()
sub_risk <- subset(risk, select = c(id, riskScore))
train_phenotype3 <- merge(sub_risk,
                          train_phenotype2,
                          by = "id")
train_phenotype3$Group<-ifelse(train_phenotype3$riskScore>median(train_phenotype3$riskScore),'High risk','Low risk')

write.table(train_phenotype3,file = "clinical_risk.xls",row.names = F,sep = "\t",quote = F)

table(train_phenotype3$Group)
Low_risk<-subset(train_phenotype3,Group=='Low risk')
High_risk<-subset(train_phenotype3,Group=='High risk')

## AGE-----
## Low_risk
age.data1 <- data.frame(group=Low_risk$Group,age=Low_risk$Age)
age.data1<-as.data.frame(sort(table(age.data1$age)))
colnames(age.data1)<-c('Age','Low_risk')
## High_risk
age.data2<-data.frame(group=High_risk$Group,age=High_risk$Age)
age.data2<-as.data.frame(sort(table(age.data2$age)))
colnames(age.data2)<-c('Age','High_risk')
age.data<-merge(age.data1,age.data2,by='Age')

# i<-1
# while(i<4){
#   age.data$All[i]<-sum(age.data$Low_risk[i],age.data$High_risk[i])
#   i<-i+1
# }

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
ggsave(filename = '01.age.pdf',age,w=7,h=6)
ggsave(filename = '01.age.png',age,w=7,h=6)

## gender-----
## Low_risk
gender.data1<-data.frame(group=Low_risk$Group,gender=Low_risk$Gender)
gender.data1<-as.data.frame(sort(table(gender.data1$gender)))
colnames(gender.data1)<-c('gender','Low_risk')
## High_risk
gender.data2<-data.frame(group=High_risk$Group,gender=High_risk$Gender)
gender.data2<-as.data.frame(sort(table(gender.data2$gender)))
colnames(gender.data2)<-c('gender','High_risk')
gender.data<-merge(gender.data1,gender.data2,by='gender')

# i<-1
# while(i<4){
#   gender.data$All[i]<-sum(gender.data$Low_risk[i],gender.data$High_risk[i])
#   i<-i+1
# }

## 堆叠图
library(reshape2)
library(plyr)
library(RColorBrewer)
mydata <- melt(gender.data,id.vars="gender",variable.name="Group",value.name="Quantity")
gender<-ggplot(mydata,aes(x=Group,y=Quantity,fill=gender))+
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
gender
ggsave(filename = '02.gender.pdf',gender,w=7,h=6)
ggsave(filename = '02.gender.png',gender,w=7,h=6)

## T.stage-----
## Low_risk
T.stage.data1<-data.frame(group=Low_risk$Group,T.stage=Low_risk$T.stage)
T.stage.data1<-as.data.frame(sort(table(T.stage.data1$T.stage)))
colnames(T.stage.data1)<-c('T.stage','Low_risk')
## High_risk
T.stage.data2<-data.frame(group=High_risk$Group,T.stage=High_risk$T.stage)
T.stage.data2<-as.data.frame(sort(table(T.stage.data2$T.stage)))
colnames(T.stage.data2)<-c('T.stage','High_risk')
T.stage.data<-merge(T.stage.data1,T.stage.data2,by='T.stage')

# i<-1
# while(i<4){
#   T.stage.data$All[i]<-sum(T.stage.data$Low_risk[i],T.stage.data$High_risk[i])
#   i<-i+1
# }

## 堆叠图
library(reshape2)
library(plyr)
library(RColorBrewer)
mydata <- melt(T.stage.data,id.vars="T.stage",variable.name="Group",value.name="Quantity")
Tstage<-ggplot(mydata,aes(x=Group,y=Quantity,fill=T.stage))+
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
Tstage
ggsave(filename = '03.T.stage.pdf',Tstage,w=7,h=6)
ggsave(filename = '03.T.stage.png',Tstage,w=7,h=6)

## N.stage-----
## Low_risk
N.stage.data1<-data.frame(group=Low_risk$Group,N.stage=Low_risk$N.stage)
N.stage.data1<-as.data.frame(sort(table(N.stage.data1$N.stage)))
colnames(N.stage.data1)<-c('N.stage','Low_risk')
## High_risk
N.stage.data2<-data.frame(group=High_risk$Group,N.stage=High_risk$N.stage)
N.stage.data2<-as.data.frame(sort(table(N.stage.data2$N.stage)))
colnames(N.stage.data2)<-c('N.stage','High_risk')
N.stage.data<-merge(N.stage.data1,N.stage.data2,by='N.stage')

# i<-1
# while(i<4){
#   N.stage.data$All[i]<-sum(N.stage.data$Low_risk[i],N.stage.data$High_risk[i])
#   i<-i+1
# }

## 堆叠图
library(reshape2)
library(plyr)
library(RColorBrewer)
mydata <- melt(N.stage.data,id.vars="N.stage",variable.name="Group",value.name="Quantity")
Nstage<-ggplot(mydata,aes(x=Group,y=Quantity,fill=N.stage))+
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
Nstage
ggsave(filename = '04.N.stage.pdf',Nstage,w=7,h=6)
ggsave(filename = '04.N.stage.png',Nstage,w=7,h=6)

## M.stage-----
## Low_risk
M.stage.data1<-data.frame(group=Low_risk$Group,M.stage=Low_risk$M.stage)
M.stage.data1<-as.data.frame(sort(table(M.stage.data1$M.stage)))
colnames(M.stage.data1)<-c('M.stage','Low_risk')
## High_risk
M.stage.data2<-data.frame(group=High_risk$Group,M.stage=High_risk$M.stage)
M.stage.data2<-as.data.frame(sort(table(M.stage.data2$M.stage)))
colnames(M.stage.data2)<-c('M.stage','High_risk')
M.stage.data<-merge(M.stage.data1,M.stage.data2,by='M.stage')

# i<-1
# while(i<4){
#   M.stage.data$All[i]<-sum(M.stage.data$Low_risk[i],M.stage.data$High_risk[i])
#   i<-i+1
# }

## 堆叠图
library(reshape2)
library(plyr)
library(RColorBrewer)
mydata <- melt(M.stage.data,id.vars="M.stage",variable.name="Group",value.name="Quantity")
Mstage<-ggplot(mydata,aes(x=Group,y=Quantity,fill=M.stage))+
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
Mstage
ggsave(filename = '05.M.stage.pdf',Mstage,w=7,h=6)
ggsave(filename = '05.M.stage.png',Mstage,w=7,h=6)

## survival-----
## Low_risk
survival.data1<-data.frame(group=Low_risk$Group,survival=Low_risk$Survival)
survival.data1<-as.data.frame(sort(table(survival.data1$survival)))
colnames(survival.data1)<-c('survival','Low_risk')
## High_risk
survival.data2<-data.frame(group=High_risk$Group,survival=High_risk$Survival)
survival.data2<-as.data.frame(sort(table(survival.data2$survival)))
colnames(survival.data2)<-c('survival','High_risk')
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
ggsave(filename = '06.survival.pdf',survival,w=7,h=6)
ggsave(filename = '06.survival.png',survival,w=7,h=6)

library(patchwork)
all_clinical_index <-age+gender + Tstage + Nstage+ Mstage+ survival+
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
ggsave('07.clinical_all.pdf',all_clinical_index,w=8,h=6)
ggsave('07.clinical_all.png',all_clinical_index,w=8,h=6)



### 高低风险组卡方检验
mydata<-train_phenotype3
##age
age.ka<-xtabs(~mydata$Age+mydata$Group,data = mydata)
chisq.test(age.ka)
##X-squared = 0.66755, df = 1, p-value = 0.4139
##gender
gender.ka<-xtabs(~mydata$Gender+mydata$Group,data = mydata)
gender.ka
chisq.test(gender.ka)
##X-squared = 0.89688, df = 1, p-value = 0.3436
##T.Stage
T.stage.ka<-xtabs(~mydata$T.stage+mydata$Group,data = mydata)
T.stage.ka
chisq.test(T.stage.ka)
##X-squared = 1.4594, df = 1, p-value = 0.227
##N.stage
N.stage.ka<-xtabs(~mydata$N.stage+mydata$Group,data = mydata)
N.stage.ka
chisq.test(N.stage.ka)
# X-squared = 11.298, df = 1, p-value = 0.000776
##M.stage
M.stage.ka<-xtabs(~mydata$M.stage+mydata$Group,data = mydata)
M.stage.ka
chisq.test(M.stage.ka)
# X-squared = 8.706, df = 1, p-value = 0.003172
##survival
survival.ka<-xtabs(~mydata$Survival+mydata$Group,data = mydata)
survival.ka
chisq.test(survival.ka)
# X-squared = 16.241, df = 1, p-value = 5.579e-05