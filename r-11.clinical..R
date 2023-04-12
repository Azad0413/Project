rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-300-8/")
if (! dir.exists("./11_clinical")){
  dir.create("./11_clinical")
}
setwd("./11_clinical")

dat.tcga<-read.delim2("../00_rawdata/dat.fpkm.xls", row.names = 1)%>% lc.tableToNum
colnames(dat.tcga)<-gsub('.','-',colnames(dat.tcga),fixed = T)
group <- read.delim2('../08_risk/risk.xls')
clinical<-read_tsv(file ='/data/nas1/luchunlin/TCGA_phenotype/TARGET-OS.clinical.tsv.gz')
#clinical<-GDCquery_clinic(project = "TARGET-OS",type = "clinical")
phenotype<-data.frame(sample=clinical$sample_id,
                      gender=clinical$Gender,
                      age=clinical$`Age at Diagnosis in Days`/365,
                      race=clinical$Race,
                      metastasis=clinical$`Disease at diagnosis`)
phenotype <- phenotype[phenotype$sample%in%colnames(dat.tcga),]
phenotype$age <- round(phenotype$age,digits = 0)
table(phenotype$race)
phenotype$race <- gsub('Unknown',NA,phenotype$race)
table(phenotype$metastasis)
phenotype$metastasis <- gsub('Metastatic (confirmed)','Metastatic',phenotype$metastasis,fixed = T)
phenotype$metastasis <- gsub('Non-metastatic (confirmed)','Non-metastatic',phenotype$metastasis,fixed = T)
phenotype$metastasis <- gsub('Non-metastatic (Confirmed)','Non-metastatic',phenotype$metastasis,fixed = T)
table(phenotype$age)
phenotype$age <- ifelse(phenotype$age>14,'>14','<=14')
group <- group[,c(1,2,6)]
colnames(group) <- c('sample','OS','Group')
group$Group <- ifelse(group$Group>median(group$Group),'High risk','Low risk')
phenotype2 <- merge(group,phenotype,by='sample')
write.table(phenotype2,file = "clinical_risk.xls",row.names = F,sep = "\t",quote = F)

table(phenotype2$Group)
Low_risk<-subset(phenotype2,Group=='Low risk')
High_risk<-subset(phenotype2,Group=='High risk')

## 卡方检验
mydata<-phenotype2
##age
age.ka<-xtabs(~mydata$age+mydata$Group,data = mydata)
age.ka

chisq.test(age.ka)
# p-value =  0.5825
##gender
gender.ka<-xtabs(~mydata$gender+mydata$Group,data = mydata)
gender.ka
chisq.test(gender.ka)

# p-value = 0.904
##race
race.ka<-xtabs(~mydata$race+mydata$Group,data = mydata)
race.ka
chisq.test(race.ka)
# p-value =  0.03233
##metastasis
metastasis.ka<-xtabs(~mydata$metastasis+mydata$Group,data = mydata)
metastasis.ka

chisq.test(metastasis.ka)
# p-value =  0.2855

# p-value = 0.05263
## 相关性（堆叠图+小提琴图）
## 堆叠图
## AGE-----
## Low_risk
age.data1 <- data.frame(group=Low_risk$Group,age=Low_risk$age)
age.data1<-as.data.frame(sort(table(age.data1$age)))
colnames(age.data1)<-c('Age','Low_risk')
## High_risk
age.data2<-data.frame(group=High_risk$Group,age=High_risk$age)
age.data2<-as.data.frame(sort(table(age.data2$age)))
colnames(age.data2)<-c('Age','High_risk')
age.data<-merge(age.data1,age.data2,by='Age')
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
## gender-----
## Low_risk
gender.data1<-data.frame(group=Low_risk$Group,gender=Low_risk$gender)
gender.data1<-as.data.frame(sort(table(gender.data1$gender)))
colnames(gender.data1)<-c('gender','Low_risk')
## High_risk
gender.data2<-data.frame(group=High_risk$Group,gender=High_risk$gender)
gender.data2<-as.data.frame(sort(table(gender.data2$gender)))
colnames(gender.data2)<-c('gender','High_risk')
gender.data<-merge(gender.data1,gender.data2,by='gender')
library(reshape2)
library(plyr)
library(RColorBrewer)
mydata <- melt(gender.data,id.vars="gender",variable.name="Group",value.name="Quantity")
gender<-ggplot(mydata,aes(x=Group,y=Quantity,fill=gender))+
  geom_bar(position = "fill",stat="identity",alpha=0.7)+
  theme_bw()+
  scale_fill_manual(values=brewer.pal(4, "Set2"))+ 
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


## race-----
## Low_risk
race.data1<-data.frame(group=Low_risk$Group,race=Low_risk$race)
race.data1<-as.data.frame(sort(table(race.data1$race)))
colnames(race.data1)<-c('race','Low_risk')
## High_risk
race.data2<-data.frame(group=High_risk$Group,race=High_risk$race)
race.data2<-as.data.frame(sort(table(race.data2$race)))
colnames(race.data2)<-c('race','High_risk')
add <- data.frame(race='Asian','High_risk'=0)
race.data2 <- rbind(race.data2,add)
race.data<-merge(race.data1,race.data2,by='race')
library(reshape2)
library(plyr)
library(RColorBrewer)
mydata <- melt(race.data,id.vars="race",variable.name="Group",value.name="Quantity")
race<-ggplot(mydata,aes(x=Group,y=Quantity,fill=race))+
  geom_bar(position = "fill",stat="identity",alpha=0.7)+
  theme_bw()+
  scale_fill_manual(values=brewer.pal(4, "Set2"))+ 
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
race

## metastasis-----
## Low_risk
metastasis.data1<-data.frame(group=Low_risk$Group,metastasis=Low_risk$metastasis)
metastasis.data1<-as.data.frame(sort(table(metastasis.data1$metastasis)))
colnames(metastasis.data1)<-c('metastasis','Low_risk')
## High_risk
metastasis.data2<-data.frame(group=High_risk$Group,metastasis=High_risk$metastasis)
metastasis.data2<-as.data.frame(sort(table(metastasis.data2$metastasis)))
colnames(metastasis.data2)<-c('metastasis','High_risk')
metastasis.data<-merge(metastasis.data1,metastasis.data2,by='metastasis')
library(reshape2)
library(plyr)
library(RColorBrewer)
mydata <- melt(metastasis.data,id.vars="metastasis",variable.name="Group",value.name="Quantity")
metastasis<-ggplot(mydata,aes(x=Group,y=Quantity,fill=metastasis))+
  geom_bar(position = "fill",stat="identity",alpha=0.7)+
  theme_bw()+
  scale_fill_manual(values=brewer.pal(4, "Set2"))+ 
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
metastasis

##小提琴图-----
risk <- read.delim2('../08_risk/risk.xls')%>%select(c('id','riskScore'))
colnames(risk) <- c('sample','riskScore')
train_phenotype3 <- merge(risk,phenotype2,by='sample')
train_phenotype3$riskScore <- as.numeric(train_phenotype3$riskScore)
## AGE
table(train_phenotype3$age)
my_comparisons <- list(c(">14","<=14"))
age_data <- data.frame(riskScore = train_phenotype3$riskScore,
                       age = factor(train_phenotype3$age,
                                    levels = c(">14","<=14")))
age_data <- na.omit(age_data)
age2<-ggplot(age_data,aes(x = age, y = riskScore, fill = age)) +
  geom_violin(trim=F) + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.3,
               position=position_dodge(0.9),
               outlier.shape = NA,fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_y_continuous(name = "riskScore",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
 # ggtitle("Age") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              test = wilcox.test,
              map_signif_level = T,
              y_position = c(3.5))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 20))+
  guides(fill='none')
age2

## gender
table(train_phenotype3$gender)
my_comparisons <- list(c("Female","Male"))
gender_data <- data.frame(riskScore = train_phenotype3$riskScore,
                       gender = factor(train_phenotype3$gender,
                                    levels = c("Female","Male")))
gender_data <- na.omit(gender_data)
gender2<-ggplot(gender_data,aes(x = gender, y = riskScore, fill = gender)) +
  geom_violin(trim=F) + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.3,
               position=position_dodge(0.9),
               outlier.shape = NA,fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_y_continuous(name = "riskScore",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
#  ggtitle("Gender") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              test = wilcox.test,
              map_signif_level = T,
              y_position = c(3.5))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 20))+
  guides(fill='none')
gender2

## race
table(train_phenotype3$race)
my_comparisons <- list(c("Asian","Black or African American"),c('Black or African American','White'),c('Asian','White'))
race_data <- data.frame(riskScore = train_phenotype3$riskScore,
                          race = factor(train_phenotype3$race,
                                          levels = c("Asian","Black or African American",'White')))
race_data <- na.omit(race_data)
race2<-ggplot(race_data,aes(x = race, y = riskScore, fill = race)) +
  geom_violin(trim=F) + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.3,
               position=position_dodge(0.9),
               outlier.shape = NA,fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_y_continuous(name = "riskScore",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
#  ggtitle("Race") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              test = wilcox.test,
              map_signif_level = T,
              y_position = c(2,3,4))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13))+
  guides(fill='none')
race2

## metastasis
table(train_phenotype3$metastasis)
my_comparisons <- list(c("Metastatic","Non-metastatic"))
metastasis_data <- data.frame(riskScore = train_phenotype3$riskScore,
                          metastasis = factor(train_phenotype3$metastasis,
                                          levels = c("Metastatic","Non-metastatic")))
metastasis_data <- na.omit(metastasis_data)
metastasis2<-ggplot(metastasis_data,aes(x = metastasis, y = riskScore, fill = metastasis)) +
  geom_violin(trim=F) + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.3,
               position=position_dodge(0.9),
               outlier.shape = NA,fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_y_continuous(name = "riskScore",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
#  ggtitle("Metastasis") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              test = wilcox.test,
              map_signif_level = T,
              y_position = c(3.5))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 20))+
  guides(fill='none')
metastasis2

library(patchwork)
age_plot <- age+age2+
  plot_layout(ncol = 2) &
  theme_bw() &
  theme(legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold",size = 12))
age_plot

ggsave(filename = '01.Age_plot.pdf',age_plot,w=7,h=4)
ggsave(filename = '01.Age_plot.png',age_plot,w=7,h=4)


library(patchwork)
gender_plot <- gender+gender2+
  plot_layout(ncol = 2) &
  theme_bw() &
  theme(legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold",size = 12))
gender_plot

ggsave(filename = '02.Gender_plot.pdf',gender_plot,w=7,h=4)
ggsave(filename = '02.Gender_plot.png',gender_plot,w=7,h=4)


library(patchwork)
race_plot <- race+race2+
  plot_layout(ncol = 2,widths = c(1,2)) &
  theme_bw() &
  theme(legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold",size = 10))
race_plot

ggsave(filename = '03.Race_plot.pdf',race_plot,w=7,h=4)
ggsave(filename = '03.Race_plot.png',race_plot,w=7,h=4)


library(patchwork)
metastasis_plot <- metastasis+metastasis2+
  plot_layout(ncol = 2) &
  theme_bw() &
  theme(legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold",size = 12))
metastasis_plot

ggsave(filename = '04.Metastasis_plot.pdf',metastasis_plot,w=7,h=4)
ggsave(filename = '04.Metastasis_plot.png',metastasis_plot,w=7,h=4)

