rm(list = ls())
setwd("/data/nas1/luchunlin/project/SJZZK-431-10/")
if (! dir.exists("./10_clinical")){
  dir.create("./10_clinical")
}
setwd("./10_clinical")
dat<-read.delim2("../00_rawdata/dat.xls", row.names = 1)%>% lc.tableToNum
group <- read.delim2('../08_risk/risk.xls')
dat <- dat[,group$id]
clinical<-read.delim2('../06_univariate_cox/E-MTAB-6389.sdrf.txt',header = T)
# 性别、血管浸润、坏死、诊断症状、饮酒、生存、感染、肝硬化
phenotype<-data.frame(sample=clinical$Source.Name,
                      gender=clinical$Characteristics.sex.,
                      vascular.invasion=clinical$Characteristics.vascular.invasion.,
                      necrosis=clinical$Characteristics.necrosis.detected.,
                      diagnosic.symptoms=clinical$Characteristics.diagnosic.symptoms.,
                      alcohol=clinical$Characteristics.alcohol.drinking.,
                      infection=clinical$Characteristics.infection.,
                      cirrhosis=clinical$Characteristics.cirrhosis.of.liver.,
                      survival=clinical$Characteristics.event.death.)
phenotype <- phenotype[phenotype$sample%in%colnames(dat),]

table(phenotype$gender)
table(phenotype$vascular.invasion)
table(phenotype$necrosis)
table(phenotype$diagnosic.symptoms)
table(phenotype$alcohol)
table(phenotype$infection)
table(phenotype$cirrhosis)
table(phenotype$survival)
phenotype$survival <- ifelse(phenotype$survival==1,'Dead','Alive')

group <- group%>%select(c('id','riskScore'))
colnames(group) <- c('sample','Group')
group$Group <- ifelse(group$Group>median(group$Group),'High risk','Low risk')
phenotype2 <- merge(group,phenotype,by='sample')

write.table(phenotype2,file = "clinical_risk.xls",row.names = F,sep = "\t",quote = F)

table(phenotype2$Group)
Low_risk<-subset(phenotype2,Group=='Low risk')
High_risk<-subset(phenotype2,Group=='High risk')

## 卡方检验-------
mydata<-phenotype2

##gender
gender.ka<-xtabs(~mydata$gender+mydata$Group,data = mydata)
gender.ka
chisq.test(gender.ka)
# X-squared = 0.65606, df = 1, p-value = 0.418

colnames(phenotype2)
##invasion
invasion.ka<-xtabs(~mydata$vascular.invasion+mydata$Group,data = mydata)
invasion.ka
chisq.test(invasion.ka)
# X-squared = 1.1003, df = 1, p-value = 0.2942
##necrosis
necrosis.ka<-xtabs(~mydata$necrosis+mydata$Group,data = mydata)
necrosis.ka

chisq.test(necrosis.ka)
# p-value = 0.4631

##diagnosic
diagnosic.ka<-xtabs(~mydata$diagnosic.symptoms+mydata$Group,data = mydata)
diagnosic.ka
chisq.test(diagnosic.ka)
# p-value = 1

##alcohol
alcohol.ka<-xtabs(~mydata$alcohol+mydata$Group,data = mydata)
alcohol.ka
chisq.test(alcohol.ka)
# p-value = 0.181

##infection
infection.ka<-xtabs(~mydata$infection+mydata$Group,data = mydata)
infection.ka
chisq.test(infection.ka)
# p-value = 1

##cirrhosis
cirrhosis.ka<-xtabs(~mydata$cirrhosis+mydata$Group,data = mydata)
cirrhosis.ka
chisq.test(cirrhosis.ka)
# p-value = 0.2917

##survival
survival.ka<-xtabs(~mydata$survival+mydata$Group,data = mydata)
survival.ka
chisq.test(survival.ka)
# p-value = 6.303e-05


##小提琴图-----
risk <- read.delim2('../08_risk/risk.xls')%>%select(c('id','riskScore'))
colnames(risk) <- c('sample','riskScore')
train_phenotype3 <- merge(risk,phenotype2,by='sample')
train_phenotype3$riskScore <- as.numeric(train_phenotype3$riskScore)

## gender
table(train_phenotype3$gender)
my_comparisons <- list(c("female","male"))
gender_data <- data.frame(riskScore = train_phenotype3$riskScore,
                          gender = factor(train_phenotype3$gender,
                                          levels = c("female","male")))
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
    ggtitle("Gender") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              test = wilcox.test,
              map_signif_level = T,
              y_position = c(-2))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 20))+
  guides(fill='none')
gender2

## vascular.invasion
table(train_phenotype3$vascular.invasion)
my_comparisons <- list(c("N","Y"))
invasion_data <- data.frame(riskScore = train_phenotype3$riskScore,
                        invasion = factor(train_phenotype3$vascular.invasion,
                                      levels = c("N","Y")))
invasion_data <- na.omit(invasion_data)
invasion2<-ggplot(invasion_data,aes(x = invasion, y = riskScore, fill = invasion)) +
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
    ggtitle("Vascular Invasion") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              test = wilcox.test,
              map_signif_level = T,
              y_position = c(-2))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13))+
  guides(fill='none')
invasion2

## necrosis
table(train_phenotype3$necrosis)
my_comparisons <- list(c("N","Y"))
necrosis_data <- data.frame(riskScore = train_phenotype3$riskScore,
                            necrosis = factor(train_phenotype3$necrosis,
                                              levels = c("N","Y")))
necrosis_data <- na.omit(necrosis_data)
necrosis2<-ggplot(necrosis_data,aes(x = necrosis, y = riskScore, fill = necrosis)) +
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
  ggtitle("Necrosis") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              test = wilcox.test,
              map_signif_level = T,
              y_position = c(-2))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13))+
  guides(fill='none')
necrosis2

## diagnosic.symptoms
table(train_phenotype3$diagnosic.symptoms)
my_comparisons <- list(c("N","Y"))
diagnosic.symptoms_data <- data.frame(riskScore = train_phenotype3$riskScore,
                            diagnosic.symptoms = factor(train_phenotype3$diagnosic.symptoms,
                                              levels = c("N","Y")))
diagnosic.symptoms_data <- na.omit(diagnosic.symptoms_data)
diagnosic.symptoms2<-ggplot(diagnosic.symptoms_data,aes(x = diagnosic.symptoms, y = riskScore, fill = diagnosic.symptoms)) +
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
  ggtitle("Diagnosic Symptoms") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              test = wilcox.test,
              map_signif_level = T,
              y_position = c(-2))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13))+
  guides(fill='none')
diagnosic.symptoms2

## alcohol
table(train_phenotype3$alcohol)
my_comparisons <- list(c("N","Y"))
alcohol_data <- data.frame(riskScore = train_phenotype3$riskScore,
                            alcohol = factor(train_phenotype3$alcohol,
                                              levels = c("N","Y")))
alcohol_data <- na.omit(alcohol_data)
alcohol2<-ggplot(alcohol_data,aes(x = alcohol, y = riskScore, fill = alcohol)) +
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
  ggtitle("Alcohol") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              test = wilcox.test,
              map_signif_level = T,
              y_position = c(-2))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13))+
  guides(fill='none')
alcohol2

## infection
table(train_phenotype3$infection)
my_comparisons <- list(c("N","Y"))
infection_data <- data.frame(riskScore = train_phenotype3$riskScore,
                            infection = factor(train_phenotype3$infection,
                                              levels = c("N","Y")))
infection_data <- na.omit(infection_data)
infection2<-ggplot(infection_data,aes(x = infection, y = riskScore, fill = infection)) +
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
  ggtitle("infection") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              test = wilcox.test,
              map_signif_level = T,
              y_position = c(-2))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13))+
  guides(fill='none')
infection2

## cirrhosis
table(train_phenotype3$cirrhosis)
my_comparisons <- list(c("N","Y"))
cirrhosis_data <- data.frame(riskScore = train_phenotype3$riskScore,
                            cirrhosis = factor(train_phenotype3$cirrhosis,
                                              levels = c("N","Y")))
cirrhosis_data <- na.omit(cirrhosis_data)
cirrhosis2<-ggplot(cirrhosis_data,aes(x = cirrhosis, y = riskScore, fill = cirrhosis)) +
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
  ggtitle("Cirrhosis") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              test = wilcox.test,
              map_signif_level = T,
              y_position = c(-2))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13))+
  guides(fill='none')
cirrhosis2

## survival
table(train_phenotype3$survival)
my_comparisons <- list(c("Dead","Alive"))
survival_data <- data.frame(riskScore = train_phenotype3$riskScore,
                            survival = factor(train_phenotype3$survival,
                                              levels = c("Dead","Alive")))
survival_data <- na.omit(survival_data)
survival2<-ggplot(survival_data,aes(x = survival, y = riskScore, fill = survival)) +
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
  ggtitle("Survival") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              test = wilcox.test,
              map_signif_level = T,
              y_position = c(-2))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13))+
  guides(fill='none')
survival2


library(patchwork)
all_plot <- gender2+invasion2+necrosis2+diagnosic.symptoms2+alcohol2+infection2+cirrhosis2+survival2+
  plot_layout(ncol = 4) &
  theme_bw() &
  theme(legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold",size = 12))
all_plot

ggsave(filename = '01.all_plot.pdf',all_plot,w=8,h=5)
ggsave(filename = '01.all_plot.png',all_plot,w=8,h=5)
