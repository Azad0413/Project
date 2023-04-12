rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/HF-0103-1/")
if (! dir.exists("./09_clinical")){
  dir.create("./09_clinical")
}
setwd("./09_clinical")

dat<-read.delim2("../00_rawdata/dat(GSE10846).xls", row.names = 1)%>% lc.tableToNum
group <- read.delim2('../07_risk/risk.xls')
phenotype<-read.delim2('../00_rawdata/phenotype.xls')
phenotype <- phenotype[phenotype$sample%in%colnames(dat),]

table(phenotype$Stage)
phenotype$Stage <- gsub('Clinical info: Stage: NA',NA,phenotype$Stage)
phenotype$Stage <- gsub('Clinical info: Stage: ','Stage',phenotype$Stage)

table(phenotype$ECOG_status)
phenotype$ECOG_status <- gsub('Clinical info: ECOG performance status: NA',NA,phenotype$ECOG_status)
phenotype$ECOG_status <- gsub('Clinical info: ECOG performance status: ','Stage',phenotype$ECOG_status)
phenotype$ECOG_status <- gsub('Stage3','Stage3/4',phenotype$ECOG_status)
phenotype$ECOG_status <- gsub('Stage4','Stage3/4',phenotype$ECOG_status)

table(phenotype$LDH_Ratio)
phenotype$LDH_Ratio <- gsub('Clinical info: LDH ratio: ','',phenotype$LDH_Ratio)
phenotype$LDH_Ratio <- gsub('NA',NA,phenotype$LDH_Ratio)
phenotype$LDH_Ratio <- as.numeric(phenotype$LDH_Ratio)
phenotype$LDH_Ratio <- ifelse(phenotype$LDH_Ratio>0.6,'>0.6','<=0.6')

table(phenotype$Age)
median(phenotype$Age)

phenotype$Age <- ifelse(phenotype$Age>60,'>60','<=60')


colnames(group)
group <- group[,c('id','OS','riskScore')]
colnames(group) <- c('sample','OS','Group')
group$Group <- as.numeric(group$Group)
group$Group <- ifelse(group$Group>median(group$Group),'High risk','Low risk')
phenotype2 <- merge(group,phenotype,by='sample')
write.table(phenotype2,file = "clinical_risk.xls",row.names = F,sep = "\t",quote = F)

table(phenotype2$Group)
Low_risk<-subset(phenotype2,Group=='Low risk')
High_risk<-subset(phenotype2,Group=='High risk')


##小提琴图-----
risk <- read.delim2('../07_risk/risk.xls')%>%select(c('id','riskScore'))
colnames(risk) <- c('sample','riskScore')
train_phenotype3 <- merge(risk,phenotype2,by='sample')
train_phenotype3$riskScore <- as.numeric(train_phenotype3$riskScore)
## AGE
table(train_phenotype3$Age)
my_comparisons <- list(c(">60","<=60"))
age_data <- data.frame(riskScore = train_phenotype3$riskScore,
                       age = factor(train_phenotype3$Age,
                                    levels = c(">60","<=60")))
age_data <- na.omit(age_data)
age<-ggplot(age_data,aes(x = age, y = riskScore, fill = age)) +
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
  ggtitle("Age") +
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
age

## gender
table(train_phenotype3$Gender)
my_comparisons <- list(c("female","male"))
gender_data <- data.frame(riskScore = train_phenotype3$riskScore,
                          gender = factor(train_phenotype3$Gender,
                                          levels = c("female","male")))
gender_data <- na.omit(gender_data)
gender<-ggplot(gender_data,aes(x = gender, y = riskScore, fill = gender)) +
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
              y_position = c(3.5))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 20))+
  guides(fill='none')
gender

## Stage
table(train_phenotype3$Stage)
my_comparisons <- list(c("Stage1","Stage2"),c('Stage2','Stage3'),c('Stage3','Stage4'),c('Stage1','Stage4'))
stage_data <- data.frame(riskScore = train_phenotype3$riskScore,
                        stage = factor(train_phenotype3$Stage,
                                      levels = c("Stage1","Stage2",'Stage3','Stage4')))
stage_data <- na.omit(stage_data)
stage<-ggplot(stage_data,aes(x = stage, y = riskScore, fill = stage)) +
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
   ggtitle("Stage") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              test = wilcox.test,
              map_signif_level = T,
              y_position = c(4,5,6,7))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13))+
  guides(fill='none')
stage

## 
table(train_phenotype3$ECOG_status)
my_comparisons <- list(c("Stage0","Stage1"),c('Stage1','Stage2'),c('Stage2','Stage3/4'),c('Stage0','Stage3/4'))
ECOG_data <- data.frame(riskScore = train_phenotype3$riskScore,
                         ECOG = factor(train_phenotype3$ECOG_status,
                                        levels = c("Stage0","Stage1",'Stage2','Stage3/4')))
ECOG_data <- na.omit(ECOG_data)
ECOG<-ggplot(ECOG_data,aes(x = ECOG, y = riskScore, fill = ECOG)) +
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
   ggtitle("ECOG Status") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              test = wilcox.test,
              map_signif_level = T,
              y_position = c(4,5,6,7))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13))+
  guides(fill='none')
ECOG

table(train_phenotype3$LDH_Ratio)
my_comparisons <- list(c("<=0.6",">0.6"))
LDH_data <- data.frame(riskScore = train_phenotype3$riskScore,
                        LDH = factor(train_phenotype3$LDH_Ratio,
                                      levels = c("<=0.6",">0.6")))
LDH_data <- na.omit(LDH_data)
LDH<-ggplot(LDH_data,aes(x = LDH, y = riskScore, fill = LDH)) +
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
   ggtitle("LDH Ratio") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              test = wilcox.test,
              map_signif_level = T,
              y_position = c(4,5,6,7))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13))+
  guides(fill='none')
LDH



library(patchwork)
all_plot <- age+gender+LDH+stage+ECOG +
  plot_layout(ncol = 3) &
  theme_bw() &
  theme(legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold",size = 12))
all_plot

ggsave(filename = '01.clinical.pdf',all_plot,w=10,h=6)
ggsave(filename = '01.clinical.png',all_plot,w=10,h=6)

