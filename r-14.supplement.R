rm(list = ls())
setwd("/data/nas1/luchunlin/project/LZZK-512/")
if (! dir.exists("./14_supplement")){
  dir.create("./14_supplement")
}
setwd("./14_supplement")

# A ----
##riskScore在不同临床指标pathologic_M（远处转移）、pathologic_N（淋巴结转移）、pathologic_T（肿瘤原发灶情况）的差异。
library(readr)
library(readxl)
survival<-read.delim2("/data/nas1/luchunlin/project/LZZK-512/03_univariate_cox/survival.xls")
phenotype<-read_tsv(file = '../08_clinical_index/TCGA-BLCA.GDC_phenotype.tsv.gz')

train_phenotype<-data.frame(sample=phenotype$submitter_id.samples,
                            T.pathologic=phenotype$pathologic_T,
                            N.pathologic=phenotype$pathologic_N,
                            M.pathologic=phenotype$pathologic_M
)
train_phenotype<-merge(train_phenotype,survival,by='sample')
write.table(train_phenotype,file = 'phenotype.xls',sep = '\t',row.names = F,quote = F)
train_phenotype2<-train_phenotype
train_phenotype2$OS<-as.numeric(train_phenotype2$OS)
train_phenotype2$OS.time<-as.numeric(train_phenotype2$OS.time)

table(train_phenotype2$T.pathologic)
train_phenotype2$T.pathologic<-gsub('a','',train_phenotype2$T.pathologic)
train_phenotype2$T.pathologic<-gsub('b','',train_phenotype2$T.pathologic)
train_phenotype2$T.pathologic<-gsub('TX',NA,train_phenotype2$T.pathologic)
train_phenotype2$T.pathologic<-gsub('T0','T0/1',train_phenotype2$T.pathologic)
train_phenotype2$T.pathologic<-gsub('T1','T0/1',train_phenotype2$T.pathologic)

table(train_phenotype2$N.pathologic)
train_phenotype2$N.pathologic<-gsub('NX',NA,train_phenotype2$N.pathologic)
train_phenotype2$N.pathologic<-gsub('N2','N2/3',train_phenotype2$N.pathologic)
train_phenotype2$N.pathologic<-gsub('N3','N2/3',train_phenotype2$N.pathologic)

table(train_phenotype2$M.pathologic)
train_phenotype2$M.pathologic<-gsub('MX',NA,train_phenotype2$M.pathologic)

library(tidyverse)
library(lance)
colnames(train_phenotype2)
colnames(train_phenotype2)<-c('id','T.pathologic','N.pathologic','M.pathologic','OS','OS.time')
risk<-read.delim2('/data/nas1/luchunlin/project/LZZK-512/06_risk/risk.xls')%>%lc.tableToNum()
sub_risk <- subset(risk, select = c(id, riskScore))
train_phenotype3 <- merge(sub_risk,
                          train_phenotype2,
                          by = "id")
train_phenotype3$Group<-ifelse(train_phenotype3$riskScore>median(train_phenotype3$riskScore),'High risk','Low risk')

write.table(train_phenotype3,file = "clinical_risk.xls",row.names = F,sep = "\t",quote = F)

library(ggpubr)
library(Ipaper)
library(ggthemes)
table(train_phenotype3$T.pathologic)
## 08-1 T-----
my_comparisons <- list(c("T0/1","T2"),c('T2','T3'),c('T3','T4'))
stage_data <- data.frame(riskScore = train_phenotype3$riskScore,
                         stage = factor(train_phenotype3$T.pathologic,
                                        levels = c("T0/1","T2",'T3','T4')))
stage_data <- na.omit(stage_data)
T.stage<-ggplot(stage_data,aes(x = stage, y = riskScore, fill = stage)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "riskScore",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("Pathologic T") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              test = wilcox.test,
              map_signif_level = T,
              y_position = c(3,4,5,6))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 20))+
  guides(fill='none')
T.stage

## 08-2 N stage-----

table(train_phenotype3$N.pathologic)
my_comparisons <- list(c("N0","N1"),c("N1","N2/3"))
N.stage_data <- data.frame(riskScore = train_phenotype3$riskScore,
                           stage = factor(train_phenotype3$N.pathologic,
                                          levels = c("N0","N1","N2/3")))
N.stage_data <- na.omit(N.stage_data)
N.stage<-ggplot(N.stage_data,aes(x = stage, y = riskScore, fill = stage)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "riskScore",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("Pathologic N") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              test = wilcox.test,
              map_signif_level = T,
              y_position = c(2,3,4))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 20))+
  guides(fill='none')
N.stage
## 08-2 M stage-----

table(train_phenotype3$M.pathologic)
my_comparisons <- list(c("M0","M1"))
M.stage_data <- data.frame(riskScore = train_phenotype3$riskScore,
                           stage = factor(train_phenotype3$M.pathologic,
                                          levels = c("M0","M1")))
M.stage_data <- na.omit(M.stage_data)
M.stage<-ggplot(M.stage_data,aes(x = stage, y = riskScore, fill = stage)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "riskScore",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("Pathologic M") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              test = wilcox.test,
              map_signif_level = T,
              y_position = c(2,3,4))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 20))+
  guides(fill='none')
M.stage


library(patchwork)
all_clinical_index <- T.stage + N.stage + M.stage +
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
ggsave('clinical.pdf',all_clinical_index,w=8,h=4)
ggsave('clinical.png',all_clinical_index,w=8,h=4)

# B GSE169455 ---------
library(GEOquery)
library(Biobase)
gset<-getGEO("GSE169455",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-read_xlsx(file = 'GSE169455_normalized_by_gene.xlsx')%>%column_to_rownames('ID_REF')
a=gset[[1]]
dat<-expr
pd<-pData(a)
group<-data.frame(sample=pd$title,group=pd$`lundtax rna subtype:ch1`)
table(group$group)
group <- group[order(group$group),]
dat <- dat[,group$sample]
write.table(group,file = 'group1.xls',sep = '\t',row.names = F,quote = F)
hubgene<-read.delim2('../04_Lasso/lasso_genes.csv',header = F)
hubgene[12,1]<-'GALNT17'
hub_exp<-dat[hubgene$V1,]
hubgene[12,1]<-'GALNTL6'
rownames(hub_exp) <- hubgene$V1
hub_exp <- na.omit(hub_exp)
hub_exp2<-hub_exp
hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))

table(group$group)
## 样本分组
hub_exp2$Group<-c(rep('Ba/Sq',196),rep('GU',287),rep('Mes-like',28),rep('NE-like',91),rep('UroA',168),rep('UroB',105),rep('UroC',168))
hub_exp2$Group <- factor(hub_exp2$Group,levels = c('Ba/Sq','GU','Mes-like','NE-like','UroA','UroB','UroC'))


##分面图形
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
exp_plot <- ggplot(hub_exp2,aes(x = Group, y = expr, fill = Symbol)) +
  #geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#4682B4","#CD3700","#FF6A6A", "#20B2AA",'#FFCCCC','#77FFCC','#D28EFF','#FA8072'), name = "Group")+
  labs(title="", x="", y = "",size=20) +
  # stat_compare_means(data = hub_exp2,
  #                    mapping = aes(group = Group),
  #                    label ="p.signif",
  #                    method = 'wilcox.test') +
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=15),
        axis.text.x=element_text(angle=0,hjust=0.5,colour="black",face="bold",size=12), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())#+facet_wrap(~Group,scales = "free",nrow = 4) 
exp_plot
ggsave('02.subtype.pdf',exp_plot,w=10,h=5)
ggsave('02.subtype.png',exp_plot,w=10,h=5)
##c------
group2 <- read_xlsx('Pathologic.xlsx')
table(group2$Pathologic.response)
colnames(group2) <- c('sample','group')

group2 <- group2[order(group2$group),]
dat2 <- dat[,group2$sample]
hubgene<-read.delim2('../04_Lasso/lasso_genes.csv',header = F)
hubgene[12,1]<-'GALNT17'
hub_exp<-dat2[hubgene$V1,]
hubgene[12,1]<-'GALNTL6'
rownames(hub_exp) <- hubgene$V1
hub_exp <- na.omit(hub_exp)
hub_exp2<-hub_exp
hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))

table(group2$group)
## 样本分组
hub_exp2$Group<-c(rep('no pR',609),rep('pCR',336),rep('pPR',98))
hub_exp2$Group <- factor(hub_exp2$Group,levels = c('no pR','pCR','pPR'))

##分面图形
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
exp_plot <- ggplot(hub_exp2,aes(x = Group, y = expr, fill = Symbol)) +
  #geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#4682B4","#CD3700", "#20B2AA",'#FFCCCC','#77FFCC','#D28EFF','#FA8072'), name = "Group")+
  labs(title="", x="", y = "",size=20) +
  # stat_compare_means(data = hub_exp2,
  #                    mapping = aes(group = Group),
  #                    label ="p.signif",
  #                    method = 'wilcox.test') +
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=15),
        axis.text.x=element_text(angle=0,hjust=0.5,colour="black",face="bold",size=12), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())#+facet_wrap(~Group,scales = "free",nrow = 1) 
exp_plot
ggsave('03.pathologic.pdf',exp_plot,w=9,h=5)
ggsave('03.pathologic.png',exp_plot,w=9,h=5)
