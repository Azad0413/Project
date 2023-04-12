## rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/BJTC-304")
if (! dir.exists("./00_rawdata")){
  dir.create("./00_rawdata")
}
setwd("./00_rawdata")
library(readr)
library(readxl)
library(tidyverse)
exp<-read_xlsx('/data/nas1/luchunlin/project/BJTC-304/00_rawdata/exp.xlsx')
colnames(exp)
exp<-column_to_rownames(exp,var = 'systematic_name')
colnames(exp)<-gsub('(normalized)','',colnames(exp),fixed = T)
colnames(exp)<-gsub('(raw)','',colnames(exp),fixed = T)
dat<-exp[,c(16:30)]
#dat<-log2(dat)
colnames(dat)

## test1 test2 control用不同的颜色标注
sample_test1 <-c('[2012--1694]', '[2017--3473]','[2019--3320]','[2020--0392]','[2021--0295]')
sample_test2<-c('[2018--3932]','[2019--4471]','[2020--1074]','[2020--2386]','[2021--1035]')
sample_control<-c('[2019--2031]','[2020--2441]','[2021--0785]','[2021--882]','[2020--2702]')
dat_test1<-dat[,sample_test1]
dat_test2<-dat[,sample_test2]
dat_control<-dat[,sample_control]

# 02 MF差异表达miRNA鉴定----------
setwd("/data/nas1/luchunlin/project/BJTC-304")
if (! dir.exists("./01_DEGs")){
  dir.create("./01_DEGs")
}
setwd("./01_DEGs")
## p<0.05, log2fc>1,limma
intrest<-data.frame(mir=c('hsa-miR-181a-5p','hsa-miR-205-5p'))
## 02-1 test1 vs control----
dat1<-cbind(dat_test1,dat_test2)
dat1<-cbind(dat1,dat_control)
## 分组矩阵

group1<-data.frame(sample=colnames(dat1),
                   group=c(rep('test',10),rep('control',5)))
type1<-group1[,2]
design1 <- model.matrix(~ -1+factor(type1,levels=c('control','test'))) 
colnames(design1)<-c('control','test')
rownames(design1)<-group1$sample
library(limma)
# 对每一个基因进行线性模型构建
fit1=lmFit(dat1,design1)
# 构建比较矩阵
contrast.matrix=makeContrasts(ControlVSMG=test-control,levels = design1)
# 构建芯片数据线性模型，计算估计的相关系数和标准差
fit2_1=contrasts.fit(fit1,contrast.matrix)
# 基于贝叶斯计算T值，F值和log-odds
fit2_1<-eBayes(fit2_1,0.01)
tempOutput1 = topTable(fit2_1, coef=1, n=Inf)
DEG1= na.omit(tempOutput1)
deintrest1<-DEG1[intrest$mir,]
# 筛选差异基因
logFC_cutoff <- 1
DEG1$change = as.factor(
  ifelse(DEG1$P.Value < 0.05 & abs(DEG1$logFC) > logFC_cutoff,
         ifelse(DEG1$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
sig_diff1 <- subset(DEG1,
                    DEG1$P.Value < 0.05 & abs(DEG1$logFC) > logFC_cutoff)

dim(DEG1)
dim(sig_diff1)
# 60  6
summary(sig_diff1$change)
# DOWN  NOT   UP 
#  29    0   31 
write.table(group1,file = "group1.xls",
            quote = F,
            sep = "\t",
            row.names = T)
write.table(DEG1,file = "DEG_all1.xls",
            quote = F,
            sep = "\t",
            row.names = T)
write.table(sig_diff1,file = "DEG_sig1.xls",
            quote = F,
            sep = "\t",
            row.names = T)
## 火山图----------
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(Ipaper)
library(scales)
library(ggrepel)

dat_rep<-DEG1[rownames(DEG1)%in%
               c(head(rownames(subset(sig_diff1,sig_diff1$logFC>3)),10),
                 head(rownames(subset(sig_diff1,sig_diff1$logFC< -2.7)),10)),]
volcano_plot<- ggplot(data = DEG1, 
                      aes(x = logFC,
                          y = -log10(P.Value), 
                          color =change)) +
  scale_color_manual(values = c("blue", "darkgray","red")) +
  scale_x_continuous(breaks = c(-10,-6,-3,-1,0,1,3,6,10)) +
  scale_y_continuous(trans = "log1p",
                     breaks = c(0,1,5,10,20,50, 100,200)) +
  geom_point(size = 1.2, alpha = 0.4, na.rm=T) +
  theme_bw(base_size = 12, base_family = "Times") +
  geom_vline(xintercept = c(-1,1),
             lty = 4,
             col = "darkgray",
             lwd = 0.6)+
  geom_hline(yintercept = -log10(0.05),
             lty = 4,
             col = "darkgray",
             lwd = 0.6)+
  theme(legend.position = "right",
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(face="bold",
                                   color="black",
                                   family = "Times",
                                   size=13),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face = "bold",
                                   color = "black",
                                   size = 15),
        axis.text.y = element_text(face = "bold",
                                   color = "black",
                                   size = 15),
        axis.title.x = element_text(face = "bold",
                                    color = "black",
                                    size = 15),
        axis.title.y = element_text(face = "bold",
                                    color = "black",
                                    size = 15)) +
  geom_label_repel(
    data = dat_rep,
    aes(label = rownames(dat_rep)),
    max.overlaps = 20,
    size = 4,
    box.padding = unit(0.5, "lines"),
    min.segment.length = 0,
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )+
  labs(x = "log (Fold Change)",
       y = "-log10 (P.Value)")
volcano_plot
ggsave('volcano(t1 vs.c).png', volcano_plot,width = 9, height = 7)
ggsave('volcano(t1 vs.c).pdf', volcano_plot,width = 9, height = 7)
### 热图--------
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(gplots)
group_rt1<-group1$group
group_rt1<-as.data.frame(group_rt1)
rt1<-dat1
colnames(group_rt1)<-'group'
rownames(group_rt1)<-group1$sample
heat1<-rt1[rownames(rt1)%in%
             c(head(rownames(subset(sig_diff1,sig_diff1$logFC>3)),10),head(rownames(subset(sig_diff1,sig_diff1$logFC< -2.7)),10)),]
#x1<-log2(heat1+1)
x1<-t(scale(t(heat1)))
ann_colors<-list(
  Group = c(Normal="lightblue",Tumor="darkorange"))
pheatmap(mat=x1,
         annotation_col = group_rt1,
         color=bluered(100),
         scale = "row",
         annotation_colors = ann_colors,
         fontsize = 10,
         show_colnames = F,
         cluster_cols = F,
         cluster_rows = T)
## 02-2 test2 vs control------
dat2<-cbind(dat_test2,dat_control)
## 分组矩阵

group2<-data.frame(sample=colnames(dat2),
                   group=c(rep('test2',5),rep('control',5)))
type2<-group2[,2]
design2 <- model.matrix(~ -1+factor(type2,levels=c('control','test2'))) 
colnames(design2)<-c('control','test2')
rownames(design2)<-group2$sample
library(limma)
# 对每一个基因进行线性模型构建
fit2=lmFit(dat2,design2)
# 构建比较矩阵
contrast.matrix=makeContrasts(ControlVSMG=test2-control,levels = design2)
# 构建芯片数据线性模型，计算估计的相关系数和标准差
fit2_2=contrasts.fit(fit2,contrast.matrix)
# 基于贝叶斯计算T值，F值和log-odds
fit2_2<-eBayes(fit2_2,0.01)
tempOutput2 = topTable(fit2_2, coef=1, n=Inf)
DEG2= na.omit(tempOutput2)
deintrest2<-DEG2[rownames(DEG2)%in%intrest$mir,]
# 筛选差异基因
logFC_cutoff <- 1
DEG2$change = as.factor(
  ifelse(DEG2$P.Value < 0.05 & abs(DEG2$logFC) > logFC_cutoff,
         ifelse(DEG2$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
sig_diff2 <- subset(DEG2,
                    DEG2$P.Value < 0.05 & abs(DEG2$logFC) > logFC_cutoff)

dim(DEG2)
dim(sig_diff2)
# 55  6
summary(sig_diff2$change)
# DOWN  NOT   UP 
# 11    0   44 
write.table(group2,file = "group2.xls",
            quote = F,
            sep = "\t",
            row.names = T)
write.table(DEG2,file = "DEG_all2.xls",
            quote = F,
            sep = "\t",
            row.names = T)
write.table(sig_diff2,file = "DEG_sig2.xls",
            quote = F,
            sep = "\t",
            row.names = T)
## 火山图----------
dat_rep<-DEG2[rownames(DEG2)%in%
                c(head(rownames(subset(sig_diff2,sig_diff2$logFC>3.49)),10),
                  head(rownames(subset(sig_diff2,sig_diff2$logFC< -1.2)),10)),]
volcano_plot<- ggplot(data = DEG2, 
                      aes(x = logFC,
                          y = -log10(P.Value), 
                          color =change)) +
  scale_color_manual(values = c("blue", "darkgray","red")) +
  scale_x_continuous(breaks = c(-10,-6,-3,-1,0,1,3,6,10)) +
  scale_y_continuous(trans = "log1p",
                     breaks = c(0,1,5,10,20,50, 100,200)) +
  geom_point(size = 1.2, alpha = 0.4, na.rm=T) +
  theme_bw(base_size = 12, base_family = "Times") +
  geom_vline(xintercept = c(-1,1),
             lty = 4,
             col = "darkgray",
             lwd = 0.6)+
  geom_hline(yintercept = -log10(0.05),
             lty = 4,
             col = "darkgray",
             lwd = 0.6)+
  theme(legend.position = "right",
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(face="bold",
                                   color="black",
                                   family = "Times",
                                   size=13),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face = "bold",
                                   color = "black",
                                   size = 15),
        axis.text.y = element_text(face = "bold",
                                   color = "black",
                                   size = 15),
        axis.title.x = element_text(face = "bold",
                                    color = "black",
                                    size = 15),
        axis.title.y = element_text(face = "bold",
                                    color = "black",
                                    size = 15)) +
  geom_label_repel(
    data = dat_rep,
    aes(label = rownames(dat_rep)),
    max.overlaps = 20,
    size = 4,
    box.padding = unit(0.5, "lines"),
    min.segment.length = 0,
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )+
  labs(x = "log (Fold Change)",
       y = "-log10 (P.Value)")
volcano_plot
ggsave('volcano(t2 vs.c).png', volcano_plot,width = 9, height = 7)
ggsave('volcano(t2 vs.c).pdf', volcano_plot,width = 9, height = 7)
## 热图------
group_rt2<-group2$group
group_rt2<-as.data.frame(group_rt2)
rt2<-dat2
colnames(group_rt2)<-'group'
rownames(group_rt2)<-group2$sample
heat2<-rt2[rownames(rt2)%in%
             c(head(rownames(subset(sig_diff2,sig_diff2$logFC>3.49)),10),head(rownames(subset(sig_diff2,sig_diff2$logFC< -1.2)),10)),]
#x1<-log2(heat1+1)
x2<-t(scale(t(heat2)))
ann_colors<-list(
  Group = c(Normal="lightblue",Tumor="darkorange"))
pheatmap(mat=x2,
         annotation_col = group_rt2,
         color=bluered(100),
         scale = "row",
         annotation_colors = ann_colors,
         fontsize = 10,
         show_colnames = F,
         cluster_cols = F,
         cluster_rows = T)
## 02-1 test2 vs test1-----
dat3<-cbind(dat_test1,dat_test2)
## 分组矩阵

group3<-data.frame(sample=colnames(dat3),
                   group=c(rep('test1',5),rep('test2',5)))
type3<-group3[,2]
design3 <- model.matrix(~ -1+factor(type3,levels=c('test1','test2'))) 
colnames(design3)<-c('test1','test2')
rownames(design3)<-group3$sample
library(limma)
# 对每一个基因进行线性模型构建
fit3=lmFit(dat3,design3)
# 构建比较矩阵
contrast.matrix=makeContrasts(ControlVSMG=test2-test1,levels = design3)
# 构建芯片数据线性模型，计算估计的相关系数和标准差
fit2_3=contrasts.fit(fit3,contrast.matrix)
# 基于贝叶斯计算T值，F值和log-odds
fit2_3<-eBayes(fit2_3,0.01)
tempOutput3 = topTable(fit2_3, coef=1, n=Inf)
DEG3= na.omit(tempOutput3)
deintrest3<-DEG3[rownames(DEG3)%in%intrest$mir,]
# 筛选差异基因
logFC_cutoff <- 1
DEG3$change = as.factor(
  ifelse(DEG3$P.Value < 0.05 & abs(DEG3$logFC) > logFC_cutoff,
         ifelse(DEG3$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
sig_diff3 <- subset(DEG3,
                    DEG3$P.Value < 0.05 & abs(DEG3$logFC) > logFC_cutoff)

dim(DEG3)
dim(sig_diff3)
# 195   6
summary(sig_diff3$change)
# DOWN  NOT   UP 
# 54    0  141 
write.table(group3,file = "group3.xls",
            quote = F,
            sep = "\t",
            row.names = T)
write.table(DEG3,file = "DEG_all3.xls",
            quote = F,
            sep = "\t",
            row.names = T)
write.table(sig_diff3,file = "DEG_sig3.xls",
            quote = F,
            sep = "\t",
            row.names = T)
## 火山图----------
dat_rep<-DEG3[rownames(DEG3)%in%rownames(heat3),]
volcano_plot<- ggplot(data = DEG3, 
                      aes(x = logFC,
                          y = -log10(P.Value), 
                          color =change)) +
  scale_color_manual(values = c("blue", "darkgray","red")) +
  scale_x_continuous(breaks = c(-10,-6,-3,-1,0,1,3,6,10)) +
  scale_y_continuous(trans = "log1p",
                     breaks = c(0,1,5,10,20,50, 100,200)) +
  geom_point(size = 1.2, alpha = 0.4, na.rm=T) +
  theme_bw(base_size = 12, base_family = "Times") +
  geom_vline(xintercept = c(-1,1),
             lty = 4,
             col = "darkgray",
             lwd = 0.6)+
  geom_hline(yintercept = -log10(0.05),
             lty = 4,
             col = "darkgray",
             lwd = 0.6)+
  theme(legend.position = "right",
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(face="bold",
                                   color="black",
                                   family = "Times",
                                   size=13),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face = "bold",
                                   color = "black",
                                   size = 15),
        axis.text.y = element_text(face = "bold",
                                   color = "black",
                                   size = 15),
        axis.title.x = element_text(face = "bold",
                                    color = "black",
                                    size = 15),
        axis.title.y = element_text(face = "bold",
                                    color = "black",
                                    size = 15)) +
  geom_label_repel(
    data = dat_rep,
    aes(label = rownames(dat_rep)),
    max.overlaps = 20,
    size = 4,
    box.padding = unit(0.5, "lines"),
    min.segment.length = 0,
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )+
  labs(x = "log (Fold Change)",
       y = "-log10 (P.Value)")
volcano_plot
ggsave('volcano(t2 vs.t1).png', volcano_plot,width = 9, height = 7)
ggsave('volcano(t2 vs.t1).pdf', volcano_plot,width = 9, height = 7)

sig.all<-sig_diff1[rownames(sig_diff1)%in%rownames(sig_diff2),]
sig.all<-sig.all[rownames(sig.all)%in%rownames(sig_diff3),]
## 热图------
group_rt3<-group3$group
group_rt3<-as.data.frame(group_rt3)
rt3<-dat3
colnames(group_rt3)<-'group'
rownames(group_rt3)<-group3$sample
heat3<-rt3[rownames(rt3)%in%
             c(head(rownames(subset(sig_diff3,sig_diff3$logFC>4.2)),13),head(rownames(subset(sig_diff3,sig_diff3$logFC< -4)),10)),]
#x1<-log2(heat1+1)
x3<-t(scale(t(heat3)))
ann_colors<-list(
  Group = c(Normal="lightblue",Tumor="darkorange"))
pheatmap(mat=x3,
         annotation_col = group_rt3,
         color=bluered(100),
         scale = "row",
         annotation_colors = ann_colors,
         fontsize = 10,
         show_colnames = F,
         cluster_cols = F,
         cluster_rows = T)
# 03 DEmiRNA靶向mRNA的预测及凋亡-增殖相关mRNA（APmRNA）的筛选----------
setwd("/data/nas1/luchunlin/project/BJTC-304")
if (! dir.exists("./02_APmRNA")){
  dir.create("./02_APmRNA")
}
setwd("./02_APmRNA")
## 使用mirDIP、TargetScan、MiRTarBase数据库预测DEmiRNA的靶向mRNA，取交集.然后与增殖凋亡基因取交集
## 03-1 test1 vs control--------
### TargetScan
## 首先要对应mir家族
mir.family<-read_xlsx('miR_Family_Info.xlsx')
mir.family<-mir.family[,c(1,4)]
diff.family1<-mir.family[mir.family$`MiRBase ID`%in%rownames(sig_diff1),]
## 筛选靶基因
targetscan<-read_xlsx('Predicted_Targets_Info.default_predictions.xlsx')
targetscan1<-targetscan[targetscan$`miR Family`%in%diff.family1$`miR family`,]
colnames(diff.family1)<-c('miR Family','MiRBase ID')
targetscan1<-merge(diff.family1,targetscan1,by='miR Family')
targetscan1<-targetscan1[,c(2,4)]

write.table(targetscan1,file = 'targetscan(test1 vs.control).xls',
            sep = '\t',
            row.names = F)
#targetscan1<-targetscan1[!duplicated(targetscan1$`Gene Symbol`),]
### 4566
### mirDIP
#mirDIP1<-read_xlsx('mirDIP(test1 vs.control).xlsx')
#length(unique(mirDIP1$`Gene Symbol`))
## mirDIP1<-mirDIP1[!duplicated(mirDIP1$`Gene Symbol`),]
### 5396
### mirTarBase
mirTarbase1<-read.csv(file = 'mirTarbase(test1 vs. control).csv')
length(unique(mirTarbase1$Target))
## mirTarbase1<-mirTarbase1[!duplicated(mirTarbase1$Target),]
### 5551
### 3个数据库取交集 mirDIP不行
## target1<-mirTarbase1[mirTarbase1$ID%in%mirDIP1$MicroRNA,]

#target1<-intersect(x=targetscan1$`Gene Symbol`,y=mirDIP1$`Gene Symbol`)
#target1<-intersect(x=target1,y=mirTarbase1$Target)%>%as.data.frame()
### 1405
#colnames(target1)<-'Symbol'


library(ggvenn)
mydata1<-list(mirDIP=mirDIP1$`Gene Symbol`,TargetScan=targetscan1$`Gene Symbol`,mirTarbase=mirTarbase1$Target)
ggvenn(mydata1,c('mirDIP','TargetScan','mirTarbase'),
       fill_color = c("#ffb2b2","#b2e7cb","#F0E442"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb","#F0E442"),
       text_color = 'black')
### 与增殖凋亡基因集取交集
hallmark<-read_xlsx('proliferation-geneset(hallmark).xlsx')
hallmark<-hallmark[-1,]
colnames(hallmark)<-'geneset'
reactome<-read_xlsx('proliferation-geneset(reactome).xlsx')
reactome<-reactome[-1,]
colnames(reactome)<-'geneset'
geneset<-rbind(hallmark,reactome)
geneset<-geneset[!duplicated(geneset$geneset),]
## 308
APmRNA1<-data.frame(symbol=intersect(x=geneset$geneset,y=mirTarbase1$Target))
### 143
tarnet1<-mirTarbase1[mirTarbase1$Target%in%APmRNA1$symbol,]

write.table(APmRNA1,file = 'APmRNA(test1 vs.control).xls',
            sep = '\t',
            row.names = F)
write.table(tarnet1,file = 'APtarnet1.xls',
            sep = '\t',
            row.names = F)
APdata1<-list('Apoptosis'=geneset$geneset,'Target'=mirTarbase1$Target)
ggvenn(APdata1,c('Apoptosis','Target'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
## 03-2 test2 vs control--------
### TargetScan
## 首先要对应mir家族
diff.family2<-mir.family[mir.family$`MiRBase ID`%in%rownames(sig_diff2),]
## 筛选靶基因
targetscan2<-targetscan[targetscan$`miR Family`%in%diff.family2$`miR family`,]
colnames(diff.family2)<-c('miR Family','MiRBase ID')
targetscan2<-merge(diff.family2,targetscan2,by='miR Family')
targetscan2<-targetscan2[,c(2,4)]

write.table(targetscan2,file = 'targetscan(test2 vs.control).xls',
            sep = '\t',
            row.names = F)
length(unique(targetscan2$`Gene Symbol`))
#targetscan2<-targetscan2[!duplicated(targetscan2$`Gene Symbol`),]
### 3303
### mirDIP
mirDIP2<-read_xlsx('mirDIP(test2 vs.control).xlsx')
length(unique(mirDIP2$`Gene Symbol`))
##mirDIP2<-mirDIP2[!duplicated(mirDIP2$`Gene Symbol`),]
### 5084
### mirTarBase
mirTarbase2<-read.csv(file = 'mirTarbase(test2 vs.control).csv')
length(unique(mirTarbase2$Target))
#mirTarbase2<-mirTarbase2[!duplicated(mirTarbase2$Target),]
### 5423
### 3个数据库取交集
target2<-intersect(x=targetscan2$`Gene Symbol`,y=mirDIP2$`Gene Symbol`)
target2<-intersect(x=target2,y=mirTarbase2$Target)%>%as.data.frame()
### 1114
colnames(target2)<-'Symbol'

### 与增殖凋亡基因集取交集
APmRNA2<-data.frame(symbol=intersect(x=geneset$geneset,y=mirTarbase2$Target))
### 126
tarnet2<-mirTarbase2[mirTarbase2$Target%in%APmRNA2$symbol,]
write.table(APmRNA2,file = 'APmRNA(test2 vs.control).xls',
            sep = '\t',
            row.names = F)
write.table(tarnet2,file = 'APtarnet2.xls',
            sep = '\t',
            row.names = F)
APdata2<-list('Apoptosis'=geneset$geneset,'Target'=mirTarbase2$Target)
ggvenn(APdata2,c('Apoptosis','Target'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
## 03-3 test2 vs test1----------
### TargetScan
## 首先要对应mir家族
diff.family3<-mir.family[mir.family$`MiRBase ID`%in%rownames(sig_diff3),]
## 筛选靶基因
targetscan3<-targetscan[targetscan$`miR Family`%in%diff.family3$`miR family`,]
colnames(diff.family3)<-c('miR Family','MiRBase ID')
targetscan3<-merge(diff.family3,targetscan3,by='miR Family')
targetscan3<-targetscan3[,c(2,4)]

write.table(targetscan3,file = 'targetscan(test2 vs.test1).xls',
            sep = '\t',
            row.names = F)
length(unique(targetscan3$`Gene Symbol`))
#targetscan3<-targetscan3[!duplicated(targetscan3$`Gene Symbol`),]
### 5398
### mirDIP
mirDIP3<-read_xlsx('mirDIP(test2 vs.test1).xlsx')
length(unique(mirDIP3$`Gene Symbol`))
##mirDIP3<-mirDIP3[!duplicated(mirDIP3$`Gene Symbol`),]
### 9180
### mirTarBase
mirTarbase3<-read.csv(file = 'mirTarbase(test2 vs. test1).csv')
length(unique(mirTarbase3$Target))
##mirTarbase3<-mirTarbase3[!duplicated(mirTarbase3$Target),]
### 8291
### 3个数据库取交集
target3<-intersect(x=targetscan3$`Gene Symbol`,y=mirDIP3$`Gene Symbol`)
target3<-intersect(x=target3,y=mirTarbase3$Target)%>%as.data.frame()
### 2827
colnames(target3)<-'Symbol'

### 与增殖凋亡基因集取交集
APmRNA3<-data.frame(symbol=intersect(x=geneset$geneset,y=mirTarbase3$Target))
### 188
tarnet3<-mirTarbase3[mirTarbase3$Target%in%APmRNA3$symbol,]
write.table(APmRNA3,file = 'APmRNA(test2 vs.test1).xls',
            sep = '\t',
            row.names = F)
write.table(tarnet3,file = 'APtarnet3.xls',
            sep = '\t',
            row.names = F)
APdata3<-list('Apoptosis'=geneset$geneset,'Target'=mirTarbase3$Target)
ggvenn(APdata3,c('Apoptosis','Target'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')

# 04 lncRNA----------
setwd("/data/nas1/luchunlin/project/BJTC-304")
if (! dir.exists("./03_lncRNA")){
  dir.create("./03_lncRNA")
}
setwd("./03_lncRNA")
## 使用miRNet、starBase和lncbase数据库预测DEmiRNA的靶向lncRNA
## 04-1 test1 vs.control-------
## miRNet
mirnet1<-read.csv(file = 'mir2lnc(test1 vs.control).csv')
mirnet1<-mirnet1[,c(1,3)]
length(unique(mirnet1$Target))

## 353
## lncBase
#lncbase1<-read_xlsx('lncBase(test1 vs.control).xlsx')
#lncbase1<-lncbase1[,c(3,1)]
#colnames(lncbase1)<-c('ID','Target')
#length(unique(lncbase1$Target))
## 3726
## starBase
#starbase1<-read_xlsx('Starbase(test1 vs.control).xlsx')
#starbase1<-starbase1[,c(2,4)]
#length(unique(starbase1$geneName))
#colnames(starbase1)<-c('ID','Target')
## 989
#mir1<-lncbase1[lncbase1$ID%in%starbase1$ID,]
#lnc1<-mir1[mir1$Target%in%starbase1$Target,]
#length(unique(lnc1$Target))
### 364
##lncnet1<-starbase[starbase$geneName%in%lnc1,]
write.table(mirnet1,file = 'lncnet1.xls',
            sep = '\t',
            row.names = F)
## 04-2 test2 vs.control------
## miRNet
mirnet2<-read.csv(file = 'mir2lnc(test2 vs.control).csv')
length(unique(mirnet2$Target))
## 343
## lncBase
#lncbase2<-read_xlsx('lncBase(test2 vs.control).xlsx')
#lncbase2<-lncbase2[,c(3,1)]
#colnames(lncbase2)<-c('ID','Target')

## 4981
## starBase
#starbase2<-read_xlsx('Starbase(test2 vs.control).xlsx')
#length(unique(starbase2$geneName))
## 894
#starbase2<-starbase2[,c(2,4)]
#length(unique(starbase2$geneName))
#colnames(starbase2)<-c('ID','Target')
#mir2<-lncbase2[lncbase2$ID%in%starbase2$ID,]
#lnc2<-mir2[mir2$Target%in%starbase2$Target,]
#length(unique(lnc2$Target))
write.table(mirnet2,file = 'lncnet2.xls',
            sep = '\t',
            row.names = F)
## 04-3 test2 vs.test1--------

mirnet3<-read.csv(file = 'mir2lnc(test2 vs. test1).csv')
length(unique(mirnet3$Target))
## 600
#lncbase3<-read_xlsx('lncBase(test2 vs.test1).xlsx')
#lncbase3<-lncbase3[,c(3,1)]
#colnames(lncbase3)<-c('ID','Target')
#length(unique(lncbase3$Target))
## 3726
## starBase
#starbase3<-read_xlsx('Starbase(test2 vs.test1).xlsx')
#starbase3<-starbase3[,c(2,4)]
#length(unique(starbase3$geneName))
#colnames(starbase3)<-c('ID','Target')
## 989
#mir3<-lncbase3[lncbase3$ID%in%starbase3$ID,]
#lnc3<-mir3[mir3$Target%in%starbase3$Target,]
#length(unique(lnc3$Target))
### 717
##lncnet1<-starbase[starbase$geneName%in%lnc1,]
write.table(mirnet3,file = 'lncnet3.xls',
            sep = '\t',
            row.names = F)
# 05 network----------
setwd("/data/nas1/luchunlin/project/BJTC-304")
if (! dir.exists("./04_network")){
  dir.create("./04_network")
}
setwd("./04_network")
##使用DEmiRNA及其靶向的APmRNA、lncRNA构建mRNA-miRNA-lncRNA调控网
##degree top10 
## 05-1 test1 vs.control-------
## 核心miRNA
cemiRNA1<-c('hsa-let-7d-5p','hsa-miR-181a-5p','hsa-miR-181c-5p','hsa-miR-186-5p',
             'hsa-miR-193b-3p','hsa-miR-205-5p','hsa-miR-211-5p','hsa-miR-216a-5p','hsa-miR-4766-5p')
library(readxl)
cemiRNA1<-read_xlsx('hubmiRNA(test1 vs.control).xlsx')
## 核心mRNA
cemRNA1<-c('BCL2','BTG2','CDKN1A','HMGB1','MCL1','PMAIP1','SOD2','XIAP','YWHAZ')
cemRNA1<-read_xlsx('hubmRNA(test1 vs.control).xlsx')
## 核心lncRNA

celncRNA1<-c('HELLPAR','KCNQ1OT1','NEAT1','NORAD','XIST')
celncRNA1<-read_xlsx('hublncRNA(test1 vs.control).xlsx')

## top degree 直方图
cemiRNA1<-cemiRNA1[order(cemiRNA1$degree),]
cemiRNA1$group<-c(rep('A',3),rep('B',1),rep('C',1),rep('D',1),rep('E',1),rep('F',1),rep('G',1))
degree<-ggplot(data = cemiRNA1,aes(x=miRNAs,y=degree,fill=group))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values= c('#FFDDAA',"#FFBB66","#FFAA33",'#FF8800','#FF7744','#FF5511','#FF0000'))+
  coord_flip()+
 # ggtitle('Hub gene')+
  guides(fill='none')+
  theme(axis.text.x=element_text(hjust=1,colour="black",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",size=10))+
  labs(x='.')
degree
cemRNA1<-cemRNA1[order(cemRNA1$degree),]
cemRNA1$group<-c(rep('A',5),rep('B',2),rep('C',2))
degree<-ggplot(data = cemRNA1,aes(x=mRNA,y=degree,fill=group))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values= c('#FF8800','#FF7744','#FF5511'))+
  coord_flip()+
  # ggtitle('Hub gene')+
  guides(fill='none')+
  theme(axis.text.x=element_text(hjust=1,colour="black",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",size=10))+
  labs(x='.')
degree
celncRNA1<-celncRNA1[order(celncRNA1$degree),]
celncRNA1$group<-c(rep('A',2),rep('B',1),rep('C',1),rep('D',1))
degree<-ggplot(data = celncRNA1,aes(x=lncRNA,y=degree,fill=group))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values= c("#FFBB66",'#FF8800','#FF7744','#FF5511'))+
  coord_flip()+
  # ggtitle('Hub gene')+
  guides(fill='none')+
  theme(axis.text.x=element_text(hjust=1,colour="black",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",size=10))+
  labs(x='.')
degree

## 05-2 test2 vs.control-----
cemiRNA2<-c('hsa-let-7c-5p','hsa-let-7c-5p','hsa-miR-128-3p','hsa-miR-216a-5p','hsa-miR-5586-5p','hsa-miR-664b-3p','hsa-miR-761')
cemiRNA2<-read_xlsx('hubmiRNA(test2 vs.control).xlsx')
## 核心mRNA
cemRNA2<-c('BCL2L1','CCND1','CDKN1A','SOD2','XIAP','WEE1','YWHAZ')
cemRNA2<-read_xlsx('hubmRNA(test2 vs.control).xlsx')
## 核心lncRNA
celncRNA2<-c('KCNQ1OT1','MIR29B2CHG','NEAT1','XIST','SNHG16')
celncRNA2<-read_xlsx('hublncRNA(test2 vs.control).xlsx')

## top degree 直方图
cemiRNA2<-cemiRNA2[order(cemiRNA2$degree),]
cemiRNA2$group<-c(rep('A',1),rep('B',1),rep('C',1),rep('D',1),rep('E',1),rep('F',1),rep('G',1))
degree<-ggplot(data = cemiRNA2,aes(x=mirRNA,y=degree,fill=group))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values= c('#FFDDAA',"#FFBB66","#FFAA33",'#FF8800','#FF7744','#FF5511','#FF0000'))+
  coord_flip()+
  # ggtitle('Hub gene')+
  guides(fill='none')+
  theme(axis.text.x=element_text(hjust=1,colour="black",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",size=10))+
  labs(x='.')
degree
cemRNA2<-cemRNA2[order(cemRNA2$degree),]
cemRNA2$group<-c(rep('A',4),rep('B',2),rep('C',1))
degree<-ggplot(data = cemRNA2,aes(x=mRNA,y=degree,fill=group))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values= c('#FF8800','#FF7744','#FF5511'))+
  coord_flip()+
  # ggtitle('Hub gene')+
  guides(fill='none')+
  theme(axis.text.x=element_text(hjust=1,colour="black",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",size=10))+
  labs(x='.')
degree
celncRNA2<-celncRNA2[order(celncRNA2$degree),]
celncRNA2$group<-c(rep('A',2),rep('B',1),rep('C',1),rep('D',1))
degree<-ggplot(data = celncRNA1,aes(x=lncRNA,y=degree,fill=group))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values= c("#FFBB66",'#FF8800','#FF7744','#FF5511'))+
  coord_flip()+
  # ggtitle('Hub gene')+
  guides(fill='none')+
  theme(axis.text.x=element_text(hjust=1,colour="black",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",size=10))+
  labs(x='.')
degree


## 05-3 test2 vs.test1-----
cemiRNA3<-c('hsa-let-7c-5p','hsa-miR-128-3p','hsa-miR-181c-5p','hsa-miR-181d-5p','hsa-miR-216a-3p','hsa-miR-3612','hsa-miR-512-3p','hsa-miR-545-3p','hsa-miR-650')
cemiRNA3<-read_xlsx('hubmiRNA(test2 vs.test1).xlsx')
## 核心mRNA
cemRNA3<-c('BCL2L1','BTG2','CCND1','CDKN1A','MAPK1','SOD2','XIAP','WEE1','YWHAZ')
cemRNA3<-read_xlsx('hubmRNA(test2 vs.test1).xlsx')
## 核心lncRNA
celncRNA3<-c('KCNQ1OT1','LINC00963','OIP5-AS1','XIST','NEAT1')
celncRNA3<-read_xlsx('hublncRNA(test2 vs.test1).xlsx')

## top degree 直方图
cemiRNA3<-cemiRNA3[order(cemiRNA3$degree),]
cemiRNA3$group<-c(rep('A',2),rep('B',2),rep('C',1),rep('D',1),rep('E',1),rep('F',1),rep('G',2))
degree<-ggplot(data = cemiRNA3,aes(x=miRNA,y=degree,fill=group))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values= c('#FFDDAA',"#FFBB66","#FFAA33",'#FF8800','#FF7744','#FF5511','#FF0000'))+
  coord_flip()+
  # ggtitle('Hub gene')+
  guides(fill='none')+
  theme(axis.text.x=element_text(hjust=1,colour="black",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",size=10))+
  labs(x='.')
degree
cemRNA3<-cemRNA3[order(cemRNA3$degree),]
cemRNA3$group<-c(rep('A',2),rep('B',1),rep('C',2),rep('D',2),rep('E',1),rep('F',1))
degree<-ggplot(data = cemRNA3,aes(x=mRNA,y=degree,fill=group))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values= c('#FFDDAA',"#FFBB66","#FFAA33",'#FF8800','#FF7744','#FF5511'))+
  coord_flip()+
  # ggtitle('Hub gene')+
  guides(fill='none')+
  theme(axis.text.x=element_text(hjust=1,colour="black",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",size=10))+
  labs(x='.')
degree
celncRNA3<-celncRNA3[order(celncRNA3$degree),]
celncRNA3$group<-c(rep('A',1),rep('B',1),rep('C',1),rep('D',1),rep('E',1))
degree<-ggplot(data = celncRNA3,aes(x=lncRNA,y=degree,fill=group))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values= c('#FFDDAA',"#FFBB66",'#FF8800','#FF7744','#FF5511'))+
  coord_flip()+
  # ggtitle('Hub gene')+
  guides(fill='none')+
  theme(axis.text.x=element_text(hjust=1,colour="black",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",size=10))+
  labs(x='.')
degree
# 06 ceRNA----------
setwd("/data/nas1/luchunlin/project/BJTC-304")
if (! dir.exists("./04_hub")){
  dir.create("./04_hub")
}
setwd("./04_hub")

## 06-1 关键mRNA-miNRAs调控网络
hubgene<-c('BCL2L1','CDKN1A','MCL1','PMAIP1','XIAP')
## 06-2 GO/KEGG富集
## 06-2-1 test1 vs.control--------
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(GOplot)
gene_transform <- bitr(hubgene,
                       fromType = "SYMBOL",
                       toType = "ENTREZID",
                       OrgDb = "org.Hs.eg.db")
ego <- enrichGO(gene = gene_transform$ENTREZID,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
write.table(ego,file = "GO.xls",sep = "\t",quote = F,row.names = F)
# 展示富集最显著的 GO term
go_bar <- barplot(ego, showCategory=5, split="ONTOLOGY") +
  facet_grid(ONTOLOGY ~ ., scales = "free")
go_bar
ggsave('GO.png',go_bar,width = 8,height = 10)
ggsave('GO.pdf',go_bar,width = 8,height = 10)

## KEGG富集分析（气泡图）
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "hsa",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 2)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk@result,file = "KEGG.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk, showCategory=15)
kk_dot
ggsave('KEGG.png',kk_dot,width = 6,height = 6)
ggsave('KEGG.pdf',kk_dot,width = 6,height = 6)
## 06-2-2 test2 vs.control--------
gene_transform <- bitr(hubgene2,
                       fromType = "SYMBOL",
                       toType = "ENTREZID",
                       OrgDb = "org.Hs.eg.db")
ego <- enrichGO(gene = gene_transform$ENTREZID,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
write.table(ego,file = "GO(test2 vs.control).xls",sep = "\t",quote = F,row.names = F)
# 展示富集最显著的 GO term
go_bar <- barplot(ego, showCategory=5, split="ONTOLOGY") +
  facet_grid(ONTOLOGY ~ ., scales = "free")
go_bar
ggsave('GO(test2 vs.control).png',go_bar,width = 8,height = 10)
ggsave('GO(test2 vs.control).pdf',go_bar,width = 8,height = 10)

##KEGG富集分析（气泡图）
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "hsa",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 2)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk@result,file = "KEGG(test2 vs.control).xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk, showCategory=15)
kk_dot
ggsave('KEGG(test2 vs.control).png',kk_dot,width = 6,height = 6)
ggsave('KEGG(test2 vs.control).pdf',kk_dot,width = 6,height = 6)
## 06-2-3 test2 vs.test1----------
gene_transform <- bitr(hubgene3,
                       fromType = "SYMBOL",
                       toType = "ENTREZID",
                       OrgDb = "org.Hs.eg.db")
ego <- enrichGO(gene = gene_transform$ENTREZID,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
write.table(ego,file = "GO(test2 vs.test1).xls",sep = "\t",quote = F,row.names = F)
# 展示富集最显著的 GO term
go_bar <- barplot(ego, showCategory=5, split="ONTOLOGY") +
  facet_grid(ONTOLOGY ~ ., scales = "free")
go_bar
ggsave('GO(test2 vs.test1).png',go_bar,width = 8,height = 10)
ggsave('GO(test2 vs.test1).pdf',go_bar,width = 8,height = 10)

##KEGG富集分析（气泡图）
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "hsa",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 2)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk@result,file = "KEGG(test2 vs.test1).xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk, showCategory=15)
kk_dot
ggsave('KEGG(test2 vs.test1).png',kk_dot,width = 6,height = 6)
ggsave('KEGG(test2 vs.test1).pdf',kk_dot,width = 6,height = 6)
## 06-3 miRNA表达分析----------
setwd("/data/nas1/luchunlin/project/BJTC-304")
if (! dir.exists("./06_exp")){
  dir.create("./06_exp")
}
setwd("./06_exp")
## 06-3-1 test1 vs.control-------
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
## 提取表达矩阵
cemiRNA$miRNA<-gsub('mir','miR',cemiRNA$miRNA)
mir.exp1<-sig_diff1[rownames(sig_diff1)%in%cemiRNA$miRNA,]
mir.exp1<-dat1[rownames(dat1)%in%rownames(mir.exp1),]
##7个 
mir.exp1$Symbol<-rownames(mir.exp1)
mir.exp1<-gather(mir.exp1,
                 key = sample,
                 value = expr,
                 -c('Symbol'))
mir.exp1$Group<-ifelse(mir.exp1$sample%in%sample_test1,'test1','control')

##分面图形

# 分面图形
exp_boxplot<-ggplot(mir.exp1,aes(x = Symbol, y = expr, fill = Group)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "Expression",expand = c(0.1,0.1))+
  scale_x_discrete(name = "group") +
  ggtitle("Expression of miRNAs") +
  stat_compare_means(data = mir.exp1,
                     mapping = aes(group = Group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F) +
  theme_bw() +
  #  geom_signif(comparisons = list(c("control","POAG")),
  #              test = t.test,
  #              map_signif_level = T)+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top")
exp_boxplot

## 06-3-2 test2 vs.control-------
mir.exp2<-dat2[cemiRNA2,]
mir.exp2$Symbol<-rownames(mir.exp2)
mir.exp2<-gather(mir.exp2,
                 key = sample,
                 value = expr,
                 -c('Symbol'))
mir.exp2$Group<-ifelse(mir.exp2$sample%in%sample_test2,'test2','control')

##分面图形

# 分面图形
exp_boxplot<-ggplot(mir.exp2,aes(x = Symbol, y = expr, fill = Group)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "Expression",expand = c(0.1,0.1))+
  scale_x_discrete(name = "group") +
  ggtitle("Expression of miRNAs") +
  stat_compare_means(data = mir.exp2,
                     mapping = aes(group = Group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F) +
  theme_bw() +
  #  geom_signif(comparisons = list(c("control","POAG")),
  #              test = t.test,
  #              map_signif_level = T)+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top")
exp_boxplot

## 06-3-3 test2 vs.test1---------
mir.exp3<-dat3[cemiRNA3,]
mir.exp3$Symbol<-rownames(mir.exp3)
mir.exp3<-gather(mir.exp3,
                 key = sample,
                 value = expr,
                 -c('Symbol'))
mir.exp3$Group<-ifelse(mir.exp3$sample%in%sample_test2,'test2','test1')

# 分面图形
exp_boxplot<-ggplot(mir.exp3,aes(x = Symbol, y = expr, fill = Group)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "Expression",expand = c(0.1,0.1))+
  scale_x_discrete(name = "group") +
  ggtitle("Expression of miRNAs") +
  stat_compare_means(data = mir.exp3,
                     mapping = aes(group = Group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F) +
  theme_bw() +
  #  geom_signif(comparisons = list(c("control","POAG")),
  #              test = t.test,
  #              map_signif_level = T)+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top")
exp_boxplot
## 06-4 核心miRNA相关性分析--------
setwd("/data/nas1/luchunlin/project/BJTC-304")
if (! dir.exists("./07_cor")){
  dir.create("./07_cor")
}
setwd("./07_cor")
## 06-4-1 test1 vs.control-------
library(ggcorrplot)
library(corrplot)
## 提取表达矩阵
dat.exp1<-dat1[cemiRNA1,]
hub_corr1<-round(cor(t(dat.exp1)),3)
## 检验基因之间的相关性p值
## 计算相关性系数并显示基因之间的相关性。相关性系数大于0为正相关，小于0为负相关。
## p小于0.05认为相关性显著
hub_p.mat1<-round(cor_pmat(t(dat.exp1)),3)
col1 <- colorRampPalette(c("#4169E1","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF",
                           "#FDDBC7","#F4A582","#D6604D", "#B2182B","#CC0000","#990000"))
hub_corr_plot<-corrplot(hub_corr1,
                        method = "circle",
                        is.corr = T,
                        type = "lower",
                        p.mat = hub_p.mat1,
                        insig = "blank",
                        outline = "white",
                        addCoef.col ="black",
                        col = col1(200))
write.table(hub_corr1,file = 'hub_corr(test1 vs.control).xls',
            sep = '\t',
            row.names = T)
## 06-4-2 test2 vs.control------
## 提取表达矩阵
dat.exp2<-dat2[cemiRNA2,]
hub_corr2<-round(cor(t(dat.exp2)),3)
## 检验基因之间的相关性p值
## 计算相关性系数并显示基因之间的相关性。相关性系数大于0为正相关，小于0为负相关。
## p小于0.05认为相关性显著
hub_p.mat2<-round(cor_pmat(t(dat.exp2)),3)

hub_corr_plot<-corrplot(hub_corr2,
                        method = "circle",
                        is.corr = T,
                        type = "lower",
                        p.mat = hub_p.mat2,
                        insig = "blank",
                        outline = "white",
                        addCoef.col ="black",
                        col = col1(200))
write.table(hub_corr2,file = 'hub_corr(test2 vs.control).xls',
            sep = '\t',
            row.names = T)
## 06-4-3 test2 vs.test1-----
dat.exp3<-dat3[cemiRNA3,]
hub_corr3<-round(cor(t(dat.exp3)),3)
## 检验基因之间的相关性p值
## 计算相关性系数并显示基因之间的相关性。相关性系数大于0为正相关，小于0为负相关。
## p小于0.05认为相关性显著
hub_p.mat3<-round(cor_pmat(t(dat.exp3)),3)

hub_corr_plot<-corrplot(hub_corr3,
                        method = "circle",
                        is.corr = T,
                        type = "lower",
                        p.mat = hub_p.mat3,
                        insig = "blank",
                        outline = "white",
                        addCoef.col ="black",
                        col = col1(200))
write.table(hub_corr3,file = 'hub_corr(test2 vs.test1).xls',
            sep = '\t',
            row.names = T)
## 06-5 核心lncRNA亚细胞定位----------
## LNCipedia用于获取DElncRNA序列，lncLocator数据库用于识别DElncRNA的亚细胞定位。
setwd("/data/nas1/luchunlin/project/BJTC-304")
if (! dir.exists("./08_lncRNA")){
  dir.create("./08_lncRNA")
}
setwd("./08_lncRNA")
HELLPAR.loc<-data.frame(location=c('Cytoplasm','Nucleus','Ribosome','Cytosol','Exosome'),
                        score=c('0.890','0.023','0.018','0.067','0.001'))
### 直方图
colnames(HELLPAR.loc)
HELLPAR<-ggplot(data = HELLPAR.loc,aes(x=location,y=score,fill=location))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values= c('#FFAA33',"#FFA488",'#77DDFF','#99FF99','#D1BBFF'))+
  # coord_flip()+
  ggtitle('HELLPAR subcellular locations')+
  guides(fill='none')+
  theme(axis.text.x=element_text(hjust=0.5,colour="black",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",size=10))+
  labs(x='.')
HELLPAR
KCNQ1OT1.loc<-data.frame(location=c('Cytoplasm','Nucleus','Ribosome','Cytosol','Exosome'),
                        score=c('0.507','0.441','0.004','0.042','0.006'))
### 直方图
KCNQ1OT1<-ggplot(data = KCNQ1OT1.loc,aes(x=location,y=score,fill=location))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values= c('#FFAA33',"#FFA488",'#77DDFF','#99FF99','#D1BBFF'))+
  # coord_flip()+
  ggtitle('KCNQ1OT1 subcellular locations')+
  guides(fill='none')+
  theme(axis.text.x=element_text(hjust=0.5,colour="black",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",size=10))+
  labs(x='.')
KCNQ1OT1

OIP5_AS1.loc<-data.frame(location=c('Cytoplasm','Nucleus','Ribosome','Cytosol','Exosome'),
                         score=c('0.087','0.010','0.195','0.702','0.006'))
### 直方图
OIP5.AS1<-ggplot(data = OIP5_AS1.loc,aes(x=location,y=score,fill=location))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values= c('#FFAA33',"#FFA488",'#77DDFF','#99FF99','#D1BBFF'))+
  # coord_flip()+
  ggtitle('OIP5-AS1 subcellular locations')+
  guides(fill='none')+
  theme(axis.text.x=element_text(hjust=0.5,colour="black",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",size=10))+
  labs(x='.')
OIP5.AS1

MALAT1.loc<-data.frame(location=c('Cytoplasm','Nucleus','Ribosome','Cytosol','Exosome'),
                         score=c('0.716','0.190','0.013','0.031','0.049'))
### 直方图
MALAT1<-ggplot(data = MALAT1.loc,aes(x=location,y=score,fill=location))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values= c('#FFAA33',"#FFA488",'#77DDFF','#99FF99','#D1BBFF'))+
  # coord_flip()+
  ggtitle('MALAT1 subcellular locations')+
  guides(fill='none')+
  theme(axis.text.x=element_text(hjust=0.5,colour="black",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",size=10))+
  labs(x='.')
MALAT1

MIR29B2CHG.loc<-data.frame(location=c('Cytoplasm','Nucleus','Ribosome','Cytosol','Exosome'),
                       score=c('0.714','0.243','0.006','0.029','0.008'))
### 直方图
MIR29B2CHG<-ggplot(data = MIR29B2CHG.loc,aes(x=location,y=score,fill=location))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values= c('#FFAA33',"#FFA488",'#77DDFF','#99FF99','#D1BBFF'))+
  # coord_flip()+
  ggtitle('MIR29B2CHG subcellular locations')+
  guides(fill='none')+
  theme(axis.text.x=element_text(hjust=0.5,colour="black",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",size=10))+
  labs(x='.')
MIR29B2CHG

NEAT1.loc<-data.frame(location=c('Cytoplasm','Nucleus','Ribosome','Cytosol','Exosome'),
                       score=c('0.723','0.206','0.008','0.051','0.012'))
### 直方图
NEAT1<-ggplot(data = NEAT1.loc,aes(x=location,y=score,fill=location))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values= c('#FFAA33',"#FFA488",'#77DDFF','#99FF99','#D1BBFF'))+
  # coord_flip()+
  ggtitle('NEAT1 subcellular locations')+
  guides(fill='none')+
  theme(axis.text.x=element_text(hjust=0.5,colour="black",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",size=10))+
  labs(x='.')
NEAT1

NORAD.loc<-data.frame(location=c('Cytoplasm','Nucleus','Ribosome','Cytosol','Exosome'),
                      score=c('0.122','0.843','0.004','0.010','0.020'))

### 直方图
NORAD<-ggplot(data = NORAD.loc,aes(x=location,y=score,fill=location))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values= c('#FFAA33',"#FFA488",'#77DDFF','#99FF99','#D1BBFF'))+
  # coord_flip()+
  ggtitle('NORAD subcellular locations')+
  guides(fill='none')+
  theme(axis.text.x=element_text(hjust=0.5,colour="black",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",size=10))+
  labs(x='.')
NORAD

XIST.loc<-data.frame(location=c('Cytoplasm','Nucleus','Ribosome','Cytosol','Exosome'),
                        score=c('0.714','0.171','0.010','0.100','0.006'))
### 直方图
XIST<-ggplot(data = XIST.loc,aes(x=location,y=score,fill=location))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values= c('#FFAA33',"#FFA488",'#77DDFF','#99FF99','#D1BBFF'))+
  # coord_flip()+
  ggtitle('XIST subcellular locations')+
  guides(fill='none')+
  theme(axis.text.x=element_text(hjust=0.5,colour="black",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",size=10))+
  labs(x='.')
XIST
LINC00963.loc<-data.frame(location=c('Cytoplasm','Nucleus','Ribosome','Cytosol','Exosome'),
                     score=c('0.022','0.009','0.015','0.945','0.009'))

### 直方图
LINC00963<-ggplot(data = LINC00963.loc,aes(x=location,y=score,fill=location))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values= c('#FFAA33',"#FFA488",'#77DDFF','#99FF99','#D1BBFF'))+
  # coord_flip()+
  ggtitle('LINC00963 subcellular locations')+
  guides(fill='none')+
  theme(axis.text.x=element_text(hjust=0.5,colour="black",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",size=10))+
  labs(x='.')
LINC00963

SNHG16.loc<-data.frame(location=c('Cytoplasm','Nucleus','Ribosome','Cytosol','Exosome'),
                          score=c('0.203','0.086','0.099','0.541','0.070'))

### 直方图
SNHG16<-ggplot(data = SNHG16.loc,aes(x=location,y=score,fill=location))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values= c('#FFAA33',"#FFA488",'#77DDFF','#99FF99','#D1BBFF'))+
  # coord_flip()+
  ggtitle('SNHG16 subcellular locations')+
  guides(fill='none')+
  theme(axis.text.x=element_text(hjust=0.5,colour="black",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",size=10))+
  labs(x='.')
SNHG16
