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
dat1<-cbind(dat_test1,dat_control)

## 分组矩阵

group1<-data.frame(sample=colnames(dat1),
                   group=c(rep('test1',5),rep('control',5)))
type1<-group1[,2]
design1 <- model.matrix(~ -1+factor(type1,levels=c('control','test1'))) 
colnames(design1)<-c('control','test1')
rownames(design1)<-group1$sample
library(limma)
# 对每一个基因进行线性模型构建
fit1=lmFit(dat1,design1)
# 构建比较矩阵
contrast.matrix=makeContrasts(ControlVSMG=test1-control,levels = design1)
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
                  head(rownames(subset(sig_diff1,sig_diff1$logFC< -2.84)),9)),]
dat_rep<-rbind(deintrest1,dat_rep)
dat_rep<-dat_rep[-7,]
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
heat1<-rt1[rownames(rt1)%in%rownames(dat_rep),]
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

# 筛选差异基因
logFC_cutoff <- 1
DEG3$change = as.factor(
  ifelse(DEG3$P.Value < 0.05 & abs(DEG3$logFC) > logFC_cutoff,
         ifelse(DEG3$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
sig_diff3 <- subset(DEG3,
                    DEG3$P.Value < 0.05 & abs(DEG3$logFC) > logFC_cutoff)
deintrest3<-DEG3[rownames(DEG3)%in%intrest$mir,]
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
dat_rep<-DEG3[rownames(DEG3)%in%
                c(head(rownames(subset(sig_diff3,sig_diff3$logFC>4.4)),8),
                  head(rownames(subset(sig_diff3,sig_diff3$logFC< -4)),10)),]
dat_rep<-rbind(deintrest3,dat_rep)
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

## 热图------
group_rt3<-group3$group
group_rt3<-as.data.frame(group_rt3)
rt3<-dat3
colnames(group_rt3)<-'group'
rownames(group_rt3)<-group3$sample
heat3<-rt3[rownames(rt3)%in%rownames(dat_rep),]
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

### 取并集---------
library(tidyverse)
sig1<-rownames(sig_diff1)%>%as.data.frame()
sig2<-rownames(sig_diff2)%>%as.data.frame()
sig.all<-rbind(sig1,sig2)
sig.all<-sig.all[!duplicated(sig.all$.),]%>%as.data.frame()
colnames(sig.all)<-'miRNAs'

library(ggvenn)
mydata<-list('Test1 vs.control'=rownames(sig_diff1),'Test2 vs.control'=rownames(sig_diff2))
ggvenn(mydata,c('Test1 vs.control','Test2 vs.control'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')


## 使用MiRTarBase数据库预测DEmiRNA的靶向mRNA，
mirTarbase<-read.csv(file = 'mir2gene.csv')
length(unique(mirTarbase$Target))
##7783
## 然后与增殖凋亡基因取交集
library(readxl)
library(readr)
hallmark<-read_xlsx('proliferation-geneset(hallmark).xlsx')
hallmark<-hallmark[-1,]
colnames(hallmark)<-'geneset'
reactome<-read_xlsx('proliferation-geneset(reactome).xlsx')
reactome<-reactome[-1,]
colnames(reactome)<-'geneset'
geneset<-rbind(hallmark,reactome)%>%as.data.frame()
geneset<-geneset[!duplicated(geneset$geneset),]%>%as.data.frame()
colnames(geneset)<-'geneset'
## 308
APmRNA<-data.frame(symbol=intersect(x=geneset$geneset,y=mirTarbase$Target))
##178
tarnet<-mirTarbase[mirTarbase$Target%in%APmRNA$symbol,]
##474
write.table(APmRNA,file = 'APmRNA.xls',
            sep = '\t',
            row.names = F)
write.table(tarnet,file = 'APtarnet.xls',
            sep = '\t',
            row.names = F)
library(ggvenn)
APdata<-list('Apoptosis(hallmark)'=hallmark$geneset,'Target'=mirTarbase$Target,'Apoptosis(reactome)'=reactome$geneset)
ggvenn(APdata,c('Apoptosis(hallmark)','Target','Apoptosis(reactome)'),
       fill_color = c("#ffb2b2","#b2e7cb","#F0E442"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb","#F0E442"),
       text_color = 'black')

# 05 network----------
setwd("/data/nas1/luchunlin/project/BJTC-304")
if (! dir.exists("./03_network")){
  dir.create("./03_network")
}
setwd("./03_network")
#lncRNA----------
## 使用miRNet预测DEmiRNA的靶向lncRNA
mirnet<-read.csv(file = 'mir2lnc.csv')
mirnet<-mirnet[,c(1,3)]
length(unique(mirnet$Target))
##1480  483
##使用DEmiRNA及其靶向的APmRNA、lncRNA构建mRNA-miRNA-lncRNA调控网
##degree top10 
## 核心miRNA
cemiRNA<-read_xlsx('cemiRNA.xlsx')
## 核心mRNA
cemRNA<-read_xlsx('cemRNA.xlsx')
## 核心lncRNA
celncRNA<-read_xlsx('celncRNA.xlsx')

## top degree 直方图
cemiRNA<-cemiRNA[order(cemiRNA$degree),]
cemiRNA$group<-c(rep('A',1),rep('B',1),rep('C',1),rep('D',1),rep('E',1),rep('F',1),rep('G',1),rep('H',1),rep('I',1),rep('J',1))
cemiRNA$miRNA<-factor(cemiRNA$miRNA,levels = cemiRNA$miRNA)
degree<-ggplot(data = cemiRNA,aes(x=miRNA,y=degree,fill=group))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values= c('#FFDDAA',"#FFBB66",'#FFA488',"#FFAA33",'#FF8800','#FFA488','#FF7744','#FF5511','#FF3333','#FF0000'))+
  coord_flip()+
  # ggtitle('Hub gene')+
  guides(fill='none')+
  theme(axis.text.x=element_text(hjust=1,colour="black",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",size=10))+
  labs(x='.')
degree
cemRNA<-cemRNA[order(cemRNA$degree),]
cemRNA$group<-c(rep('A',2),rep('B',3),rep('C',3),rep('D',1))
cemRNA$mRNA<-factor(cemRNA$mRNA,levels = cemRNA$mRNA)
degree<-ggplot(data = cemRNA,aes(x=mRNA,y=degree,fill=group))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values= c("#FFAA33",'#FF8800','#FF7744','#FF5511'))+
  coord_flip()+
  # ggtitle('Hub gene')+
  guides(fill='none')+
  theme(axis.text.x=element_text(hjust=1,colour="black",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",size=10))+
  labs(x='.')
degree
celncRNA<-celncRNA[order(celncRNA$degree),]
celncRNA$group<-c(rep('A',2),rep('B',2),rep('C',1),rep('D',1),rep('E',2))
celncRNA$lncRNA<-factor(celncRNA$lncRNA,levels = celncRNA$lncRNA)
degree<-ggplot(data = celncRNA,aes(x=lncRNA,y=degree,fill=group))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values= c("#FFBB66","#FFAA33",'#FF8800','#FF7744','#FF5511'))+
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

## 06-3 miRNA表达分析----------
setwd("/data/nas1/luchunlin/project/BJTC-304")
if (! dir.exists("./05_exp")){
  dir.create("./05_exp")
}
setwd("./05_exp")

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
mir.exp2<-sig_diff2[rownames(sig_diff2)%in%cemiRNA$miRNA,]
mir.exp2<-dat2[rownames(dat2)%in%rownames(mir.exp2),]
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
                     method = 't.test',
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
if (! dir.exists("./06_cor")){
  dir.create("./06_cor")
}
setwd("./06_cor")
## 06-4-1 test1 vs.control-------
library(ggcorrplot)
library(corrplot)
## 提取表达矩阵
dat.exp1<-sig_diff1[rownames(sig_diff1)%in%cemiRNA$miRNA,]
dat.exp1<-dat1[rownames(dat1)%in%rownames(dat.exp1),]

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
dat.exp2<-sig_diff2[rownames(sig_diff2)%in%cemiRNA$miRNA,]
dat.exp2<-dat2[rownames(dat2)%in%rownames(dat.exp2),]
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
## 06-5 核心lncRNA亚细胞定位----------
## LNCipedia用于获取DElncRNA序列，lncLocator数据库用于识别DElncRNA的亚细胞定位。
setwd("/data/nas1/luchunlin/project/BJTC-304")
if (! dir.exists("./07_lncRNA")){
  dir.create("./07_lncRNA")
}
setwd("./07_lncRNA")
celncRNA
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