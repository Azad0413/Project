## rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/BJTC-282")
if (! dir.exists("./00_rawdata")){
  dir.create("./00_rawdata")
}
setwd("./00_rawdata")
##（1）GSE6477 (2)GSE7888（3）TCGA  MMRF-CoMMpass (预后用) (4)验证集GSE57317
library(GEOquery)
library(Biobase)
Sys.setenv("VROOM_CONNECTION_SIZE"=99999999)
## 01-1 GSE6477 -----
gset1<-getGEO("GSE6477",
              destdir = '.',
              GSEMatrix = T,
              getGPL = F)
expr1<-as.data.frame(exprs(gset1[[1]]))
a1=gset1[[1]]

gpl1<-getGEO("GPL96",destdir = '.')
a1=gset1[[1]]
gpl1<-Table(gpl1)    
colnames(gpl1)
probe2symobl1<-gpl1 %>%
  select('ID','Gene Symbol')%>%
  filter('Gene Symbol'!='')%>%
  separate('Gene Symbol',c('symbol','drop'),sep = '///')%>%
  select(-drop)
probe2symobl1=probe2symobl1[probe2symobl1$symbol!='',]

dat1<-expr1
dat1$ID<-rownames(dat1)
dat1$ID<-as.character(dat1$ID)
probe2symobl1$ID<-as.character(probe2symobl1$ID)
dat1<-dat1 %>%
  inner_join(probe2symobl1,by='ID')%>% 
  select(-ID)%>%     ## 去除多余信息
  select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
### 将需要的样本提取出来
pd1<-pData(a1)
pd1<-pd1[c(2,19)]
pd1<-pd1[order(pd1$description),]
pd1<-pd1[-c(1:21),]
pd1<-pd1[order(pd1$geo_accession),]
dat_final1<-dat1[,pd1$geo_accession]
dat_final1<-log2(dat_final1+1)
write.table(dat_final1,file = "exp1.xls",
            quote = F,
            sep = '\t',
            row.names = T)

##01-1-1 GSE109651---------(弃用)
gset1<-getGEO("GSE109651",
              destdir = '.',
              GSEMatrix = T,
              getGPL = F)
expr1<-as.data.frame(exprs(gset1[[1]]))
gpl1<-getGEO("GPL570",destdir = '.')
a1=gset1[[1]]
gpl1<-Table(gpl1)    
colnames(gpl1)
probe2symobl1<-gpl1 %>%
  select('ID','Gene Symbol')%>%
  filter('Gene Symbol'!='')%>%
  separate('Gene Symbol',c('symbol','drop'),sep = '///')%>%
  select(-drop)
probe2symobl1=probe2symobl1[probe2symobl1$symbol!='',]

dat1<-expr1
dat1$ID<-rownames(dat1)
dat1$ID<-as.character(dat1$ID)
probe2symobl1$ID<-as.character(probe2symobl1$ID)
dat1<-dat1 %>%
  inner_join(probe2symobl1,by='ID')%>% 
  select(-ID)%>%     ## 去除多余信息
  select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
## 将需要的样本提取出来  
pd1<-pData(a1)
dat_final1<-dat1
dat_final1<-log2(dat_final1+1)
write.table(dat_final1,file = "exp1.xls",
            quote = F,
            sep = '\t',
            row.names = T)

## 01-2 GSE7888-----
gset2<-getGEO("GSE7888",
              destdir = '.',
              GSEMatrix = T,
              getGPL = F)
expr2<-as.data.frame(exprs(gset2[[1]]))
gpl2<-getGEO("GPL570",destdir = '.')
a2=gset2[[1]]

gpl2<-Table(gpl2)    
colnames(gpl2)
probe2symobl2<-gpl2 %>%
  select('ID','Gene Symbol')%>%
  filter('Gene Symbol'!='')%>%
  separate('Gene Symbol',c('symbol','drop'),sep = '//')%>%
  select(-drop)
probe2symobl2=probe2symobl2[probe2symobl2$symbol!='',]

dat2<-expr2
dat2$ID<-rownames(dat2)
dat2$ID<-as.character(dat2$ID)
probe2symobl2$ID<-as.character(probe2symobl2$ID)
dat2<-dat2 %>%
  inner_join(probe2symobl2,by='ID')%>% 
  select(-ID)%>%     ## 去除多余信息
  select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
## 将需要的样本提取出来  6个early stage  6ge senescing stage
pd2<-pData(a2)
pd2<-pd2[,c(1,2)]
pd2<-pd2[-c(7:17),]
dat_final2<-dat2[,pd2$geo_accession]
dat_final2<-log2(dat_final2+1)
write.table(dat_final2,file = "exp2.xls",
            quote = F,
            sep = '\t',
            row.names = T)

 ## 01-3 MMRF-------
expr3<-read_tsv(file = 'MMRF-COMMPASS.htseq_fpkm.tsv')
expr3<-as.data.frame(expr3)
rownames(expr3)<-expr3[,1]
expr3<-expr3[,-1]
expr_counts<-read_tsv(file = 'MMRF-COMMPASS.htseq_counts.tsv') %>%as.data.frame()
rownames(expr_counts)<-expr_counts[,1]
#expr_counts<-expr_counts[,-1]
genecode<-read.table(file = 'gencode.v22.annotation.gene.probeMap')
probe2symbol3<-genecode[,(1:2)]
probe2symbol3<-probe2symbol3[-1,]
colnames(probe2symbol3)<-c('ID','symbol')
dat3<-expr3
dat3$ID<-rownames(dat3)
dat3$ID<-as.character(dat3$ID)
probe2symbol3$ID<-as.character(probe2symbol3$ID)
dat3<-dat3 %>%
  inner_join(probe2symbol3,by='ID')%>% 
  select(-ID)%>%     ## 去除多余信息
  select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除

#expr_counts$ID<-rownames(expr_counts)
#expr_counts$ID<-as.character(expr_counts$ID)
#probe2symbol3$ID<-as.character(probe2symbol3$ID)
#expr_counts<-expr_counts %>%
#  inner_join(probe2symbol3,by='ID')%>% 
#  select(-ID)%>%     ## 去除多余信息
#  select(symbol,everything())%>%     ## 重新排列
#  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
#  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
#  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
#  select(-rowMean)%>%     ## 反向选择去除rowMean这一列
#  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
#dataNorm.MMPF <- TCGAanalyze_Normalization(tabDF = expr_counts,
#                                           geneInfo = geneInfo,
#                                           method = "gcContent")
# 将标准化后的数据再过滤，得到最终的数据

#dataFilt.MMPF.final <- TCGAanalyze_Filtering(tabDF = dataNorm.MMPF,
#                                             method = "quantile", 
#                                             qnt.cut =  0.25)
#dim(dataFilt.MMPF.final)
dat_final3<-dat3
#dat_final3<-dat3[rownames(dat3)%in%rownames(dataFilt.MMPF.final),]

phenotype<-read_tsv(file = 'MMRF-COMMPASS.Xena_phenotype.tsv')
phenotype<-data.frame(sample=phenotype$samples.submitter_id,
                      gender=phenotype$demographic.gender,
                      race=phenotype$demographic.race,
                      stage=phenotype$diagnoses.tumor_stage)
#phenotype<-phenotype[phenotype$sample%in%colnames(dat_final3),]

survival<-read_tsv(file = 'MMRF-COMMPASS.survival.tsv')
survival<-survival[,-3]
survival<-survival[survival$sample%in%colnames(dat_final3),]

phenotype<-merge(survival,phenotype,by='sample')

## 01-4 GSE2658-------（弃用）
gset_va<-getGEO("GSE2658",
                destdir = '.',
                GSEMatrix = T,
                getGPL = F)
expr_va<-as.data.frame(exprs(gset_va[[1]]))
a_va=gset_va[[1]]
pd_va<-pData(a_va)
dat_va<-expr_va
dat_va$ID<-rownames(dat_va)
dat_va$ID<-as.character(dat_va$ID)
probe2symobl2$ID<-as.character(probe2symobl2$ID)
dat_va<-dat_va %>%
  inner_join(probe2symobl2,by='ID')%>% 
  select(-ID)%>%     ## 去除多余信息
  select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
dat_final_va<-dat_va
survival_va<-data.frame(sample=pd_va$geo_accession,
                        OS=pd_va$characteristics_ch1,
                        OS.time=pd_va$characteristics_ch1.2)
survival_va$OS.time
survival_va$OS<-gsub('[SURIND=','',survival_va$OS,fixed = T)
survival_va$OS<-gsub(' (Indicator of disease-related death; integer, 0=alive or death by other cause, 1=disease related death, na=death cause undetermined)]','',survival_va$OS,fixed = T)
survival_va$OS.time<-gsub('[SURTIM=','',survival_va$OS.time,fixed = T)
survival_va$OS.time<-gsub(' (Follow-up time in months from Pre-Treatment baseline; integer)]','',survival_va$OS.time,fixed = T)

survival_va$OS<-as.numeric(survival_va$OS)
survival_va$OS.time<-as.numeric(survival_va$OS.time)
survival_va$OS.time<-survival_va$OS.time*30

##01-4 验证集GSE57317--------
gset_va<-getGEO("GSE57317",
                destdir = '.',
                GSEMatrix = T,
                getGPL = F)
expr_va<-as.data.frame(exprs(gset_va[[1]]))
a_va=gset_va[[1]]
pd_va<-pData(a_va)
dat_va<-expr_va
dat_va$ID<-rownames(dat_va)
dat_va$ID<-as.character(dat_va$ID)
probe2symobl2$ID<-as.character(probe2symobl2$ID)
dat_va<-dat_va %>%
  inner_join(probe2symobl2,by='ID')%>% 
  select(-ID)%>%     ## 去除多余信息
  select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
dat_final_va<-dat_va
survival_va<-data.frame(sample=pd_va$geo_accession,
                         OS=as.numeric(pd_va$`os censored:ch1`),
                         OS.time=as.numeric(pd_va$`OS time:ch1`))

# 02 DEGs-----------
setwd("/data/nas1/luchunlin/project/BJTC-282")
if (! dir.exists("./01_DEGs")){
  dir.create("./01_DEGs")
}
setwd("./01_DEGs")
## 02-1 GSE6477差异分析----

## 分组矩阵

group1<-data.frame(sample=colnames(dat_final1),
                   group=c(rep('MM',129),rep('control',12)))
type1<-group1[,2]
design1 <- model.matrix(~ -1+factor(type1,levels=c('control','MM'))) 
colnames(design1)<-c('control','MM')
rownames(design1)<-group1$sample
library(limma)
# 对每一个基因进行线性模型构建
fit1=lmFit(dat_final1,design1)
# 构建比较矩阵
contrast.matrix=makeContrasts(ControlVSMG=MM-control,levels = design1)
# 构建芯片数据线性模型，计算估计的相关系数和标准差
fit2_1=contrasts.fit(fit1,contrast.matrix)
# 基于贝叶斯计算T值，F值和log-odds
fit2_1<-eBayes(fit2_1,0.01)
tempOutput1 = topTable(fit2_1, coef=1, n=Inf)
DEG1= na.omit(tempOutput1)

# 筛选差异基因
logFC_cutoff <- 0
DEG1$change = as.factor(
  ifelse(DEG1$adj.P.Val < 0.05 & abs(DEG1$logFC) > logFC_cutoff,
         ifelse(DEG1$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
sig_diff1 <- subset(DEG1,
                    DEG1$adj.P.Val < 0.05 & abs(DEG1$logFC) > logFC_cutoff)

dim(DEG1)
dim(sig_diff1)
# 3561    7
summary(sig_diff1$change)
# DOWN  NOT   UP 
# 1530    0 2031 
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
### 火山图---------
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(Ipaper)
library(scales)
library(ggrepel)
library(EnhancedVolcano)
dat_rep1<-DEG1[rownames(DEG1)%in%
               c(head(rownames(subset(sig_diff1,sig_diff1$logFC>2.33)),10),
                 head(rownames(subset(sig_diff1,sig_diff1$logFC< -4.6)),10)),]
volcano_plot<- ggplot(data = DEG1, 
                      aes(x = logFC,
                          y = -log10(adj.P.Val), 
                          color =change)) +
  scale_color_manual(values = c("blue", "darkgray","red")) +
  scale_x_continuous(breaks = c(-10,-6,-3,-1,0,1,3,6,10)) +
  scale_y_continuous(trans = "log1p",
                     breaks = c(0,1,5,10,20,50, 100,200)) +
  geom_point(size = 1.2, alpha = 0.4, na.rm=T) +
  theme_bw(base_size = 12, base_family = "Times") +
  #  geom_vline(xintercept = c(-1,1),
  #             lty = 4,
  #             col = "darkgray",
  #             lwd = 0.6)+
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
    data = dat_rep1,
    aes(label = rownames(dat_rep1)),
    max.overlaps = 20,
    size = 4,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE,
    min.segment.length = 0)+
  labs(x = "log (Fold Change)",
       y = "-log10 (adj.P.Val)")
volcano_plot
### 热图--------
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(gplots)
group_rt1<-group1$group
group_rt1<-as.data.frame(group_rt1)
rt1<-dat_final1
colnames(group_rt1)<-'group'
rownames(group_rt1)<-group1$sample
heat1<-rt1[rownames(rt1)%in%
           c(head(rownames(subset(sig_diff1,sig_diff1$logFC>1)),10),head(rownames(subset(sig_diff1,sig_diff1$logFC< -2)),10)),]
#x<-log2(heat+1)
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
## 02-2 GSE7888差异分析----
group2<-data.frame(sample=colnames(dat_final2),
                   group=c(rep('MSCs_early_stage',6),rep('MSCs_senescing_stage',6)))

type2<-group2[,2]
design2 <- model.matrix(~ -1+factor(type2,levels=c('MSCs_early_stage','MSCs_senescing_stage'))) 
colnames(design2)<-c('MSCs_early_stage','MSCs_senescing_stage')
rownames(design2)<-group2$sample

# 对每一个基因进行线性模型构建
fit2=lmFit(dat_final2,design2)
# 构建比较矩阵
contrast.matrix=makeContrasts(ControlVSMG=MSCs_senescing_stage-MSCs_early_stage,levels = design2)
# 构建芯片数据线性模型，计算估计的相关系数和标准差
fit2_2=contrasts.fit(fit2,contrast.matrix)
# 基于贝叶斯计算T值，F值和log-odds

fit2_2<-eBayes(fit2_2,0.01)
tempOutput2 = topTable(fit2_2, coef=1, n=Inf)
DEG2= na.omit(tempOutput2)

# 筛选差异基因
logFC_cutoff <- 0
DEG2$change = as.factor(
  ifelse(DEG2$adj.P.Val < 0.05 & abs(DEG2$logFC) > logFC_cutoff,
         ifelse(DEG2$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
sig_diff2 <- subset(DEG2,
                    DEG2$adj.P.Val < 0.05 & abs(DEG2$logFC) > logFC_cutoff)

dim(DEG2)
dim(sig_diff2)
# 217   7
summary(sig_diff2$change)
# DOWN  NOT   UP 
# 89    0  128
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
## 火山图-----
dat_rep2<-DEG2[rownames(DEG2)%in%
                 c(head(rownames(subset(sig_diff2,sig_diff2$logFC>3.7)),10),
                   head(rownames(subset(sig_diff2,sig_diff2$logFC< -2.6)),10)),]
volcano_plot2<- ggplot(data = DEG2, 
                      aes(x = logFC,
                          y = -log10(adj.P.Val), 
                          color =change)) +
  scale_color_manual(values = c("blue", "darkgray","red")) +
  scale_x_continuous(breaks = c(-10,-6,-3,-1,0,1,3,6,10)) +
  scale_y_continuous(trans = "log1p",
                     breaks = c(0,1,5,10,20,50, 100,200)) +
  geom_point(size = 1.2, alpha = 0.4, na.rm=T) +
  theme_bw(base_size = 12, base_family = "Times") +
  #  geom_vline(xintercept = c(-1,1),
  #             lty = 4,
  #             col = "darkgray",
  #             lwd = 0.6)+
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
    data = dat_rep2,
    aes(label = rownames(dat_rep2)),
    max.overlaps = 20,
    size = 4,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE,
    min.segment.length = 0)+
  labs(x = "log (Fold Change)",
       y = "-log10 (adj.P.Val)")
volcano_plot2
###热图------

group_rt2<-group2$group
group_rt2<-as.data.frame(group_rt2)
rt2<-dat_final2
colnames(group_rt2)<-'group'
rownames(group_rt2)<-group2$sample
heat2<-rt2[rownames(rt2)%in%
             c(head(rownames(subset(sig_diff2,sig_diff2$logFC>3)),10),head(rownames(subset(sig_diff2,sig_diff2$logFC< -2)),10)),]

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



## 02-2 GSE109651差异分析----（弃用）
group3<-data.frame(sample=colnames(dat_final1),
                   group=c(rep('Main_population',7),rep('Side_population',7)))

type3<-group3[,2]
design3 <- model.matrix(~ -1+factor(type3,levels=c('Main_population','Side_population'))) 
colnames(design3)<-c('Main_population','Side_population')
rownames(design3)<-group3$sample

# 对每一个基因进行线性模型构建
fit3=lmFit(dat_final1,design3)
# 构建比较矩阵
contrast.matrix=makeContrasts(ControlVSMG=Main_population-Side_population,levels = design3)
# 构建芯片数据线性模型，计算估计的相关系数和标准差
fit2_3=contrasts.fit(fit3,contrast.matrix)
# 基于贝叶斯计算T值，F值和log-odds
fit2_3<-eBayes(fit2_3,0.01)
tempOutput3 = topTable(fit2_3, coef=1, n=Inf)
DEG3= na.omit(tempOutput3)

# 筛选差异基因
logFC_cutoff <- 0
DEG3$change = as.factor(
  ifelse(DEG3$adj.P.Val < 0.05 & abs(DEG3$logFC) > logFC_cutoff,
         ifelse(DEG3$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
sig_diff3 <- subset(DEG3,
                    DEG3$adj.P.Val < 0.05 & abs(DEG3$logFC) > logFC_cutoff)

dim(DEG3)
dim(sig_diff3)
# 2062    7
summary(sig_diff3$change)
# DOWN  NOT   UP 
# 934    0 1128 
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

###火山图
dat_rep2<-DEG2[rownames(DEG2)%in%
                 c(head(rownames(subset(sig_diff2,sig_diff2$logFC>3)),10),
                   head(rownames(subset(sig_diff2,sig_diff2$logFC< -3)),10)),]

volcano_plot2<-EnhancedVolcano(DEG2,lab = rownames(DEG2),
                               x='logFC',
                               y='adj.P.Val',
                               pCutoff = 0.05,
                               FCcutoff = 1,
                               maxoverlapsConnectors = Inf,
                               legendLabels = c('NS','Log2FC','adj.P.Value',
                                                'adj.P.Value & LogFC'),
                               selectLab = c(rownames(dat_rep2)),
                               drawConnectors = T,
                               boxedLabels = T,
                               ylim = c(0,4)
)+
  labs(y = "-log10 (adj.P.Value)")
volcano_plot2

## 02-3 MM中的MSC衰老差异基因-------
## 取交集
sig_MM_MSC<-sig_diff1[rownames(sig_diff1)%in%rownames(sig_diff2),]
## 40
write.table(sig_MM_MSC,file = "sig_MM_MSC.xls",
            quote = F,
            sep = '\t',
            row.names = T)
## 


# 03 GO/KEGG富集----------
setwd("/data/nas1/luchunlin/project/BJTC-282")
if (! dir.exists("./02_GO_KEGG")){
  dir.create("./02_GO_KEGG")
}
setwd("./02_GO_KEGG")
diff_gene_names <- rownames(sig_MM_MSC)
gene_transform <- bitr(diff_gene_names,
                       fromType = "SYMBOL",
                       toType = c("ENTREZID"),
                       OrgDb = "org.Hs.eg.db")

ego <- enrichGO(gene = gene_transform$ENTREZID,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "ALL",
                pAdjustMethod = "none",
                pvalueCutoff = 0.05,
                readable = TRUE)
write.table(ego,file = "GO.xls",sep = "\t",quote = F,row.names = F)
# 展示富集最显著的 GO term
go_bar <- barplot(ego, showCategory=5, split="ONTOLOGY") +
  facet_grid(ONTOLOGY ~ ., scales = "free")
go_bar
##05-2 KEGG富集分析（气泡图）
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "hsa",
                 pAdjustMethod = "none",
                 pvalueCutoff = 0.05)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk@result,file = "KEGG.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk, showCategory=15)
kk_dot
# 04 PPI----------
setwd("/data/nas1/luchunlin/project/BJTC-282")
if (! dir.exists("./03_PPI")){
  dir.create("./03_PPI")
}
setwd("./03_PPI")
## top degree 直方图
top10degree<-read_xlsx('/data/nas1/luchunlin/project/BJTC-282/03_PPI/TOP DEGREE.xlsx')
colnames(top10degree)<-c('Symbol','Degree')
top10degree<-as.data.frame(top10degree)
top10degree<-top10degree[order(top10degree$Degree),]
top10degree$group<-c(rep('low',7),rep('mid',1),rep('high',1))
degree<-ggplot(data = top10degree,aes(x=Symbol,y=Degree,fill=group))+
  geom_bar(stat = 'identity')+
  scale_fill_manual(values= c('#FF4500',"#F4A460","#FF7F50"))+
  coord_flip()+
  ggtitle('Hub gene')+
  guides(fill='none')+
  theme(axis.text.x=element_text(hjust=1,colour="black",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",size=10))+
  labs(x='.')

degree

# 05 风险模型构建---------
## 05-1 单因素cox回归----
setwd("/data/nas1/luchunlin/project/BJTC-282")
if (! dir.exists("./04_univariate_cox")){
  dir.create("./04_univariate_cox")
}
setwd("./04_univariate_cox")
## 匹配生存数据
train_data<-t(dat_final3)
train_data<-as.data.frame(train_data)
train_data<-train_data[,colnames(train_data)%in%rownames(sig_MM_MSC)]
train_data$sample<-rownames(train_data)
train_data<-merge(survival,train_data,by='sample')
rownames(train_data)<-train_data$sample
train_data<-train_data[,-1]
### 单因素cox
library(survival)
library(survminer)
colnames_sum <- colnames(train_data)
colnames_sum <- gsub("-","_",colnames_sum)
colnames_sum <- gsub(" ","_",colnames_sum)
colnames(train_data) <- colnames_sum

covariates <- colnames_sum[-which(colnames_sum %in% c("OS", "OS.time"))]
#Surv()函数产生一个生存对象  生存时间对生存的影响 对每一个变量构建生存分析公式
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste("Surv(OS.time, OS)~", x)))  #as.formula(). 将字符串转换成公式。构建formula对象
# coxph函数用于计算cox模型 循环对每一个特征做cox回归分析
univ_models <- lapply(univ_formulas,
                      function(x) {coxph(x, data = train_data)})

univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         #获取p值
                         p.value<-signif(x$wald["pvalue"], digits=3)
                         #获取HR
                         HR <-signif(x$coef[2], digits=3);
                         #获取95%置信区间
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"],3)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],3)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", 
                                      HR.confint.upper, ")")
                         res<-c(p.value,HR)
                         names(res)<-c("p.value","HR (95% CI for HR)")
                         return(res)
                       })

## coef是公式中的回归系数b（有时候也叫beta值）。exp(coef)是cox模型中的风险比（HR）
## z代表wald统计量，是coef除以其标准误se(coef)。ower .95 upper .95则是exp(coef)的95%置信区间，可信区间越窄，可信度越高，你的实验越精确，越是真理。
res_mod <- t(as.data.frame(univ_results, check.names = FALSE))
res_mod <- as.data.frame(res_mod)
res_results_0.05 <- res_mod[which(as.numeric(res_mod$p.value) < 0.05),]
res_results_0.05 <- na.omit(res_results_0.05)
write.table(res_results_0.05,
            file = "univariate_cox_result_0.05.xls",
            quote = F,
            sep = '\t',
            row.names = T)
dim(res_results_0.05)
## 16 2
library(tidyr)
res_results_0.05_2 <- separate(res_results_0.05, "HR (95% CI for HR)",
                               into = c("HR", "HR.95L", "HR.95H"),
                               sep = " ")
res_results_0.05_2 <- separate(res_results_0.05_2, "HR.95L",
                               into = c("HR.95L", "HR.95H"),
                               sep = "\\-")
res_results_0.05_2$HR.95L <- gsub("\\(", "", res_results_0.05_2$HR.95L)
res_results_0.05_2$HR.95H <- gsub("\\)", "", res_results_0.05_2$HR.95H)

res_results_0.05_2[,1:ncol(res_results_0.05_2)] <- as.numeric(unlist(res_results_0.05_2[,1:ncol(res_results_0.05_2)]))
res_results_0.05_2 <- res_results_0.05_2[order(res_results_0.05_2$HR),]
hz <- paste(round(res_results_0.05_2$HR,3),
            "(",round(res_results_0.05_2$HR.95L,3),
            "-",round(res_results_0.05_2$HR.95H,3),")",sep = "")


tabletext <- cbind(c(NA,"Gene",rownames(res_results_0.05_2)),
                   c(NA,"P value",ifelse(res_results_0.05_2$p.value<0.001,
                                         "< 0.001",
                                         round(res_results_0.05_2$p.value,3))),
                   c(NA,"Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1
library(forestplot)

pdf(file = "univariate_cox_forest.pdf", height = 10, width = 10, onefile = F)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE, TRUE,rep(FALSE, 70)),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,NA,res_results_0.05_2$HR),
           lower=c(NA,NA,res_results_0.05_2$HR.95L), #95%置信区间下限
           upper=c(NA,NA,res_results_0.05_2$HR.95H), #95%置信区间上限
           boxsize=0.2,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(0, 1,2,3,4,5,6,7,8,9), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(1.2,"cm"), #固定行高
           graphwidth = unit(.5,"npc"), #图在表中的宽度比例
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines = list("3" = gpar(col = "black", lty = 1, lwd = 2)),
           # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
           #                 "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
           #                 "59" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
           # mar=unit(rep(0.01, times = 4), "cm"),#图形页边距
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Hazard Ratio",
           grid = T) # 垂直于x轴的网格线，对应每个刻度
dev.off()
png(filename = "univariate_cox_forest.png", height = 700, width = 700)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE, TRUE,rep(FALSE, 70)),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,NA,res_results_0.05_2$HR),
           lower=c(NA,NA,res_results_0.05_2$HR.95L), #95%置信区间下限
           upper=c(NA,NA,res_results_0.05_2$HR.95H), #95%置信区间上限
           boxsize=0.2,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(0, 1,2,3,4,5,6,7,8,9), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(1.2,"cm"), #固定行高
           graphwidth = unit(.5,"npc"), #图在表中的宽度比例
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines = list("3" = gpar(col = "black", lty = 1, lwd = 2)),
           # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
           #                 "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
           #                 "59" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
           # mar=unit(rep(0.01, times = 4), "cm"),#图形页边距
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Hazard Ratio",
           grid = T) # 垂直于x轴的网格线，对应每个刻度

dev.off()

## 05-2 Lasso回归-----
setwd("/data/nas1/luchunlin/project/BJTC-282")
if (! dir.exists("./05_Lasso")){
  dir.create("./05_Lasso")
}
setwd("./05_Lasso")
library(glmnet)
x_all <- subset(train_data, select = -c(OS, OS.time))
x_all <- x_all[,rownames(res_results_0.05)]
y_all <- subset(train_data, select = c(OS, OS.time))

# 拟合模型
fit <- glmnet(as.matrix(x_all), Surv(y_all$OS.time,y_all$OS), 
              family = "cox") 
#dev.new()
png(filename = "lasso_model.png", height = 400, width = 500)
plot(fit, xvar = "lambda",label = TRUE, las=1)
dev.off()
pdf(file = "lasso_model.pdf", height = 5)
plot(fit, xvar = "lambda",label = TRUE, las=1)
dev.off()
# 交叉验证拟合模型
set.seed(9)
cvfit = cv.glmnet(as.matrix(x_all),
                  Surv(y_all$OS.time,y_all$OS),nfold=50,
                  family = "cox") 

png(filename = "lasso_verify.png", height = 400, width = 500)
plot(cvfit, las =1)
dev.off()
pdf(file = "lasso_verify.pdf", height = 5)
plot(cvfit, las =1)
dev.off()

# 提取指定lambda时特征的系数
coef.min = coef(cvfit, s = "lambda.min")  ## lambda.min & lambda.1se 取一个
#coef.min
#coef.min<-as.data.frame(as.matrix(coef.min))[lasso_geneids$lasso_geneids,]
#rownames(coef.min)<-lasso_geneids$lasso_geneids
#class(coef.min)
cvfit$lambda.min
# [1] 0.008623835
# 找出那些回归系数没有被惩罚为0的
active.min = which(coef.min@i != 0)

# 提取基因名称
lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i+1]
lasso_geneids
#coef.min<-as.data.frame(coef.min)
## 14
## [1] "COBLL1"  "CCND1"   "MANSC1"  "MAN2A1"  "EGR1"    "SLC31A2" "HSPB8"   "SCN3A"   "SCN9A"   "SOX11"   "RGS7"    "NTF3"    "OCLN"    "ZBED1"  
write(lasso_geneids, "lasso_genes.csv")
write.csv(x_all,file = "Lasso_x.csv",quote = F)
write.csv(y_all,file = "Lasso_y.csv",quote = F)

# 06 风险模型的构建与验证------
setwd("/data/nas1/luchunlin/project/BJTC-282")
if (! dir.exists("./06_risk")){
  dir.create("./06_risk")
}
setwd("./06_risk")
#lasso_geneids<-as.data.frame(lasso_geneids)

## 06-1 计算每个患者的风险评分，展示生存状态分布-----
##手动算
#train_data2<-train_data[,lasso_geneids$lasso_geneids]
#risk<-data.frame(train_data2)
#risk$riskscore<-NA
#risk$risk<-NA
#cnt<-1
#while (cnt < 859) {
#  risk$riskscore[cnt]<-sum(coef.min*train_data2[cnt,])
#  cnt = cnt + 1
#}

#dim(train_data2)
#cnt<-1
#while (cnt < 859) {
#  risk$risk[cnt]=as.vector(ifelse( risk$riskscore[cnt]>median(risk$riskscore),0,1))
#  cnt = cnt + 1
#}

riskScore=predict(cvfit,newx = as.matrix(x_all),s=cvfit$lambda.min)
coef.min

riskScore<-as.numeric(riskScore)
class(riskScore)
coxGene=lasso_geneids
outCol=c("OS","OS.time",coxGene)
risk=as.vector(ifelse(riskScore>median(riskScore),0,1))
risk <- as.data.frame(c(cbind(id=rownames(cbind(train_data[,outCol],
                                                riskScore,
                                                risk)),
                              cbind(train_data[,outCol],
                                    riskScore,
                                    risk))))
table(risk$risk)
#0  1 
#429 429 
library(ggplot2)
library(ggthemes)
library(Ipaper)
median(riskScore)
# [1] -1.86919
risk_dis <- ggplot(risk, aes(x=reorder(id, riskScore), 
                             y=riskScore, 
                             color = factor(risk, 
                                            levels = c(0, 1), 
                                            labels = c("High Risk", "Low Risk")))) +
  geom_point() +
  scale_color_manual(values = c("#A73030FF", "#0073C2FF")) + 
  scale_x_discrete(breaks = risk[order(risk$riskScore),]$id[c(1,100,200,300,400,500,600,700,800,900)],
                   labels = c(1,100,200,300,400,500,600,700,800,900),
                   expand = c(0.02,0)) +
  geom_vline(xintercept = nrow(risk[which(risk$risk==1),]) + 0.5,
             lty = 2) +
  geom_hline(yintercept = median(riskScore),
             lty =2) +
  labs(x = "Patients(increasing risk score)",
       y = "Risk Score",
       title = "Train Risk Score Distribution") + 
  theme_base() +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.justification = c(0,1),
        legend.position = c(0,1),
        #        legend.margin = margin(c(-5,4,4,3)),
        legend.background = element_rect(color = "black", size = .3),
        plot.title = element_text(size = 15, hjust = 0.5))
risk_dis
surv_stat <- ggplot(risk, aes(x=reorder(id, riskScore),
                              y=OS.time/365,
                              color = factor(OS,
                                             levels = c(0,1),
                                             labels = c("Alive", "Dead")))) +
  geom_point() +
  scale_color_manual(values = c("#0073C2FF","#A73030FF")) +
  scale_x_discrete(breaks = risk[order(risk$riskScore),]$id[c(1,100,200,300,400,500,600,700,800,900)],
                   labels = c(1,100,200,300,400,500,600,700,800,900),
                   expand = c(0.02,0)) +
  ylim(c(0,10))+
  geom_vline(xintercept = nrow(risk[which(risk$risk==1),]) + 0.5,
             lty = 2) +
  labs(x = "Patients(increasing risk score)",
       y = "Progression-free Interval (days)",
       title = "Train survival state distribution") + 
  theme_base() +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.justification = c(0,1),
        legend.position = c(0,1),
        #        legend.margin = margin(c(-5,4,4,3)),
        legend.background = element_rect(color = "black", size = .3),
        plot.title = element_text(size = 15, hjust = 0.5))

surv_stat

## 06-2 KM曲线和ROC曲线------
kmfit<-survfit(Surv(OS.time, OS) ~ risk, data =  risk)
train_survival_median <- ggsurvplot(kmfit,
                                    pval = TRUE, 
                                    conf.int = F,
                                    legend.labs=c("High risk","Low risk" ),
                                    legend.title="Risk score",
                                    title="Train KM",
                                    font.main = c(15,"bold"),
                                    risk.table = TRUE, 
                                    risk.table.col = "strata", 
                                    linetype = "strata", 
                                    surv.median.line = "hv", 
                                    ggtheme = theme_bw(), 
                                    palette = c("#A73030FF", "#0073C2FF"))
train_survival_median

riskscore <- function(survival_cancer_df,
                      candidate_genes_for_cox,
                      cox_report){
  library("dplyr")
  risk_score_table <- survival_cancer_df[, candidate_genes_for_cox]
  for (each_sig_gene in colnames(risk_score_table)){
    risk_score_table$each_sig_gene <- risk_score_table[,each_sig_gene]*
      (summary(cox_report)$coefficients[each_sig_gene,1])
  }
  risk_score_table <- cbind(risk_score_table,
                            "total_risk_score"=exp(rowSums(risk_score_table))) %>%
    cbind(survival_cancer_df[,c("OS.time", "OS")])
  risk_score_table <- risk_score_table[,c("OS.time",
                                          "OS",
                                          candidate_genes_for_cox,
                                          "total_risk_score")]
}
candidate_genes_for_cox2 <- c(rownames(cox_table)[cox_table[,3]>0.05])
risk_score_table_multi_cox2 <- riskscore(train_data,
                                         candidate_genes_for_cox2,
                                         cox_more_2)

multi_ROC <- function(time_vector, risk_score_table){
  library(survivalROC)
  single_ROC <- function(single_time){
    for_ROC <- survivalROC(Stime=risk_score_table$OS.time,
                           status=risk_score_table$OS,
                           marker=risk_score_table$riskScore,
                           predict.time=single_time,
                           method = 'KM')
    data.frame('True_positive'=for_ROC$TP, 'False_positive'=for_ROC$FP,
               'Cut_values'=for_ROC$cut.values, 'Time'=rep(single_time, length(for_ROC$TP)),
               'AUC'=rep(for_ROC$AUC, length(for_ROC$TP)))
  }
  multi_ROC_list <- lapply(time_vector, single_ROC)
  do.call(rbind, multi_ROC_list)
}
for_multi_ROC <- multi_ROC(time_vector = c(365*seq(1,5,2)), 
                           risk_score_table = risk)
for_multi_ROC$Time <- factor(for_multi_ROC$Time)
# 画ROC曲线 
library(scales)
library(geomROC)
library(plotROC)

auc_y1 <- round(for_multi_ROC[which(for_multi_ROC$Time==365),5][1],2)
auc_y3 <- round(for_multi_ROC[which(for_multi_ROC$Time==1095),5][1],2)
auc_y5 <- round(for_multi_ROC[which(for_multi_ROC$Time==1825),5][1],2)


ROC <- ggplot(for_multi_ROC, aes(x=False_positive,
                                 y=True_positive, 
                                 label=Cut_values, 
                                 color=Time)) + 
  scale_color_manual(breaks = c("365", "1095", "1825"),
                     labels = c("1 years", "3 years", "5 years"),
                     values = c("#00468b", "#A73030FF", "#42b540")) +
  geom_roc(labels = F, stat = 'identity') + 
  style_roc() + 
  geom_abline(slope = 1, intercept = 0, color = 'gray', linetype=2) +
  theme_bw() +
  labs(title = "Train ROC") +
  theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
        panel.grid = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold")) +
  # annotate('text', x=.75, y=.25, label=paste('AUC of 1 years =', round(auc_y1,2))) + 
  # annotate('text', x=.75, y=.15, label=paste('AUC of 3 years =', round(auc_y3,2))) + 
  # annotate('text', x=.75, y=.05, label=paste('AUC of 5 years =', round(auc_y5,2))) +
  annotate("text", x=0.75, y=c(0.25, 0.15, 0.05),
           label = c(paste('AUC of 1 years =', format(auc_y1,nsmall=2)),
                     paste('AUC of 3 years =', format(auc_y3,nsmall=2)),
                     paste('AUC of 5 years =', format(auc_y5,nsmall=2))))
ROC
train_roc<-roc(risk$OS,risk$riskScore,
               levels=c(0,1),)

plot(train_roc,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     print.auc.pattern = "AUC=%.2f",
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="Train ROC",
     col="#FF2E63",
     legacy.axes=F)
## 06-3 外部数据库验证（KM，年ROC）---------
setwd("/data/nas1/luchunlin/project/BJTC-282")
if (! dir.exists("./07_External_va")){
  dir.create("./07_External_va")
}
setwd("./07_External_va")

test_data<-t(dat_final_va)%>%as.data.frame()
#test_data<-t(scale(t(test_data)))%>%as.data.frame()
#test_data<-log2(test_data+1)
test_data$sample<-rownames(test_data)
test_data<-merge(survival_va,test_data,by='sample')
rownames(test_data)<-test_data$sample
test_data<-test_data[,-1]

test_data<-test_data[,colnames(test_data)%in%colnames(train_data)]


# 开始验证
x_all_out <- subset(test_data, select = -c(OS, OS.time))
x_all_out <- x_all_out[,rownames(res_results_0.05)]
riskScore_out=predict(cvfit,newx = as.matrix(x_all_out),s=cvfit$lambda.min)
riskScore_out<-as.numeric(riskScore_out)

risk_out=as.vector(ifelse(riskScore_out>median(riskScore_out),0,1))
risk_out <- as.data.frame(c(cbind(id=rownames(cbind(test_data[,outCol],
                                                    riskScore_out,
                                                    risk_out)),
                                  cbind(test_data[,outCol],
                                        riskScore_out,
                                        risk_out))))

library(ggplot2)
library(ggthemes)
median(riskScore_out)

risk_dis_out <- ggplot(risk_out, aes(x=reorder(id, riskScore_out), 
                                     y=riskScore_out, 
                                     color = factor(risk_out, 
                                                    levels = c(0, 1), 
                                                    labels = c("High Risk", "Low Risk")))) +
  geom_point() +
  scale_color_manual(values = c("#A73030FF", "#0073C2FF")) + 
  scale_x_discrete(breaks = risk_out[order(risk_out$riskScore_out),]$id[c(1,10,20,30,40,50,60,70,80)],
                   labels = c(1,10,20,30,40,50,60,70,80),
                   expand = c(0.02,0)) +
  geom_vline(xintercept = nrow(risk_out[which(risk_out$risk_out==1),]) + 0.5,
             lty = 2) +
  geom_hline(yintercept = median(riskScore_out),
             lty =2) +
  labs(x = "Patients(increasing risk score)",
       y = "Risk Score",
       title = "Validation Risk Score Distribution") + 
  theme_base() +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.justification = c(0,1),
        legend.position = c(0,1),
        #        legend.margin = margin(c(-5,4,4,3)),
        legend.background = element_rect(color = "black", size = .3),
        plot.title = element_text(size = 15, hjust = 0.5))
risk_dis_out
surv_stat_out <- ggplot(risk_out, aes(x=reorder(id, riskScore_out),
                                      y=OS.time/365,
                                      color = factor(OS,
                                                     levels = c(0,1),
                                                     labels = c("Alive", "Dead")))) +
  geom_point() +
  scale_color_manual(values = c("#0073C2FF", "#A73030FF")) +
  scale_x_discrete(breaks = risk_out[order(risk_out$riskScore_out),]$id[c(1,10,20,30,40,50,60,70,80)],
                   labels = c(1,10,20,30,40,50,60,70,80),
                   expand = c(0.02,0)) +
  ylim(x=c(0,7)) +
  geom_vline(xintercept = nrow(risk_out[which(risk_out$risk_out==1),]) + 0.5,
             lty = 2) +
  labs(x = "Patients(increasing risk score)",
       y = "Survival time (years)",
       title = "Validation Survival State Distribution") + 
  theme_base() +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.justification = c(0,1),
        legend.position = c(0,1),
        #        legend.margin = margin(c(-5,4,4,3)),
        legend.background = element_rect(color = "black", size = .3),
        plot.title = element_text(size = 15, hjust = 0.5))
surv_stat_out
kmfit_out <- survfit(Surv(OS.time, OS) ~ risk_out, data =  risk_out)
verify_survival_median <- ggsurvplot(kmfit_out,
                                     pval = TRUE, 
                                     conf.int = F,
                                     legend.labs=c("High risk","Low risk" ),
                                     legend.title="Risk score",
                                     title="Validition KM",
                                     font.main = c(15,"bold"),
                                     risk.table = TRUE, 
                                     risk.table.col = "strata", 
                                     linetype = "strata", 
                                     surv.median.line = "hv", 
                                     ggtheme = theme_bw(), 
                                     palette = c("#A73030FF", "#0073C2FF"))
verify_survival_median


riskscore <- function(survival_cancer_df,
                      candidate_genes_for_cox,
                      cox_report){
  library("dplyr")
  risk_score_table <- survival_cancer_df[, candidate_genes_for_cox]
  for (each_sig_gene in colnames(risk_score_table)){
    risk_score_table$each_sig_gene <- risk_score_table[,each_sig_gene]*
      (summary(cox_report)$coefficients[each_sig_gene,1])
  }
  risk_score_table <- cbind(risk_score_table,
                            "total_risk_score"=exp(rowSums(risk_score_table))) %>%
    cbind(survival_cancer_df[,c("OS.time", "OS")])
  risk_score_table <- risk_score_table[,c("OS.time",
                                          "OS",
                                          candidate_genes_for_cox,
                                          "total_risk_score")]
}
candidate_genes_for_cox2 <- c(rownames(cox_table)[cox_table[,3]>0.05])
risk_score_table_multi_cox2 <- riskscore(test_data,
                                         candidate_genes_for_cox2,
                                         cox_more_2)

multi_ROC <- function(time_vector, risk_score_table){
  library(survivalROC)   
  single_ROC <- function(single_time){
    for_ROC <- survivalROC(Stime=risk_score_table$OS.time,
                           status=risk_score_table$OS,
                           marker=risk_score_table$riskScore_out,
                           predict.time=single_time,method = 'KM')
    data.frame('True_positive'=for_ROC$TP, 'False_positive'=for_ROC$FP,
               'Cut_values'=for_ROC$cut.values, 'Time'=rep(single_time, length(for_ROC$TP)),
               'AUC'=rep(for_ROC$AUC, length(for_ROC$TP)))
  }
  multi_ROC_list <- lapply(time_vector, single_ROC)
  do.call(rbind, multi_ROC_list)
}

for_multi_ROC <- multi_ROC(time_vector = c(365*seq(1,5,2)), 
                           risk_score_table = risk_out)
for_multi_ROC$Time <- factor(for_multi_ROC$Time)
# 画ROC曲线 
library(scales)
library(geomROC)
library(plotROC)
auc_y1 <- round(for_multi_ROC[which(for_multi_ROC$Time==365),5][1],2)
auc_y3 <- round(for_multi_ROC[which(for_multi_ROC$Time==1095),5][1],2)
auc_y5 <- round(for_multi_ROC[which(for_multi_ROC$Time==1825),5][1],2)


ROC <- ggplot(for_multi_ROC, aes(x=False_positive,
                                 y=True_positive, 
                                 label=Cut_values, 
                                 color=Time)) + 
  scale_color_manual(breaks = c("365", "1095", "1825"),
                     labels = c("1 years", "3 years", "5 years"),
                     values = c("#00468b", "#A73030FF", "#42b540")) +
  geom_roc(labels = F, stat = 'identity') + 
  style_roc() + 
  geom_abline(slope = 1, intercept = 0, color = 'gray', linetype=2) +
  theme_bw() +
  labs(title = "Train ROC") +
  theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
        panel.grid = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold")) +
  # annotate('text', x=.75, y=.25, label=paste('AUC of 1 years =', round(auc_y1,2))) + 
  # annotate('text', x=.75, y=.15, label=paste('AUC of 3 years =', round(auc_y3,2))) + 
  # annotate('text', x=.75, y=.05, label=paste('AUC of 5 years =', round(auc_y5,2))) +
  annotate("text", x=0.75, y=c(0.25, 0.15, 0.05),
           label = c(paste('AUC of 1 years =', format(auc_y1,nsmall=2)),
                     paste('AUC of 3 years =', format(auc_y3,nsmall=2)),
                     paste('AUC of 5 years =', format(auc_y5,nsmall=2))))
ROC

test_roc<-roc(risk_out$OS,risk_out$riskScore_out,
         levels=c(0,1))

plot(test_roc,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     print.auc.pattern = "AUC=%.2f",
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="Test ROC",
     col="#FF2E63",
     legacy.axes=F)
# 07 预后模型的构建与评价--------
setwd("/data/nas1/luchunlin/project/BJTC-282")
if (! dir.exists("./12_prog_model")){
  dir.create("./12_prog_model")
}
setwd("./12_prog_model")

## 07-1 单因素Cox----------
train_phenotype<-phenotype
train_phenotype$OS<-as.numeric(train_phenotype$OS)
train_phenotype$OS.time<-as.numeric(train_phenotype$OS.time)
train_phenotype$gender<-ifelse(train_phenotype$gender=='male',1,2)
train_phenotype$stage<-gsub('III','3',train_phenotype$stage)
train_phenotype$stage<-gsub('II','2',train_phenotype$stage)
train_phenotype$stage<-gsub('I','1',train_phenotype$stage)
train_phenotype$stage<-gsub('Unknown',NA,train_phenotype$stage)
train_phenotype$stage<-gsub('unknown',NA,train_phenotype$stage)
train_phenotype$race<-gsub('not reported',NA,train_phenotype$race)
train_phenotype2$race<-gsub('NA',NA,train_phenotype2$race)
train_phenotype2$race<-gsub('not allowed to collect',NA,train_phenotype2$race)
train_phenotype2$race<-gsub('other',NA,train_phenotype2$race)
train_phenotype$race<-ifelse(train_phenotype$race=='asian',1,ifelse(train_phenotype$race=='white',2,3))
colnames(train_phenotype)<-c('id','OS','OS.time','gender','race','stage')

sub_risk <- subset(risk, select = c(id, riskScore))

train_risk_clinical <- merge(train_phenotype,
                             sub_risk,
                             by = "id")
rownames(train_risk_clinical) <- train_risk_clinical$id
train_risk_clinical = subset(train_risk_clinical, select = -c(id))
dim(train_risk_clinical)
colnames_train <- colnames(train_risk_clinical)
covariates_train <- colnames_train[-which(colnames_train %in% c("OS", "OS.time"))]

train_risk_clinical$stage<-factor(train_risk_clinical$stage)
train_risk_clinical$race<-factor(train_risk_clinical$race)
library(survival)
res.risk = coxph(Surv(time = OS.time, event = OS) ~ riskScore, data = train_risk_clinical) %>% summary
res.risk = c(res.risk$conf.int[-2], res.risk$coefficients[5])
res.race = coxph(Surv(time = OS.time, event = OS) ~ race, data = train_risk_clinical) %>% summary
res.race = cbind(res.race$conf.int[,-2], res.race$coefficients[,5])
res.gender = coxph(Surv(time = OS.time, event = OS) ~ gender, data = train_risk_clinical) %>% summary
res.gender = c(res.gender$conf.int[-2], res.gender$coefficients[5])
res.stage = coxph(Surv(time = OS.time, event = OS) ~ stage, data = train_risk_clinical) %>% summary
res.stage = cbind(res.stage$conf.int[,-2], res.stage$coefficients[,5])

res.ref = c(1,1,1,NA)
res = rbind(res.risk, res.gender,res.ref,res.race, res.ref, res.stage) %>% as.data.frame()
rownames(res)

res$Indicators = c("riskScore","Gender","asian(Reference)","white","black or african american","Stage1(Reference)","Stage2","Stage3")
colnames(res) = c("hr","low","up","pv","Indicator")
res$p = signif(res$pv, 2) %>% paste0("p = ", .)
res$p[is.na(res$pv)] = NA
res$Indicator = factor(res$Indicator, levels = rev(res$Indicator))
rownames(res) <- res$Indicator
res2 <- data.frame(p.value=res$pv,
                   HR=res$hr,
                   HR.95L=res$low,
                   HR.95H=res$up,
                   Indicator=res$Indicator)
rownames(res2) <- res2$Indicator
write.table(res2, file = "univariate_cox_prog_forest.xls", sep = "\t", quote = F)
res2 <- subset(res2, select = -c(Indicator))
library(tidyr)
hz <- paste(round(res2$HR,3),
            "(",round(res2$HR.95L,3),
            "-",round(res2$HR.95H,3),")",sep = "")
hz
hz[c(3,6)] <- ""

tabletext <- cbind(c(NA,rownames(res2)),
                   c("P value",ifelse(res2$p.value<0.0001,
                                      "< 0.0001",
                                      round(res2$p.value,4))),
                   c("Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1
library(forestplot)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,TRUE,FALSE,rep(FALSE, 7)),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,res2$HR),
           lower=c(NA,res2$HR.95L), #95%置信区间下限
           upper=c(NA,res2$HR.95H), #95%置信区间上限
           boxsize=0.1,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(0,5,10,15,20,25,30,35), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(1.2,"cm"), #固定行高
           graphwidth = unit(.6,"npc"), #图在表中的宽度比例
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),
           # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
           #                 "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
           #                 "59" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
           # mar=unit(rep(0.01, times = 4), "cm"),#图形页边距
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Hazard Ratio",
           grid = T,
           title = "Univariate") # 垂直于x轴的网格线，对应每个刻度
## 07-2 多因素Cox----------
res.mul = coxph(Surv(time = OS.time, event = OS) ~ riskScore + gender + stage, data = train_risk_clinical)%>% summary
res.mul = cbind(res.mul$conf.int[,-2], res.mul$coefficients[,5]) %>% as.data.frame()
res.mul = rbind(res.mul[1:2,], res.ref, res.mul[c(3:4),])
res.mul$Indicators = c("riskScore","Gender","stage1(Reference)", "stage2","stage3")
colnames(res.mul) = c("hr","low","up","pv","Indicator")
res.mul$p = signif(res.mul$pv, 2) %>% paste0("p = ", .)
res.mul$p[is.na(res.mul$pv)] = NA
res.mul$Indicator = factor(res.mul$Indicator, levels = rev(res.mul$Indicator))
rownames(res.mul) <- res.mul$Indicator

multi_res <- data.frame(p.value=res.mul$pv,
                        HR=res.mul$hr,
                        HR.95L=res.mul$low,
                        HR.95H=res.mul$up,
                        Indicator=res.mul$Indicator)
rownames(multi_res) <- multi_res$Indicator
multi_res
write.csv(multi_res,
          file = "multivariate_cox_prog_result.csv",
          quote = F,
          row.names = T)
multi_res <- subset(multi_res, select = -c(Indicator))

library(tidyr)
hz <- paste(round(multi_res$HR,3),
            "(",round(multi_res$HR.95L,3),
            "-",round(multi_res$HR.95H,3),")",sep = "")
hz

hz[c(3)] <- ""
tabletext <- cbind(c(NA,rownames(multi_res)),
                   c("P value",ifelse(multi_res$p.value<0.0001,
                                      "< 0.0001",
                                      round(multi_res$p.value,4))),
                   c("Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,multi_res$HR),
           lower=c(NA,multi_res$HR.95L), #95%置信区间下限
           upper=c(NA,multi_res$HR.95H), #95%置信区间上限
           boxsize=0.1,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(0,0.5,1,1.5,2,2.5,3,3.5,4,4.5,5,5.5), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(1.2,"cm"), #固定行高
           graphwidth = unit(.6,"npc"), #图在表中的宽度比例
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),
           # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
           #                 "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
           #                 "59" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
           # mar=unit(rep(0.01, times = 4), "cm"),#图形页边距
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Hazard Ratio",
           grid = T,
           title = "Multivariate") # 垂直于x轴的网格线，对应每个刻度
##07-3 构建COX模型，绘制列线图---------
multi_cov<-c('riskScore',"stage")
cox_data_prog <- as.formula(paste0('Surv(OS.time, OS)~',
                                   paste(multi_cov,
                                         sep = '',
                                         collapse = '+')))
cox_more_prog <- coxph(cox_data_prog,
                  data = as.data.frame(train_risk_clinical))

# Nomogram
library(rms)
ddist <- datadist(train_risk_clinical)
options(datadist='ddist')

# 构建COX模型，绘制列线图

res.cox <- psm(cox_data_prog,
               data = train_risk_clinical, dist = 'lognormal')
surv <- Survival(res.cox) # 构建生存概率函数
function(x) surv(365, x) # 1年事件发生概率
function(x) surv(1095, x) # 3年事件发生概率
function(x) surv(1825, x) # 5年事件发生概率

nom.cox <- nomogram(res.cox,
                    fun = list(function(x) surv(365, x),
                               function(x) surv(1095, x),
                               function(x) surv(1825, x)),
                    funlabel=c("1-year Survival Probability", "3-year Survival Probability", "5-year Survival Probability"),
                    maxscale = 10,
                    fun.at = c(0.01,seq(0.1,0.9,by=0.2),0.95,0.99),
                    lp=F)


plot(nom.cox, cex.axis  = 1.5, cex.var = 1.7)
png(filename = "nomogram_line_points.png", height = 600, width = 1200)
plot(nom.cox, cex.axis  = 1.5, cex.var = 1.7)
dev.off()
pdf(file = "nomogram_line_points.pdf", height = 8, width = 17)
plot(nom.cox, cex.axis  = 1.5, cex.var = 1.7)
dev.off()
##07-4 构建校准曲线---------

##绘制1年生存期校准曲线
coxm_1 <- cph(cox_data_prog,
              data=train_risk_clinical,
              surv=T,
              x=T,
              y=T,
              time.inc = 365)
cal_1<-calibrate(coxm_1,u=365,cmethod='KM',m=100,B=1000)

##绘制3年生存期校曲线
##time.in 和 u 要是一样的，都是要评价的时间节点
coxm_3 <- cph(cox_data_prog,
              data=train_risk_clinical,
              surv=T,
              x=T,y=T,
              time.inc = 3*365)
cal_3 <-calibrate(coxm_3,u=3*365,cmethod='KM',m=200,B=1000)

coxm_5 <- cph(cox_data_prog,
              data=train_risk_clinical,
              surv=T,
              x=T,y=T,
              time.inc = 5*365)
cal_5 <-calibrate(coxm_5,u=5*365,cmethod='KM',m=100)

par(mar=c(7,4,4,3),cex=1.5)
plot(cal_1,
     subtitles = F,
     lwd=2,lty=1, ##设置线条形状和尺寸
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 1-5 year Progression-free Interval',#便签
     ylab='Actual 1-5 year Progression-free Interval(Proportion)',#标签
     col="#00468b",#设置一个颜色
     xlim = c(0,1),ylim = c(0,1)) ##x轴和y轴范围
plot(cal_3,
     add = T,
     subtitles = F,
     lwd=2,lty=1,  ##设置线条宽度和线条类型
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 1-5 year Progression-free Interval',#便签
     ylab='Actual 1-5 year Progression-free Interval(Proportion)',#标签
     col="#ed0000",#设置一个颜色
     xlim = c(0,1),ylim = c(0,1)) ##x轴和y轴范围
plot(cal_5,
     add = T,
     subtitles = F,
     lwd=2,lty=1, ##设置线条形状和尺寸
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 1-5year Progression-free Interval',#便签
     ylab='Actual 1-year Progression-free Interval(Proportion)',#标签
     col="#42b540",#设置一个颜色
     xlim = c(0,1),ylim = c(0,1)) ##x轴和y轴范围



#加上图例
legend("bottomright", legend=c("1-year", "3-year", "5-year"), 
       col=c("#00468b", "#ed0000", "#42b540"), 
       lwd=2)
#调整对角线
abline(0,1,lty=5,lwd=2,col="grey")

## 07-5 KM----------
res.mul = coxph(Surv(time = OS.time, event = OS)~ riskScore + gender + stage, data = train_risk_clinical)
train_risk_clinical$riskScore_prog = predict(res.mul, newdata = train_risk_clinical, type = "lp")

train_risk_clinical$risk_prog = ifelse(train_risk_clinical$riskScore_prog > median(train_risk_clinical$riskScore_prog, na.rm = T), "high", "low")
train_risk_clinical$risk_prog = factor(train_risk_clinical$risk_prog, levels = c("high", "low"), labels = c("High risk", "Low risk"))
surv.fit = survfit(Surv(time = OS.time,event = OS) ~ risk_prog, data = train_risk_clinical)
prog_survival_median <- ggsurvplot(surv.fit,
                                   pval = TRUE, 
                                   conf.int = F,
                                   legend.labs=c("High risk","Low risk" ),
                                   legend.title="Risk score",
                                   title="Prognosis Model KM",
                                   font.main = c(15,"bold"),
                                   risk.table = TRUE, 
                                   risk.table.col = "strata", 
                                   linetype = "strata", 
                                   surv.median.line = "hv",
                                   ggtheme = theme_bw(), 
                                   palette = c("#A73030FF", "#0073C2FF"))
prog_survival_median

## 07-6 ROC--------
library(survivalROC)
library(tidyverse)
###计算每个患者的风险评分，展示生存状态分布
prog_roc<-roc(train_risk_clinical2$OS,train_risk_clinical2$riskScore_prog,
                levels=c(0,1))

plot(prog_roc,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     print.auc.pattern = "AUC=%.2f",
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="Prognosis Model ROC",
     col="#FF2E63",
     legacy.axes=F)
# 开始验证
train_risk_clinical2 <- train_risk_clinical

library(survival)
library(survminer)

riskscore <- function(survival_cancer_df,
                      candidate_genes_for_cox,
                      cox_report){
  library("dplyr")
  risk_score_table <- survival_cancer_df[, candidate_genes_for_cox]
  for (each_sig_gene in colnames(risk_score_table)){
    risk_score_table$each_sig_gene <- risk_score_table[,each_sig_gene]*
      (summary(cox_report)$coefficients[each_sig_gene,1])
  }
  risk_score_table <- cbind(risk_score_table,
                            "total_risk_score"=exp(rowSums(risk_score_table))) %>%
    cbind(survival_cancer_df[,c("OS.time", "OS")])
  risk_score_table <- risk_score_table[,c("OS.time",
                                          "OS",
                                          candidate_genes_for_cox,
                                          "total_risk_score")]
}


# candidate_genes_for_cox2 <- c(rownames(cox_table)[cox_table[,3]>0.05])

risk_score_table_multi_cox2<-train_risk_clinical2[,-10]
multi_ROC_out <- function(time_vector, risk_score_table){
  library(survivalROC)
  single_ROC <- function(single_time){
    for_ROC <- survivalROC(Stime=risk_score_table$OS.time,
                           status=risk_score_table$OS,
                           marker=risk_score_table$riskScore_prog,
                           predict.time=single_time,method = 'KM')
    data.frame('True_positive'=for_ROC$TP, 'False_positive'=for_ROC$FP,
               'Cut_values'=for_ROC$cut.values, 'Time'=rep(single_time, length(for_ROC$TP)),
               'AUC'=rep(for_ROC$AUC, length(for_ROC$TP)))
  }
  multi_ROC_list <- lapply(time_vector, single_ROC)
  do.call(rbind, multi_ROC_list)
}

for_multi_ROC <- multi_ROC_out(time_vector = c(365*seq(1,5,2)), 
                               risk_score_table = train_risk_clinical2)
for_multi_ROC$Time <- factor(for_multi_ROC$Time)

# 画ROC曲线 
library(scales)
library(geomROC)
library(plotROC)
auc_y1 <- round(for_multi_ROC[which(for_multi_ROC$Time==365),5][1],2)
auc_y3 <- round(for_multi_ROC[which(for_multi_ROC$Time==1095),5][1],2)
auc_y5 <- round(for_multi_ROC[which(for_multi_ROC$Time==1825),5][1],2)

ROC <- ggplot(for_multi_ROC, aes(x=False_positive,
                                 y=True_positive, 
                                 label=Cut_values, 
                                 color=Time)) + 
  scale_color_manual(breaks = c("365", "1095", "1825"),
                     labels = c("1 years", "3 years", "5 years"),
                     values = c("#00468b", "#A73030FF", "#42b540")) +
  geom_roc(labels = F, stat = 'identity') + 
  style_roc() + 
  geom_abline(slope = 1, intercept = 0, color = 'gray', linetype=2) +
  theme_bw() +
  labs(title = "ROC for Prognosis Model") +
  theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
        panel.grid = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold")) +
  # annotate('text', x=.75, y=.25, label=paste('AUC of 1 years =', round(auc_y1,2))) + 
  # annotate('text', x=.75, y=.15, label=paste('AUC of 2 years =', round(auc_y3,2))) + 
  # annotate('text', x=.75, y=.05, label=paste('AUC of 3 years =', round(auc_y5,2))) +
  annotate("text", x=0.75, y=c(0.25, 0.15, 0.05),
           label = c(paste('AUC of 1 years =', format(auc_y1,nsmall=2)),
                     paste('AUC of 3 years =', format(auc_y3,nsmall=2)),
                     paste('AUC of 5 years =', format(auc_y5,nsmall=2))))
ROC
#08 风险评分与临床指标相关性分析-------------
setwd("/data/nas1/luchunlin/project/BJTC-282")
if (! dir.exists("./08_clinical_index")){
  dir.create("./08_clinical_index")
}
setwd("./08_clinical_index")

train_phenotype2<-phenotype
train_phenotype2$OS<-as.numeric(train_phenotype2$OS)
train_phenotype2$OS.time<-as.numeric(train_phenotype2$OS.time)

train_phenotype2$stage<-gsub('III','Stage3',train_phenotype2$stage)
train_phenotype2$stage<-gsub('II','Stage2',train_phenotype2$stage)
train_phenotype2$stage<-gsub('I','Stage1',train_phenotype2$stage)
train_phenotype2$stage<-gsub('Unknown',NA,train_phenotype2$stage)
train_phenotype2$stage<-gsub('unknown',NA,train_phenotype2$stage)
train_phenotype2$race<-gsub('not reported',NA,train_phenotype2$race)
train_phenotype2$race<-gsub('NA',NA,train_phenotype2$race)
train_phenotype2$race<-gsub('not allowed to collect',NA,train_phenotype2$race)
train_phenotype2$race<-gsub('other',NA,train_phenotype2$race)

colnames(train_phenotype2)<-c('id','OS','OS.time','gender','race','stage')
train_phenotype3 <- merge(sub_risk,
                          train_phenotype2,
                          by = "id")
write.table(train_phenotype3,
            file = "clinical_risk.csv",
            row.names = T,
            sep = "\t",
            quote = F)

library(ggpubr)
library(Ipaper)
library(ggthemes)

## 08-1 stage-----
my_comparisons <- list(c("Stage1","Stage2"),c("Stage2","Stage3"),c("Stage1","Stage3"))
stage_data <- data.frame(riskScore = train_phenotype3$riskScore,
                         stage = factor(train_phenotype3$stage,
                                        levels = c("Stage1","Stage2", "Stage3")))
stage_data <- na.omit(stage_data)

stage<-ggplot(stage_data,aes(x = stage, y = riskScore, fill = stage)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "riskScore",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("Stage") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              test = t.test,
              map_signif_level = T,
              y_position = c(4,5,6))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 20))+
  guides(fill='none')
stage
## 08-2 gender-----
my_comparisons <- list(c("male", "female"))
gender<-ggplot(train_phenotype3,aes(x = gender, y = riskScore, fill = gender)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "riskScore",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("Gender") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              test = t.test,
              map_signif_level = T,
              y_position = c(4,5,6))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 20))+
  guides(fill='none')
gender
## 08-2 race-----
my_comparisons <- list(c("white","asian"),c("asian","black or african american"),c("white","black or african american"))
stage_race <- data.frame(riskScore = train_phenotype3$riskScore,
                         race = factor(train_phenotype3$race,
                                        levels = c("white","asian", "black or african american")))
stage_race <- na.omit(stage_race)

race<-ggplot(stage_race,aes(x = race, y = riskScore, fill = race)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "riskScore",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("Race") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              test = t.test,
              map_signif_level = T,
              y_position = c(4,5,6))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 18),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13))+
  guides(fill='none')
race
library(patchwork)
all_clinical_index <- stage + gender + 
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


# 09 免疫检验点分子-----
setwd("/data/nas1/luchunlin/project/BJTC-282")
if (! dir.exists("./09_checkpoint")){
  dir.create("./09_checkpoint")
}
setwd("./09_checkpoint")

checkpoint <- read.table("checkpoint.txt",
                         header = F)
checkpoint <- checkpoint$V1
length(checkpoint)
fit <- lmFit(gsea_exp, design_risk)
fit2 <- contrasts.fit(fit ,compare_risk)
fit3 <- eBayes(fit2)
allGeneSets <- topTable(fit3, coef = 1, number = Inf)

checkpoint_DEG <- allGeneSets[which(rownames(allGeneSets)%in%checkpoint),]
dim(checkpoint_DEG)
# 46 6
logFC_cutoff <- 0
checkpoint_DEG$change = as.factor(
  ifelse(checkpoint_DEG$adj.P.Val < 0.05 & abs(checkpoint_DEG$logFC) > logFC_cutoff,
         ifelse(checkpoint_DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)

sig_checkpoint <- rownames(checkpoint_DEG[which(checkpoint_DEG$adj.P.Val < 0.05),])
length(sig_checkpoint)
#20
sig_checkpoint
#[1] "TNFSF9"   "LAIR1"    "TNFRSF14" "CD70"     "CD274"    "C10orf54" "CD40"     "CD200"    "TNFRSF8"  "BTNL2"   
#[15] "TIGIT"    "CD28"     "CTLA4"    "KIR3DL1"  "BTLA"     "TMIGD2"   "TNFRSF4"  "IDO2"     "ADORA2A"  "LAG3"   
write.table(sig_checkpoint, 
            'sig_checkpoint.xls',
            quote = F, sep="\t")
sig_expr <- as.data.frame(gsea_exp[sig_checkpoint,])
sig_expr <- log2(sig_expr + 1)
sig_expr$gene <- rownames(sig_expr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(Ipaper)
head(sig_expr[,1:3])
violin_dat <- gather(sig_expr, key=sample, value='log2(expr+1)', -c("gene"))
head(violin_dat)
violin_dat$group <- ifelse(violin_dat$sample %in% high_sample,
                           "High", "Low") 
head(violin_dat)

violin_plot <- ggplot(violin_dat, aes(x=gene, 
                                      y=`log2(expr+1)`,
                                      fill=group)) +
  #  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#E63F00", "#009FCC"), name = "Group")+
  labs(title="Immune Checkpoint", x="", y = "log2(expr+1)",size=20) +
  stat_compare_means(data = violin_dat,
                     mapping = aes(group = group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F) +
  #  geom_signif(comparisons = my_comparisons,
  #              test = t.test,
  #              map_signif_level = T)+
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
violin_plot

# 10 TIDE---------
setwd("/data/nas1/luchunlin/project/BJTC-282")
if (! dir.exists("./10_TIDE")){
  dir.create("./10_TIDE")
}
setwd("./10_TIDE")
tide_dat <- dat_final3
# fpkm转TPM
FPKM2TPM <- function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

# tide_dat[which(tide_dat<0)] <- 0
tide_dat <- apply(tide_dat,2,FPKM2TPM)
tide_dat <- log2(tide_dat + 1)
rownmean <- apply(tide_dat,1,mean)
tide_dat2 <- sweep(tide_dat, 1, rownmean)
dim(tide_dat2)
write.table(tide_dat2,
            file ="tide_dat.txt",
            sep = "\t",
            quote = F,
            row.names = T)
tide_result <- read.csv("MMRF_tide_result.csv",header = T)
View(tide_result)
tide_result2 <- subset(tide_result, select = c("Patient", "TIDE"))
tide_plot_dat <- data.frame(Patient=risk$id,
                            riskScore=risk$riskScore,
                            risk_group=risk$risk)
tide_plot_dat <- merge(tide_plot_dat, tide_result2, by = "Patient")
tide_plot_dat$risk_group <- factor(tide_plot_dat$risk_group,
                                   level = c(0, 1),
                                   labels = c("High Risk", "Low Risk"))

dim(tide_plot_dat)
library(ggplot2)
library(ggpubr)
install.packages('ggside')
# riskScore
{
  library(ggstatsplot)
  tide_cor <- ggscatterstats(data = tide_plot_dat,
                             y = riskScore,
                             x = TIDE,
                             centrality.para = "mean",
                             margins = "both",
                             xfill = "#A73030FF",
                             yfill = "#0073C2FF",
                             type = "pearson",
                             xlab = "TIDE Prediction Score",
                             marginal.type = "histogram",
                             title = "Relationship between TIDE Prediction Score and riskScore"
  )+theme(axis.title.x = element_text(size = 20),
          axis.title.y = element_text(size = 20),
          plot.title = element_text(size = 18))
  tide_cor
  ggsave(filename = "tide_riskscore_cor.png", height = 10, width = 12,tide_cor)
  ggsave(filename = "tide_riskscore_cor.pdf", height = 10, width = 12,tide_cor)
  
  tide_box <- ggboxplot(tide_plot_dat,
                        x = "risk_group",
                        y = "TIDE",
                        fill = "risk_group",
                        palette =c("#A73030FF", "#0073C2FF")) +
    stat_compare_means(label.y = 3.2) +
    theme_bw() +
    labs(title = "", x = "", y = "TIDE Prediction Score") +
    theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18),
          axis.text.x=element_text(colour="black",face="bold",size=15), 
          axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12),
          axis.title.y=element_text(size=18,face="bold"),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  tide_box
  ggsave(filename = "tide_riskscore_boxplot.png", height = 5, width = 4,tide_box)
  ggsave(filename = "tide_riskscore_boxplot.pdf", height = 5, width = 4,tide_box)
}
# 11 IPS---------
setwd("/data/nas1/luchunlin/project/BJTC-282")
if (! dir.exists("./11_IPS")){
  dir.create("./11_IPS")
}
setwd("./11_IPS")
ips_dat<-tide_dat

## IOBR包
if (!requireNamespace("IOBR", quietly = TRUE))
  devtools::install_github("IOBR/IOBR")
library(IOBR)

ips<-deconvo_tme(eset = ips_dat, method = "ips", plot= FALSE)
# ips<-IPS_calculation(eset = ips_dat,plot = F)
head(ips)
ips_result<-subset(ips,select=c('ID','IPS_IPS'))
colnames(ips_result)<-c('ID','IPS')
ips_plot_dat<-data.frame(ID=risk$id,
                         riskScore=risk$riskScore,
                         risk_group=risk$risk)
ips_plot_dat<-merge(ips_plot_dat,ips_result,by='ID')
ips_plot_dat$risk_group<-factor(ips_plot_dat$risk_group,
                                levels = c(0,1),
                                labels = c("High Risk", "Low Risk"))

library(ggridges)
ggplot(ips_plot_dat, aes(x = IPS, y = risk_group)) +
  geom_density_ridges_gradient(aes(fill = risk_group),
                               scale = 3, size = 0.3) +
  theme(legend.position = "none")+
  labs(y = "")+
  theme(axis.text.y=element_text(hjust=0.5,colour="black",size=12),
        axis.text.x=element_text(hjust=0.5,colour="black",size=12))
# 12 化疗药物敏感性---------
setwd('F:/luchunlin/project/BJTC-282/')
if (! dir.exists("./12_Medicinal_Sensity")){
  dir.create("./12_Medicinal_Sensity")
}
setwd("./12_Medicinal_Sensity")

#install.packages("pRRophetic_0.5.tar.gz", repos = NULL, dependencies = TRUE)
library(pRRophetic)
library(ggplot2)
set.seed(12345)
model_expr<-dat_final3[lasso_geneids, risk$id]
riskscore<-data.frame(risk$id,risk$risk)
colnames(riskscore)<-c('sample','risk')
riskscore$risk[which(riskscore$risk==1)] <-'Low risk'
riskscore$risk[which(riskscore$risk==0)] <-'High risk'
head(riskscore)
drug<-read.table(file = 'drugs.txt',sep='\t',header=F)
ic50<-data.frame(riskscore$sample)

a<-data.frame(row.names=riskscore$sample,riskscore$risk)

colnames(a)<-'risk'

cnt<-1

while (cnt < 139) {
  
  predictedPtype <- pRRopheticPredict(as.matrix(model_expr), drug[cnt,],selection=1)
  
  Tipifarnib<-data.frame(predictedPtype)
  
  colnames(Tipifarnib)<-drug[cnt,]
  
  a<-cbind(a,Tipifarnib)
  
  cnt = cnt + 1
}

write.table(a,'IC50.xls',sep='\t',quote=F)

b<-a
b[b<0]<-NA
# 先写成函数的形式，方便调用
removeRowsAllNa  <- function(x){x[apply(x, 1, function(y) any(!is.na(y))),]}
removeColsAllNa  <- function(x){x[, apply(x, 2, function(y) any(!is.na(y)))]}
c<-removeColsAllNa(b)
na_flag <- apply(is.na(c), 2, sum)
x <- c[, which(na_flag == 0)]
View(x)
dim(x)
# [1] 858  98
medicinal_result <- t(subset(x, select = -risk)) 
high_group <- risk$id[which(risk$risk==0)]
low_group <- risk$id[which(risk$risk==1)]
pvalue = padj = log2FoldChange <- matrix(0, nrow(medicinal_result), 1)
for (i in 1:nrow(medicinal_result)){
  pvalue[i, 1] = p.value = wilcox.test(medicinal_result[i, high_group],
                                       medicinal_result[i, low_group])$p.value
  log2FoldChange[i, 1] = mean(medicinal_result[i, high_group]) - 
    mean(medicinal_result[i, low_group])
}
padj <- p.adjust(as.vector(pvalue), "fdr", n = length(pvalue))
rTable <- data.frame(log2FoldChange, 
                     pvalue, 
                     padj,
                     row.names = rownames(medicinal_result))
high_group_res <- signif(apply(medicinal_result[rownames(rTable), high_group], 
                               1,
                               median), 4)
low_group_res <- signif(apply(medicinal_result[rownames(rTable), low_group], 
                              1, 
                              median), 4)
rTable <- data.frame(high_group_res, 
                     low_group_res,
                     rTable[, c("padj", "pvalue", "log2FoldChange")])
rTable$drugs <- rownames(rTable)
rTable$sig <- ifelse(rTable$padj < 0.05,
                     ifelse(rTable$padj < 0.01, 
                            ifelse(rTable$padj < 0.001,
                                   ifelse(rTable$padj < 0.0001,
                                          paste(rTable$drugs, "****",  sep = ""),
                                          paste(rTable$drugs, "***", sep = "")),
                                   paste(rTable$drugs, "**", sep = "")),
                            paste(rTable$drugs, "*",  sep = "")), 
                     rTable$drugs)

write.table(rTable,
            file = "drugs_wilcox_test.xls",
            quote = F,
            row.names = F,
            sep = '\t')

### 发散条形图绘制

#install.packages('ggprism')
library(ggprism)
## 横坐标药物，纵坐标：IC50(H)/IC50(L)-1
dat_plot<-data.frame(drug=rownames(rTable),
                     'IC50(H)/IC50(L)-1'=(rTable$high_group_res/rTable$low_group_res-1),
                     pvalue=rTable$pvalue,
                     padj=rTable$padj)
dat_plot$threshold=factor(ifelse(dat_plot$padj<0.05&dat_plot$pvalue<0.05,'P<0.05 & FDR<0.05',ifelse(dat_plot$pvalue<0.05&dat_plot$padj>0.05,'P<0.05 & FDR>0.05','P>0.05 & FDR>0.05')))
dat_plot<-dat_plot%>%arrange(desc(dat_plot$IC50.H..IC50.L..1))

dat_plot$drug<-factor(dat_plot$drug,levels=dat_plot$drug)

p <- ggplot(data = dat_plot,aes(x = drug,y = IC50.H..IC50.L..1,fill = threshold)) +
  geom_col()+
  scale_fill_manual(values = c('P<0.05 & FDR<0.05'= '#5F9EA0','P>0.05 & FDR>0.05'='#cccccc','P<0.05 & FDR>0.05'='#FFD700')) +
  xlab('') + 
  ylab('Exp(Median IC50(H))/Exp(Median IC50(L))-1') + 
  theme_prism(border = T) +
  theme(
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 65,size = 7,
                               hjust = 1,vjust = 1),
    axis.text.y = element_text(size = 13),
    legend.position = c(0.85,0.85),
    legend.text = element_text(size = 8,face = 'bold')
  )
p







