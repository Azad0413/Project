##rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/BJTC-276")
if (! dir.exists("./00_rawdata")){
  dir.create("./00_rawdata")
}
setwd("./00_rawdata")
##01-1 TCGA OV--------
library(TCGAbiolinks)
library(readr)
library(readxl)
library(tidyverse)
## 读取从xena下载的数据
expr<-read_tsv(file = 'TCGA-OV.htseq_counts.tsv')
expr<-as.data.frame(expr)
rownames(expr)<-expr[,1]
expr<-expr[,-1]
## xena下载的数据经过了log2+1转化，需要将其还原
expr<-2^expr-1
## 对数据进行id转化
genecode<-read.table(file = 'gencode.v22.annotation.gene.probeMap')
probe2symbol<-genecode[,(1:2)]
colnames(probe2symbol)<-c('ID','symbol')
probe2symbol<-probe2symbol[-1,]
dat<-expr
dat$ID <- rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symbol$ID<-as.character(probe2symbol$ID)
dat<-dat %>%
  inner_join(probe2symbol,by='ID')%>% 
  select(-ID)%>%     ## 去除多余信息
  select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('TCGA',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
dim(dat)
dat_tcga<-dat
## 加载注释文件
#library("rtracklayer")
#gtf_data = import('gencode.v40.annotation.gtf') #gtf的路径
#gtf_data = as.data.frame(gtf_data)
#提取表达谱内的lncRNA,记得先把表达谱内的Ensembl_ID改成"gene_id")
#protein_coding=gtf_data%>%
#  dplyr::filter(type=="gene",gene_type=="protein_coding")%>%
#  dplyr::select(gene_id,gene_type,gene_name)


#keep<-rowSums(dat>0)>=floor(0.75*ncol(dat))
#dat<-dat[keep,]
#dat_tcga<-dat[rownames(dat)%in%protein_coding$gene_name,]
## 15866 

##fpkm
expr_fpkm<-read_tsv(file = 'TCGA-OV.htseq_fpkm.tsv')
expr_fpkm<-as.data.frame(expr_fpkm)
rownames(expr_fpkm)<-expr_fpkm[,1]
expr_fpkm<-expr_fpkm[,-1]
## xena下载的数据经过了log2+1转化，需要将其还原
expr_fpkm<-2^expr_fpkm-1
## 对数据进行id转化
genecode<-read.table(file = 'gencode.v22.annotation.gene.probeMap')
probe2symbol<-genecode[,(1:2)]
colnames(probe2symbol)<-c('ID','symbol')
probe2symbol<-probe2symbol[-1,]
dat_fpkm<-expr_fpkm
dat_fpkm$ID <- rownames(dat_fpkm)
dat_fpkm$ID<-as.character(dat_fpkm$ID)
probe2symbol$ID<-as.character(probe2symbol$ID)
dat_fpkm<-dat_fpkm %>%
  inner_join(probe2symbol,by='ID')%>% 
  select(-ID)%>%     ## 去除多余信息
  select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('TCGA',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
dim(dat_fpkm)
# fpkm转TPM
FPKM2TPM <- function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

dat_tpm <- apply(dat_fpkm,2,FPKM2TPM)

#dat_fpkm<-dat_fpkm[rownames(dat_fpkm)%in%rownames(dat_tcga),]
## 01-2 差异基因集------
library(GEOquery)
library(Biobase)
##  GSE54388, GSE69428, GSE14407, GSE10971联合作为训练集
### 01-2-1GSE23554；GSE14764------
gset1<-getGEO("GSE54388",
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
pd1<-pData(a1)
pd1<-data.frame(sample=pd1$geo_accession,
                statu=pd1$`diagnosis:ch1`)
pd1$group<-ifelse(pd1$statu=='healthy (normal)','control','OV')
## 16个癌症，6个对照
### 01-2-2 GSE69428--------
gset2<-getGEO("GSE69428",
              destdir = '.',
              GSEMatrix = T,
              getGPL = F)
expr2<-as.data.frame(exprs(gset2[[1]]))
a2=gset2[[1]]
dat2<-expr2
dat2$ID<-rownames(dat2)
dat2$ID<-as.character(dat2$ID)
dat2<-dat2 %>%
  inner_join(probe2symobl1,by='ID')%>% 
  select(-ID)%>%     ## 去除多余信息
  select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
pd2<-pData(a2)
pd2<-pd2[c(1:20),]
dat2<-dat2[,pd2$geo_accession]
pd2<-data.frame(sample=pd2$geo_accession,
                statu=pd2$title)
pd2$group<-c(rep('OV',10),rep('control',10))
## 10个肿瘤10个对照 提取处理

### 把2个数据集整合在一起
dat.all<-merge(dat1,dat2,by=0)%>%
  column_to_rownames(var = 'Row.names')
#dat.all<-merge(dat.all,dat3,by=0)%>%
#  column_to_rownames(var = 'Row.names')
## 分组
group<-rbind(pd1,pd2)
table(group$group)
## control      OV 
##      16     26 

### 01-2-5 去除批次效应-----
## 批次信息 一列为样本名字与表达矩阵一致，另一列为批次信息，用1,2表示
library(sva)
library(bladderbatch)
## 转化为矩阵
sif<-data.frame(batch=c(rep('1',22),rep('2',20)))
dat.all<-dat.all[,group$sample]
dat.all<-t(dat.all)
dat_pca<-cbind(sif,dat.all)
dat_pca<-as.data.frame(lapply(dat_pca,as.numeric))
dat_pca<-t(dat_pca)
modcombat=model.matrix(~1,data=sif)
combat=ComBat(dat = dat_pca,
              batch = sif$batch,
              mod = modcombat,
              par.prior = T,
              prior.plots = F)
colnames(combat)<-rownames(dat.all)
dat.final<-combat[-1,]
group<-group[order(group$group),]
dat.final<-dat.final[,group$sample]
write.table(dat.final,'dat_final.xls',
            sep = '\t',
            quote = F)


# 02 差异基因鉴定----------
setwd("/data/nas1/luchunlin/project/BJTC-276")
if (! dir.exists("./01_DEGs")){
  dir.create("./01_DEGs")
}
setwd("./01_DEGs")
library(limma)
group<-group[,c(1,3)]

write.table(group,file = 'group.xls',
            sep = '\t',
            quote = F,
            row.names = T)
type<-group[,2]
design <- model.matrix(~ -1+factor(type,levels=c('control','OV'))) 
colnames(design)<-c('control','OV')
rownames(design)<-group$id
# 计算logFC和FDR
library(limma)
# 对每一个基因进行线性模型构建
fit=lmFit(dat.final,design)
# 构建比较矩阵
contrast.matrix=makeContrasts(ControlVSMG=OV-control,levels = design)
# 构建芯片数据线性模型，计算估计的相关系数和标准差
fit2=contrasts.fit(fit,contrast.matrix)
# 基于贝叶斯计算T值，F值和log-odds
fit2<-eBayes(fit2,0.01)
tempOutput = topTable(fit2, coef=1, n=Inf)
DEG= na.omit(tempOutput)
write.table(DEG,file = "DEG.xls",
            quote = F,
            sep = "\t",
            row.names = T)
# 筛选差异基因
logFC_cutoff <- 0

DEG$change = as.factor(
  ifelse(DEG$adj.P.Val < 0.05 & abs(DEG$logFC) > logFC_cutoff,
         ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
sig_diff <- subset(DEG,
                   DEG$adj.P.Val < 0.05 & abs(DEG$logFC) > logFC_cutoff)
dim(DEG)
dim(sig_diff)
## 3872    7
summary(sig_diff$change)

# DOWN  NOT   UP 
# 1438    0 2434 
# 输出结果
write.table(DEG,file = "DEG_all.xls",
            quote = F,
            sep = "\t",
            row.names = T)
write.table(sig_diff,file = "DEG_sig.xls",
            quote = F,
            sep = "\t",
            row.names = T)
### 火山图---------
# devtools::install_github("kongdd/Ipaper")
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(Ipaper)
library(scales)
library(ggrepel)

dat_rep<-DEG[rownames(DEG)%in%
               c(head(rownames(subset(sig_diff,sig_diff$logFC>3.17)),10),
                 head(rownames(subset(sig_diff,sig_diff$logFC< -3.24)),10)),]
colnames(dat_rep)
volcano_plot<- ggplot(data = DEG, 
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
    data = dat_rep,
    aes(label = rownames(dat_rep)),
    max.overlaps = 20,
    size = 4,
    box.padding = unit(0.5, "lines"),
    min.segment.length = 0,
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )+
  labs(x = "log (Fold Change)",
       y = "-log10 (adj.P.Val)")
volcano_plot
ggsave('volcano.png', volcano_plot,width = 8, height = 7)
ggsave('volcano.pdf', volcano_plot,width = 8, height = 7)
### 热图--------
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(gplots)
group_rt<-group
rt<-dat.final[,group_rt$sample]
group_rt<-group_rt$group%>%as.data.frame()
colnames(group_rt)<-'group'
rownames(group_rt)<-colnames(rt)
heat<-rt[rownames(rt)%in%rownames(dat_rep),]
#x<-log2(heat+1)
x<-t(scale(t(heat)))
ann_colors<-list(
  group = c(control="lightblue",OV="darkorange"))
pheatmap(mat=x,
         annotation_col = group_rt,
         color=bluered(100),
         scale = "row",
         annotation_colors = ann_colors,
         fontsize = 10,
         show_colnames = F,
         cluster_cols = F,
         cluster_rows = T)

# 03 差异免疫细胞鉴定-----------
setwd("/data/nas1/luchunlin/project/BJTC-276")
if (! dir.exists("./02_ssGSEA")){
  dir.create("./02_ssGSEA")
}
setwd("./02_ssGSEA")
## 03-1 ssGSEA--------
library(GSVA)
gene_set <- read.table("/data/nas1/luchunlin/pipeline/ssGSEA/mmc3.txt",
                       header = T,
                       sep ="\t")
dat.final2 <- dat.final[,rownames(group_rt)]
gene_list <- split(as.matrix(gene_set)[,1],
                   gene_set[,2])

ssgsea_score = gsva(dat.final2, gene_list, 
                    method = "ssgsea", 
                    ssgsea.norm = TRUE, 
                    verbose = TRUE)
write.table(ssgsea_score,
            file = "ssgsea_result.xls",
            sep = "\t",
            quote = F)
## 富集分数画热图
annotation_col<-as.data.frame(group_rt$group)
colnames(annotation_col)='Group'
rownames(annotation_col)=colnames(ssgsea_score)
color.key<-c("#3300CC", "#3399FF", "white", "#FF3333", "#CC0000")
ann_colors<-list(
  Group = c(control="#00FFFF",OV="#FFAEB9"))
pheatmap(
  ssgsea_score,
  color = colorRampPalette(color.key)(50),
  border_color = 'darkgrey',
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  labels_row = NULL,
  clustering_method = 'ward.D2',
  show_rownames = T,
  show_colnames = F,
  fontsize_col = 5,
  cluster_cols = F,
  cluster_rows = T)
## 03-2 差异免疫细胞鉴定-------
tiics_result <- ssgsea_score
pvalue = padj = log2FoldChange <- matrix(0, nrow(tiics_result), 1)
group_OV<-group[group$group=='OV',]
group_OV<-as.character(group_OV$sample)
group_control<-group[group$group=='control',]
group_control<-as.character(group_control$sample)

for (i in 1:nrow(tiics_result)){
  pvalue[i, 1] = p.value = wilcox.test(tiics_result[i, group_OV],
                                       tiics_result[i, group_control])$p.value
  log2FoldChange[i, 1] = mean(tiics_result[i, group_OV]) - 
    mean(tiics_result[i, group_control])
}
padj <- p.adjust(as.vector(pvalue), "fdr", n = length(pvalue))
rTable <- data.frame(log2FoldChange, 
                     pvalue, 
                     padj,
                     row.names = rownames(tiics_result))
control <- signif(apply(tiics_result[rownames(rTable), group_control], 
                        1,
                        mean), 4)
OV <- signif(apply(tiics_result[rownames(rTable), group_OV], 
                       1, 
                       mean), 4)
rTable <- data.frame(control, 
                     OV,
                     rTable[, c("padj", "pvalue", "log2FoldChange")])
rTable$immune_cell <- rownames(rTable)
rTable$sig <- ifelse(rTable$padj < 0.05,
                     ifelse(rTable$padj < 0.01, 
                            ifelse(rTable$padj < 0.001,
                                   ifelse(rTable$padj < 0.0001,
                                          paste(rTable$immune_cell, "****",  sep = ""),
                                          paste(rTable$immune_cell, "***", sep = "")),
                                   paste(rTable$immune_cell, "**", sep = "")),
                            paste(rTable$immune_cell, "*",  sep = "")), 
                     rTable$immune_cell)

diff_Table<-rTable[which(rTable$padj<0.05),]
## 15
write.table(rTable,
            file = "tiics_wilcox_test.xls",
            quote = F,
            row.names = F,
            sep = '\t')
write.table(diff_Table,
            file = "diff_tiics_wilcox_test.xls",
            quote = F,
            row.names = F,
            sep = '\t')
### 箱线图----
library(tidyr)
library(ggplot2)
library(ggpubr)
library(Ipaper)

xCell2 <- data.frame(Immune_Cell=rownames(tiics_result), 
                     tiics_result, 
                     pvalue=rTable$padj)
xCell3 <- xCell2[which(xCell2$pvalue<0.05),]
#xCell3<-xCell2
diff_tiics <- rownames(xCell3)
violin_dat <- gather(xCell3, key=Group, value=score, -c("Immune_Cell","pvalue"))

violin_dat$Group <- ifelse(gsub("\\.","-",violin_dat$Group) %in% group_control,
                           "control", "OV") 
head(violin_dat)
boxplot_diff_TIICs <- ggplot(violin_dat, aes(x=Immune_Cell, 
                                             y=score,
                                             fill=Group)) +
  # geom_violin(trim=T,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)#"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar",
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  geom_point(aes(fill = Group),
             size = 0.05,
             position = position_dodge(0.9))+
  scale_fill_manual(values= c("#48D1CC", "#FF7256"))+ #设置填充的颜色
  labs(title="", x="", y = "Score",size=20) +
  stat_compare_means(data = violin_dat,
                     mapping = aes(group = Group),
                     label ="p.signif",
                     hide.ns = F) +
  theme_bw()+#把背景设置为白底
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18), # 将图表标题居中
        axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=10), #设置x轴刻度标签的字体显示倾斜角度为45度，并向下调整1(hjust = 1)，字体大小为14
        axis.text.y=element_text(hjust=0.5,colour="black",size=12), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.x=element_text(size=16,face="bold"),#设置x轴标题的字体属性
        axis.title.y=element_text(size=14,face="bold"), #设置y轴标题的字体属性
        legend.text=element_text(face="bold", hjust = 0.5,colour="black", size=11), #设置图例的子标题的字体属性
        legend.title=element_text(face="bold", colour="black", size=11),#设置图例的总标题的字体属性
        #legend.justification=c(-0.1,1.2), #可调整图例的位置。##(1,1)第一个1是调整图例在图的内外(左右移动)，第二个1是在图中上下移动。
        #legend.position=c(0, 1.04), #legend.position=c(0,1)左上角，(1,1)是在右上角。
        panel.grid.major = element_blank(), #不显示网格线
        panel.grid.minor = element_blank()) #不显示网格线
boxplot_diff_TIICs
# 04 WGCNA----------
setwd("/data/nas1/luchunlin/project/BJTC-276")
if (! dir.exists("./03_WGCNA")){
  dir.create("./03_WGCNA")
}
setwd("./03_WGCNA")
## 基于OV样本 提取表达矩阵
dat.ov<-dat.final[,group_OV]
library(WGCNA)
library(reshape2)
library(stringr)
options(stringsAsFactors = F)
enableWGCNAThreads()
exprMat<-dat.ov
dim(exprMat)
type = "unsigned"
corType = "pearson"
corFnc = ifelse(corType=="pearson", cor, bicor)
# 对二元变量，如样本性状信息计算相关性时，
# 或基因表达严重依赖于疾病状态时，需设置下面参数
maxPOutliers = ifelse(corType=="pearson",1,0.05)
# 关联样品性状的二元变量时，设置
robustY = ifelse(corType=="pearson",T,F)
dataExpr <- exprMat[rownames(dat.ov),]
## 04-1 数据筛选-----
## 筛选中位绝对偏差（MAD）前75%的基因，至少MAD大于0.01
## 筛选后会降低运算量，也会失去部分信息
## 也可不做筛选，使MAD大于0即可
m.mad <- apply(dataExpr,1,mad)
dataExprVar <- dataExpr[which(m.mad >
                                max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
## 转换为样品在行，基因在列的矩阵
dataExpr <- as.data.frame(t(dataExprVar))
## 检测缺失值
gsg = goodSamplesGenes(dataExpr, verbose = 3)
gsg$allOK
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:",
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:",
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
dim(dataExpr)
# [1] 26 13555

## 04-2 软阈值筛选----
## 样本聚类，查看是否有利群样本
hclus
sampleTree = hclust(dist(dataExpr,method = 'euclidean'), method = "complete")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")

#abline(h = 130, col = "red")
# Determine cluster under the line
# 剪枝算法，cutHeight修剪树枝的高度，minSize集群最小数
#clust = cutreeStatic(sampleTree, cutHeight = 135, minSize = 10)
#table(clust)
## 14 376
# clust 1 contains the samples we want to keep.
#keepSamples = (clust==1)
## 符合要求的数据
#dataExpr = dataExpr[keepSamples,]
## 提取列
nGenes = ncol(dataExpr)
## 提取行
nSamples = nrow(dataExpr)
# 设置网络构建参数选择范围，计算无尺度分布拓扑矩阵
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(dataExpr, powerVector=powers,
                        networkType=type, verbose=5)
par(mfrow = c(1,2))
cex1 = 0.9
# 横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，
# 网络越符合无标度特征 (non-scale)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
# 筛选标准。R-square=0.85
abline(h=0.85,col="red")
# Soft threshold与平均连通性
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,
     cex=cex1, col="red")
power = sft$powerEstimate
power
## 6
## 04-3一步法网络构建---------
## One-step network construction and module detection##
# power: 上一步计算的软阈值
net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                       TOMType = type, minModuleSize = 50,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType,
                       maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                       saveTOMFileBase = "DiffGene_TOM",
                       verbose = 3)
# 根据模块中基因数目的多少，降序排列，依次编号为 `1-最大模块数`。
# **0 (grey)**表示**未**分入任何模块的基因。
table(net$colors)  
#    0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20   
##  847 3084 2439 2328 1019  574  504  407  347  326  256  215  192  178  150  144  140  124  113   92   76 
## 04-4 层级聚类数展示各个模块
## 灰色的为**未分类**到模块的基因。
# Convert labels to colors for plotting
moduleLabels = net$colors
# moduleColors = labels2colors(moduleLabels)
table(moduleLabels)
table(moduleColors)
# Plot the dendrogram and the module colors underneath
# 如果对结果不满意，还可以recutBlockwiseTrees，节省计算时间
# plotdendroandcolors 函数，接受一个聚类的对象，以及该对象里面包含的所有个体所对应的颜色。
png(filename = "cluster_dendrogram.png", height = 600, width = 800)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
pdf(file = "cluster_dendrogram.pdf", height = 7, width = 10)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
## 04-5 绘制模块间的相关性热图-------
# module eigengene, 可以绘制线图，作为每个模块的基因表达趋势的展示
MEs = net$MEs
### 不需要重新计算，改下列名字就好
### 官方教程是重新计算的，起始可以不用这么麻烦
MEs_col = MEs
colnames(MEs_col) = paste0("ME", labels2colors(
  as.numeric(str_replace_all(colnames(MEs),"ME",""))))
MEs_col = orderMEs(MEs_col)
# 根据基因间表达量进行聚类所得到的各模块间的相关性图
# marDendro/marHeatmap 设置下、左、上、右的边距
sizeGrWindow(6,6)
plotEigengeneNetworks(MEs_col, 
                      setLabels = "Eigengene dendrogram and heatmap", 
                      marDendro = c(1,3,4,4),
                      marHeatmap = c(3,3,0,2), 
                      plotDendrograms = T, 
                      xLabelsAngle = 90)

# Plot the dendrogram
plotEigengeneNetworks(MEs_col, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)


## 04-6 关联表型数据------
trait<-t(tiics_result[rownames(diff_Table),])
trait<-trait[rownames(MEs_col),]


## 模块与表型数据关联
if (corType=="pearson") {
  modTraitCor = cor(MEs_col, trait, use = "p")
  modTraitP = corPvalueStudent(modTraitCor, nSamples)
} else {
  modTraitCorP = bicorAndPvalue(MEs_col, trait, robustY=robustY)
  modTraitCor = modTraitCorP$bicor
  modTraitP   = modTraitCorP$p
}
modTraitCor2 <- modTraitCor[!rownames(modTraitCor) %in% "MEgrey",]
modTraitP2 <- modTraitP[!rownames(modTraitP) %in% "MEgrey",]
# signif表示保留几位小数
textMatrix = paste(signif(modTraitCor2, 2), 
                   "\n(", signif(modTraitP2, 1), 
                   ")", 
                   sep = "")
dim(textMatrix) = dim(modTraitCor2)
par(mar = c(10, 10, 5, 2));
labeledHeatmap(Matrix = modTraitCor2, 
               xLabels = colnames(trait), 
               yLabels = colnames(MEs_col)[-16], 
               cex.lab = 1, 
               ySymbols = colnames(MEs_col)[-16], 
               colorLabels = FALSE, 
               colors = blueWhiteRed(50, 0.8),
               textMatrix = textMatrix, 
               setStdMargins = FALSE, 
               cex.text = 0.8, 
               zlim = c(-1,1),
               main = paste("Module-TIICs relationships"))

## 04-7 挑选免疫浸润细胞相关模块基因
{
## pink
module='pink'
probes=colnames(dataExpr)
inModule=(moduleColors==module)
modProbes=probes[inModule]
modGenes<-as.data.frame(modProbes)
colnames(modGenes)<-'modgene'
## 574
## 加1个cyan
module2='cyan'
probes=colnames(dataExpr)
inModule2=(moduleColors==module2)
modProbes2=probes[inModule2]
modGenes2<-as.data.frame(modProbes2)
colnames(modGenes2)<-'modgene'
modGenes<-rbind(modGenes,modGenes2)
## 1078
## 加1个salmon
module3='salmon'
probes=colnames(dataExpr)
inModule3=(moduleColors==module3)
modProbes3=probes[inModule3]
modGenes3<-as.data.frame(modProbes3)
colnames(modGenes3)<-'modgene'
modGenes<-rbind(modGenes,modGenes3)
## 1404
## 加1个greenyellow
module4='greenyellow'
probes=colnames(dataExpr)
inModule4=(moduleColors==module4)
modProbes4=probes[inModule4]
modGenes4<-as.data.frame(modProbes4)
colnames(modGenes4)<-'modgene'
modGenes<-rbind(modGenes,modGenes4)
## 1528
## 加1个black
module5='black'
probes=colnames(dataExpr)
inModule5=(moduleColors==module5)
modProbes5=probes[inModule5]
modGenes5<-as.data.frame(modProbes5)
colnames(modGenes5)<-'modgene'
modGenes<-rbind(modGenes,modGenes5)
## 1641
## 加1个tan
module6='tan'
probes=colnames(dataExpr)
inModule6=(moduleColors==module6)
modProbes6=probes[inModule6]
modGenes6<-as.data.frame(modProbes6)
colnames(modGenes6)<-'modgene'
modGenes<-rbind(modGenes,modGenes6)
## 1785
}
{## 加1个pink
module7='pink'
probes=colnames(dataExpr)
inModule7=(moduleColors==module7)
modProbes7=probes[inModule7]
modGenes7<-as.data.frame(modProbes7)
colnames(modGenes7)<-'modgene'
modGenes<-rbind(modGenes,modGenes7)
## 2132
## 加1个lightcyan
module8='lightcyan'
probes=colnames(dataExpr)
inModule8=(moduleColors==module8)
modProbes8=probes[inModule8]
modGenes8<-as.data.frame(modProbes8)
colnames(modGenes8)<-'modgene'
modGenes<-rbind(modGenes,modGenes8)
## 2272
## 加1个purple
module9='purple'
probes=colnames(dataExpr)
inModule9=(moduleColors==module9)
modProbes9=probes[inModule9]
modGenes9<-as.data.frame(modProbes9)
colnames(modGenes9)<-'modgene'
modGenes<-rbind(modGenes,modGenes9)
## 2582
## 加1个brown
module10='brown'
probes=colnames(dataExpr)
inModule10=(moduleColors==module10)
modProbes10=probes[inModule10]
modGenes10<-as.data.frame(modProbes10)
colnames(modGenes10)<-'modgene'
modGenes<-rbind(modGenes,modGenes10)
## 4856

## 加1个greenyellow
module11='greenyellow'
probes=colnames(dataExpr)
inModule11=(moduleColors==module11)
modProbes11=probes[inModule11]
modGenes11<-as.data.frame(modProbes11)
colnames(modGenes11)<-'modgene'
modGenes<-rbind(modGenes,modGenes11)
## 5071
}
write.table(modGenes,file = 'modGenes.xls',sep = '\t',quote = F,row.names = F)
# 05 EMT-免疫相关差异基因----------
setwd("/data/nas1/luchunlin/project/BJTC-276")
if (! dir.exists("./04_EMT_IMU")){
  dir.create("./04_EMT_IMU")
}
setwd("./04_EMT_IMU")
## 05-1 差异免疫细胞相关基因------
## 差异表达基因与WGCNA基因取交集
DE_mod<-sig_diff[rownames(sig_diff)%in%modGenes$modgene,]
DE_modgene<-data.frame(symbol=rownames(DE_mod))
## 1313
write.table(DE_modgene,file = 'DE_mod.xls',sep = '\t',quote = F,row.names = F)
library(ggvenn)
mydata1<-list('DEGs'=rownames(sig_diff),'MOD gene'=modGenes$modgene)
ggvenn(mydata1,c('DEGs','MOD gene'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')


## 05-2 差异免疫-EMT基因------
## 上面得到的基因与EMT基因取交集
EMT_gene1<-read_xlsx('EMT_geneset.xlsx')
colnames(EMT_gene1)<-'symbol'
EMT_gene2<-read_xlsx('EMT.xlsx')
colnames(EMT_gene2)<-'symbol'
#EMT_gene2<-EMT_gene2[-752,]
EMT_gene<-rbind(EMT_gene1,EMT_gene2)
#EMT_gene<-EMT_gene2
EMT_gene<-EMT_gene[!duplicated(EMT_gene$symbol),]
##1317
write.table(EMT_gene,file = 'EMT_gene.xls',sep = '\t',quote = F,row.names = F)
#EMT_gene<-EMT_gene[-c(476),]
##DE_emtmod<-sig_diff[rownames(sig_diff)%in%EMT_gene$HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION,]
DE_emtmod<-DE_mod[rownames(DE_mod)%in%EMT_gene$symbol,]
## 139
DE_emtmodgene<-data.frame(symbol=rownames(DE_emtmod))
write.table(DE_emtmodgene,file = 'DE_emtmod.xls',sep = '\t',quote = F,row.names = F)

mydata2<-list('Mod-DEGs'=rownames(DE_mod),'EMT'=EMT_gene$symbol)
ggvenn(mydata2,c('Mod-DEGs','EMT'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')

# 06 GO/KEGG----------
setwd("/data/nas1/luchunlin/project/BJTC-276")
if (! dir.exists("./05_GO_KEGG")){
  dir.create("./05_GO_KEGG")
}
setwd("./05_GO_KEGG")
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(GOplot)
gene_transform <- bitr(rownames(DE_emtmod),
                       fromType = "SYMBOL",
                       toType = c("ENTREZID", "ENSEMBL", "REFSEQ"),
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
go_dot<-dotplot(ego, showCategory=5, split="ONTOLOGY") +
  facet_grid(ONTOLOGY ~ ., scales = "free")
go_dot
CC <- enrichGO(gene = gene_transform$ENTREZID,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "CC",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
cc_dot<-dotplot(CC, showCategory=10)+
  ggtitle(label = 'Cellular Component')
cc_dot
BP <- enrichGO(gene = gene_transform$ENTREZID,
               OrgDb = org.Hs.eg.db,
               keyType = "ENTREZID",
               ont = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               readable = TRUE)
bp_dot<-dotplot(BP, showCategory=10)+
  ggtitle(label = 'Biological Process')
bp_dot
MF <- enrichGO(gene = gene_transform$ENTREZID,
               OrgDb = org.Hs.eg.db,
               keyType = "ENTREZID",
               ont = "MF",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               readable = TRUE)
mf_dot<-dotplot(MF, showCategory=10)+
  ggtitle(label = 'Molecular Function')
mf_dot

## KEGG富集分析（气泡图）
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "hsa",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk@result,file = "KEGG.xls",sep = "\t",quote = F,row.names = F)
kk2<-data.frame(Category = "ALL",ID = kk$ID,Term = kk$Description, Genes = gsub("/", ", ", kk$geneID), adj_pval = kk$p.adjust)
## 读取logFC文件
genelist<-data.frame(ID=rownames(DE_emtmod),
                     logFC=DE_emtmod$logFC)
rownames(genelist)<-genelist$ID
genelist<-genelist[order(genelist$logFC,decreasing = T),]
circ<-circle_dat(kk2,genelist)
circ<-circ[c(1:83),]
# GOBubble(circ, labels = 3,table.legend =F)
# GOCircle(circ,rad1=2.5,rad2=3.5,label.size=4,nsub=10) 
#termNum<-20
#geneNum<-nrow(genelist)
chord<-chord_dat(circ,genelist)
kegg_chord<-GOChord(chord,
                    gene.order = 'logFC',
                    gene.space = 0.25,
                    gene.size = 8,
                    space = 0.01,
                    lfc.col = c('red','white','blue'),
                    process.label = 10)
kegg_chord
# 07 TMB----------
setwd("/data/nas1/luchunlin/project/BJTC-276")
if (! dir.exists("./06_TMB")){
  dir.create("./06_TMB")
}
setwd("./06_TMB")
## 在TCGA-OV数据库中，对得到的EMT与免疫相关差异基因进行突变分析。
library(maftools)
mut_ov<-read_tsv(file = 'TCGA-OV.varscan2_snv.tsv')
head(mut_ov)   
colnames(mut_ov) =c( "Tumor_Sample_Barcode", "Hugo_Symbol", 
                  "Chromosome", "Start_Position", 
                  "End_Position", "Reference_Allele", "Tumor_Seq_Allele2", 
                  "HGVSp_Short" , 'effect' ,"Consequence",
                  "vaf" ) 
mut_ov$Entrez_Gene_Id =1
mut_ov$Center ='ucsc'
mut_ov$NCBI_Build ='GRCh38'
mut_ov$NCBI_Build ='GRCh38'
mut_ov$Strand ='+'
mut_ov$Variant_Classification = mut_ov$effect
tail(sort(table(mut_ov$Variant_Classification )))
mut_ov$Tumor_Seq_Allele1 = mut_ov$Reference_Allele
mut_ov$Variant_Type = ifelse(
  mut_ov$Reference_Allele %in% c('A','C','T','G') & mut_ov$Tumor_Seq_Allele2 %in% c('A','C','T','G'),
  'SNP','INDEL'
)
mut_ov2<-mut_ov[mut_ov$Hugo_Symbol%in%rownames(DE_emtmod),]
table(mut_ov$Variant_Type )
table(mut_ov2$Hugo_Symbol)
tcga.ov = read.maf(maf = mut_ov2,
                     vc_nonSyn=names(tail(sort(table(mut_ov$Variant_Classification )))))
oncoplot(maf = tcga.ov,top = 30,fontSize = 0.7)
# 08 聚类----------
setwd("/data/nas1/luchunlin/project/BJTC-276")
if (! dir.exists("./07_cluster")){
  dir.create("./07_cluster")
}
setwd("./07_cluster")
## 08-1 聚类-----
##根据EMT与免疫相关差异基因的表达情况，将TCGA-OV样本聚类。
## 提取表达矩阵
cluster_exp<-dat_tcga[rownames(DE_emtmod),]
cluster_exp<-log2(cluster_exp+1)
cluster_exp<-as.matrix(cluster_exp)
cluster_exp<-cluster_exp[,colnames(cluster_exp)%in%survival$sample]
library(ConsensusClusterPlus)
#sweep函数减去中位数进行标准化
df<-sweep(cluster_exp,1, apply(cluster_exp,1,median,na.rm=T))
maxK <-  9 #最多分成几组
results <-  ConsensusClusterPlus(df,
                                 maxK = maxK,
                                 reps = 1000,              # 抽样次数(一般1000或更多)
                                 pItem = 0.8,             # 抽样比例
                                 pFeature = 1,
                                 clusterAlg = "pam",      # 聚类方法
                                 seed = 100,
                                 title="consensus",
                                 innerLinkage="complete",
                                 plot="pdf")

icl = calcICL(results,
              title="consensus",
              plot="pdf")

## 筛选最佳聚类数
### 一致性矩阵热图白色块最干净，尽量不掺杂蓝色
### 累积分布曲线下降的坡度最平缓
### delta area 曲线的肘部点横坐标
### 聚类一致性直方图又高又平均

### 可以使用PAC标准进行筛选
Kvec = 2:maxK
x1 = 0.1; x2 = 0.9 # threshold defining the intermediate sub-interval
PAC = rep(NA,length(Kvec))
names(PAC) = paste("K=",Kvec,sep="") # from 2 to maxK

for(i in Kvec){
  M = results[[i]]$consensusMatrix
  Fn = ecdf(M[lower.tri(M)])
  PAC[i-1] = Fn(x2) - Fn(x1)
}#end for i

# The optimal K
optK = Kvec[which.min(PAC)]
optK
## [3]
cluster<-results[[3]]$consensusClass
cluster
## 不同cluster的热图
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(gplots)
cluster_rt<-as.data.frame(cluster)
cluster_rt$cluster<-ifelse(cluster_rt$cluster=='1','cluster 1',
                           ifelse(cluster_rt$cluster=='2','cluster 2','cluster 3'))
cluster_rt$sample<-rownames(cluster_rt)
cluster_rt<-cluster_rt[order(cluster_rt$cluster),]
heat_cluster<-cluster_exp[rownames(cluster_exp)%in%rownames(DE_emtmod),]
heat_cluster<-heat_cluster[,cluster_rt$sample]
cluster_rt<-cluster_rt[,1]
cluster_rt<-as.data.frame(cluster_rt)
rownames(cluster_rt)<-colnames(heat_cluster)
colnames(cluster_rt)<-'group'
x<-log2(heat_cluster+1)
#x<-t(scale(t(heat_cluster)))
ann_colors<-list(
  group = c('cluster 1'="lightblue",'cluster 2'="darkorange",'cluster 3'='pink'))
pheatmap(mat=x,
         annotation_col = cluster_rt,
         color=bluered(100),
         scale = "row",
         annotation_colors = ann_colors,
         fontsize = 8,
         show_colnames = F,
         show_rownames = F,
         cluster_cols = F,
         cluster_rows = T)

## 08-2 PCA分析--------
## OV样本主成分分析
cluster_exp<-t(cluster_exp)
pca1<-prcomp(cluster_exp[,-ncol(cluster_exp)],center = TRUE,scale. = TRUE)
df1<-pca1$x  ## 提取PC score
df1 <- as.data.frame(df1)
summ1 <- summary(pca1)
summ1
# 提取主成分的方差贡献率,生成坐标轴标题
summ1 <- summary(pca1)
summ1
xlab1 <- paste0("PC1(",round(summ1$importance[2,1]*100,2),"%)")
ylab1 <- paste0("PC2(",round(summ1$importance[2,2]*100,2),"%)")
# 绘制PCA得分图
cluster<-as.data.frame(cluster)
cluster$cluster<-ifelse(cluster$cluster=='1','cluster 1',
                        ifelse(cluster$cluster=='2','cluster 2','cluster 3'))
library(ggplot2)
p.pca1 <- ggplot(data = df1,aes(x = PC1,y = PC2,color = cluster$cluster))+
  stat_ellipse(aes(fill = cluster$cluster),
               type = "norm",geom = "polygon",alpha = 0.25,color = NA)+ # 添加置信椭圆
  geom_point(size = 3.5)+
  labs(x = xlab1,y = ylab1,color = "Condition",title = "PCA Scores Plot")+
  guides(fill = "none")+
  theme_bw()+
  scale_fill_manual(values = c("purple","orange","pink"))+
  scale_colour_manual(values = c("purple","orange","pink"))+
  theme(plot.title = element_text(hjust = 0.5,size = 15),
        axis.text = element_text(size = 11),axis.title = element_text(size = 13),
        legend.text = element_text(size = 11),legend.title = element_text(size = 13),
        plot.margin = unit(c(0.4,0.4,0.4,0.4),'cm'))
p.pca1

## 08-3 不同cluster的KM曲线-------
survival<-read_tsv(file = 'TCGA-OV.survival.tsv')       ## 在下面一个部分
survival<-survival[survival$sample%in%colnames(dat_tcga),]
cluster_dat<-as.data.frame(cluster_exp)
#cluster_dat<-as.data.frame(t(dat_tpm[rownames(dat_tpm)%in%rownames(DE_emtmod),]))
#cluster_dat<-cluster_dat[rownames(cluster_dat)%in%survival$sample,]
cluster_dat$cluster<-as.vector(cluster$cluster,)
cluster_dat$sample<-rownames(cluster_dat)
cluster_dat<-merge(survival,cluster_dat,by='sample')
rownames(cluster_dat)<-cluster_dat$sample
cluster_dat<-cluster_dat[,-c(1,3)]

kmfit<-survfit(Surv(OS.time, OS) ~ cluster, data =  cluster_dat)
cluster_survival_median <- ggsurvplot(kmfit,
                                    pval = TRUE, 
                                    conf.int = F,
                                    legend.labs=c("cluster 1","cluster 2","cluster 3" ),
                                    legend.title="cluster",
                                    title="Train KM",
                                    font.main = c(15,"bold"),
                                    risk.table = TRUE, 
                                    risk.table.col = "strata", 
                                    linetype = "strata", 
                                    surv.median.line = "hv", 
                                    ggtheme = theme_bw(), 
                                    palette = c("#A73030FF", "#0073C2FF","purple"))
cluster_survival_median
## 08-4 不同cluster的ESTIMATE评分------
group_cluster<-cluster
group_cluster_estimate<-group_cluster$cluster%>%as.factor()
design<-model.matrix(~0 + group_cluster_estimate)
rownames(design)<-rownames(group_cluster)
colnames(design)<-levels(group_cluster_estimate)
design<-as.data.frame(design)
Cluster1<-rownames(design)[which(design$`cluster 1`==1)]
Cluster2<-rownames(design)[which(design$`cluster 2`==1)]
Cluster3<-rownames(design)[which(design$`cluster 3`==1)]

length(Cluster1)
## 121
length(Cluster2)
## 154
length(Cluster3)
## 103
#install.packages("estimate", repos="http://R-Forge.R-project.org")
library(estimate)
dat_estimate<-log2(t(cluster_exp+1))
write.table(dat_estimate, 
            'dat_estimate_log2.txt', 
            col.names = T, 
            row.names = T, 
            quote = F, sep="\t")
filterCommonGenes(input.f = './dat_estimate_log2.txt', 
                  output.f = 'dat_estimate.gct', 
                  id = 'GeneSymbol')
estimateScore('dat_estimate.gct', 'dat_purity.gct', platform="affymetrix")
# [1] "1 gene set: StromalSignature  overlap= 112"
# [1] "2 gene set: ImmuneSignature  overlap= 140"
es_score <- read.table('dat_purity.gct', skip = 2, header = T)
immu_score <- es_score[,3:length(es_score)]
rownames(immu_score) <- es_score$NAME
write.table(es_score,
            file = "es_score.xls",
            sep = "\t",
            quote = F,
            row.names = F)

## 小提琴图--------
violin_dat <- data.frame(t(immu_score))
rownames(violin_dat)<-rownames(group_cluster)
violin_dat$sample <- rownames(violin_dat)
violin_dat$group <- ifelse(violin_dat$sample %in% Cluster1,
                           "cluster 1",
                           ifelse(violin_dat$sample%in%Cluster2,'cluster 2','cluster 3'))
library(tidyr)
library(ggplot2)
library(ggpubr)
library(Ipaper)
my_comparisons <- list(c("cluster 1","cluster 2"),c("cluster 1","cluster 3"),c("cluster 2","cluster 3"))

{

    p1 <- ggplot(violin_dat, aes(x=group,y=StromalScore, fill=group))+
    geom_violin() + #绘制小提琴图
    stat_boxplot(geom="errorbar",
                 width=0.1,
                 position = position_dodge(0.9)) +
    geom_boxplot(width=0.4,
                 position=position_dodge(0.9),
                 outlier.shape = NA,
                 fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
    geom_point(aes(fill = group),
               size = 0.5,
               position = position_dodge(0.9))+
    scale_fill_manual(values = c("#B0E0E6", "#FF8888","orange"), name = "Group") + 
    geom_signif(comparisons = my_comparisons,
                map_signif_level = T,
                y_position = c(40,60,50))+
  #  stat_compare_means(aes(group = group),method = 'wilcox')+ 
    theme_bw()+ #背景变为白色
    ylim(-50,70) +
    theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
          panel.grid = element_blank()) + 
    theme(axis.text.x=element_text(family="Times",size=12,face="bold"), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
          axis.text.y=element_text(family="Times",size=12,face="bold"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
          axis.title.y=element_text(family="Times",size = 15,face="bold"), #设置y轴标题的字体属性
          axis.title.x=element_text(family="Times",size = 10,face="bold"),
          plot.title = element_text(hjust = 0.5,family="Times",size = 15,face="bold"),
          legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold"))+
    labs(title="Stromal Score", x="", y="Stromal Score")
  p1
  p2 <- ggplot(violin_dat, aes(x=group,y=ImmuneScore, fill=group))+
    geom_violin() + #绘制小提琴图
    stat_boxplot(geom="errorbar",
                 width=0.1,
                 position = position_dodge(0.9)) +
    geom_boxplot(width=0.4,
                 position=position_dodge(0.9),
                 outlier.shape = NA,
                 fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
    geom_point(aes(fill = group),
               size = 0.5,
               position = position_dodge(0.9))+
    scale_fill_manual(values = c("#B0E0E6", "#FF8888",'orange'), name = "Group") + 
    geom_signif(comparisons = my_comparisons,
                map_signif_level = T,
                y_position = c(40,60,50))+
  #  stat_compare_means(aes(group = group),method = 'wilcox')+ 
    theme_bw()+ #背景变为白色
    ylim(-50,70) +
    theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
          panel.grid = element_blank()) + 
    theme(axis.text.x=element_text(family="Times",size=12,face="bold"), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
          axis.text.y=element_text(family="Times",size=12,face="bold"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
          axis.title.y=element_text(family="Times",size = 15,face="bold"), #设置y轴标题的字体属性
          axis.title.x=element_text(family="Times",size = 10,face="bold"),
          plot.title = element_text(hjust = 0.5,family="Times",size = 15,face="bold"),     
          legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold"))+
    labs(title="Immune Score", x="", y="Immune Score")
  p2
  p3 <- ggplot(violin_dat, aes(x=group, y=ESTIMATEScore, fill=group))+
    geom_violin() + #绘制小提琴图
    stat_boxplot(geom="errorbar",
                 width=0.1,
                 position = position_dodge(0.9)) +
    geom_boxplot(width=0.4,
                 position=position_dodge(0.9),
                 outlier.shape = NA,
                 fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
    geom_point(aes(fill = group),
               size = 0.5,
               position = position_dodge(0.9))+
    scale_fill_manual(values = c("#B0E0E6", "#FF8888",'orange'), name = "Group") + 
    geom_signif(comparisons = my_comparisons,
                map_signif_level = T,
                y_position = c(40,60,50))+
  #  stat_compare_means(aes(group = group),method = 'wilcox')+ 
    theme_bw()+ #背景变为白色
    theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
          panel.grid = element_blank()) + 
    ylim(-90, 70) +
    theme(axis.text.x=element_text(family="Times",size=12,face="bold"), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
          axis.text.y=element_text(family="Times",size=12,face="bold"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
          axis.title.y=element_text(family="Times",size = 15,face="bold"), #设置y轴标题的字体属性
          axis.title.x=element_text(family="Times",size = 10,face="bold"),
          plot.title = element_text(hjust = 0.5,family="Times",size = 15,face="bold"),
          legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold"))+
    labs(title="ESTIMATE Score", x="", y="ESTIMATE Score")
  p3
  p4 <- ggplot(violin_dat, aes(x=group, y=TumorPurity, fill=group))+
    geom_violin() + #绘制小提琴图
    stat_boxplot(geom="errorbar",
                 width=0.1,
                 position = position_dodge(0.9)) +
    geom_boxplot(width=0.4,
                 position=position_dodge(0.9),
                 outlier.shape = NA,
                 fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
    geom_point(aes(fill = group),
               size = 0.5,
               position = position_dodge(0.9))+
    scale_fill_manual(values = c("#B0E0E6", "#FF8888",'orange'), name = "Group") + 
    geom_signif(comparisons = my_comparisons,
                map_signif_level = T,
                y_position = c(0.835,0.845,0.84))+
   # stat_compare_means(aes(group = group),method = 'wilcox')+ 
    theme_bw()+ #背景变为白色
    theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
          panel.grid = element_blank()) + 
    ylim(0.815, 0.85) +
    theme(axis.text.x=element_text(family="Times",size=12,face="bold"), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
          axis.text.y=element_text(family="Times",size=12,face="bold"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
          axis.title.y=element_text(family="Times",size = 15,face="bold"), #设置y轴标题的字体属性
          axis.title.x=element_text(family="Times",size = 10,face="bold"),
          plot.title = element_text(hjust = 0.5,family="Times",size = 15,face="bold"),
          legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold"))+
    labs(title="Tumor Purity", x="", y="Tumor Purity")
  p4
  p5 <- cowplot::plot_grid(p1,p2,p3,p4,
                           nrow = 2, 
                           align = 'h', 
                           vjust = -0.3)
  p5
}


## 08-5 不同cluster的ssGSEA--------
library(GSVA)
gene_set <- read.table("/data/nas1/luchunlin/pipeline/ssGSEA/mmc3.txt",
                       header = T,
                       sep ="\t")
gene_list <- split(as.matrix(gene_set)[,1],
                   gene_set[,2])
dat.tcga<-as.matrix(log2(dat_tcga+1))
dat.tcga<-dat.tcga[,colnames(dat.tcga)%in%survival$sample]
ssgsea_score = gsva(dat.tcga, gene_list, 
                    method = "ssgsea", 
                    ssgsea.norm = TRUE, 
                    verbose = TRUE)
write.table(ssgsea_score,
            file = "ssgsea_result.xls",
            sep = "\t",
            quote = F)
tiics_result <- ssgsea_score

### 画热图（结合临床信息）
clinical<-read_tsv(file ='TCGA-OV.GDC_phenotype.tsv')
phenotype<-data.frame(sample=clinical$submitter_id.samples,
                      stage=clinical$clinical_stage,
                      age=clinical$age_at_index.demographic)

phenotype$stage<-gsub('A','',phenotype$stage)
phenotype$stage<-gsub('B','',phenotype$stage)
phenotype$stage<-gsub('C','',phenotype$stage)
phenotype$age<-cut(phenotype$age,breaks = c(30,40,50,60,70,80,90),labels = c('30-40','40-50','50-60','60-70','70-80','80-90'))

gsea_group<-data.frame(sample=colnames(ssgsea_score),
                       group=cluster$cluster)
gsea_group<-merge(gsea_group,phenotype,by='sample')

gsea_group<-gsea_group[order(gsea_group$group),]
gsea_rt_dat<-ssgsea_score[,gsea_group$sample]
rownames(gsea_group)<-gsea_group$sample
gsea_group<-dplyr::select(gsea_group,c('group','stage','age'))
gsea_group$group<-as.factor(gsea_group$group)
gsea_group$stage<-as.factor(gsea_group$stage)
gsea_group$age<-as.factor(gsea_group$age)
## gsea_rt_dat<-gsea_rt_dat[diff_Table$immune_cell,]
#rownames(gsea_rt_dat)<-diff_Table$sig
gsea_rt_dat<-log2(gsea_rt_dat+1)
ann_colors<-list(
  group=c('cluster 1'='#FFB7DD','cluster 2'='#77DDFF','cluster 3'='orange'),
  stage=c('Stage I'='#FF8C00','Stage II'='#20B2AA','Stage III'='#40E0D0','Stage IV'='#9ACD32'),
  age=c('30-40'='#33CCCC','40-50'='#CCCCFF','50-60'='#FF9966','60-70'='#FFCCFF','70-80'='#DDA0DD','80-90'='#FF69B4')
)
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(gplots)
pheatmap(gsea_rt_dat,
         color = bluered(100),
         border_color = NA,
         annotation_col = gsea_group,
         annotation_colors = ann_colors,
         labels_row = NULL,
         clustering_method = 'ward.D2',
         show_rownames = T,
         show_colnames = F,
         fontsize_col = 5,
         cluster_cols = F,
         cluster_rows = T)

#09 单因素cox----------
setwd("/data/nas1/luchunlin/project/BJTC-276")
if (! dir.exists("./08_univariate_cox")){
  dir.create("./08_univariate_cox")
}
setwd("./08_univariate_cox")
## 匹配生存数据
survival<-read_tsv(file = 'TCGA-OV.survival.tsv')
survival<-survival[survival$sample%in%colnames(dat_tcga),]
## 378个匹配
train_dat<-t(dat_tpm[,survival$sample])%>%as.data.frame()
#train_dat<-log2(train_dat+0.001)
#train_dat<-scale(train_dat)
train_dat<-train_dat[,rownames(DE_emtmod)]
train_dat$sample<-rownames(train_dat)
train_dat<-merge(survival,train_dat,by='sample')
rownames(train_dat)<-train_dat$sample
train_dat<-train_dat[,-c(1,3)]
### 按照生存状态（OS）进行三七分组。七分作为训练集，三分作为验证集。
library(caret)
set.seed(24)
###8 24---12
{all.expr<-train_dat
all.expr$sample<-rownames(all.expr)
sam<-createDataPartition(all.expr$OS.time,p= .7,list = F)
train_sample<-all.expr$sample[sam]
train_sample<-as.data.frame(train_sample)
test_sample<-all.expr$sample[-sam]
test_sample<-as.data.frame(test_sample)
train_data<-train_dat[rownames(train_dat)%in%train_sample$train_sample,]
test_data<-train_dat[rownames(train_dat)%in%test_sample$test_sample,]
dim(train_data)
##[1] 267  
dim(test_data)
##[1] 111 

### 单因素cox
library(survival)
library(survminer)
colnames_sum <- colnames(train_data)
#colnames_sum <- gsub("-","_",colnames_sum)
#colnames_sum <- gsub(" ","_",colnames_sum)
#colnames(train_data) <- colnames_sum
covariates <- colnames_sum[-which(colnames_sum %in% c("OS", "OS.time"))]
#Surv()函数产生一个生存对象  生存时间对生存的影响 对每一个变量构建生存分析公式
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste("Surv(OS.time, OS)~", x)))  #as.formula(). 将字符串转换成公式。构建formula对象
# coxph函数用于计算cox模型 循环对每一个特征做cox回归分析
univ_models <- lapply(univ_formulas,
                      function(x) {coxph(x, data = train_data)})
univ_models$KRT7
univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         #获取p值
                         p.value<-signif(x$wald["pvalue"],digits=3)
                         #获取HR
                         HR <-signif(x$coef[2], digits=4);
                         #获取95%置信区间
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"],4)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],4)
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
}
res_results_0.05 <- na.omit(res_results_0.05)

{
write.table(res_results_0.05,
            file = "univariate_cox_result_0.05.xls",
            quote = F,
            sep = '\t',
            row.names = T)
dim(res_results_0.05)
## 13 2
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
hz <- paste(round(res_results_0.05_2$HR,4),
            "(",round(res_results_0.05_2$HR.95L,4),
            "-",round(res_results_0.05_2$HR.95H,4),")",sep = "")


tabletext <- cbind(c(NA,"Gene",rownames(res_results_0.05_2)),
                   c(NA,"P value",ifelse(res_results_0.05_2$p.value<0.001,
                                         "< 0.001",
                                         round(res_results_0.05_2$p.value,3))),
                   c(NA,"Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1
library(forestplot)
pdf(file = "univariate_cox_forest.pdf", height = 55, width = 10, onefile = F)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE, TRUE,rep(FALSE, 300)),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,NA,res_results_0.05_2$HR),
           lower=c(NA,NA,res_results_0.05_2$HR.95L), #95%置信区间下限
           upper=c(NA,NA,res_results_0.05_2$HR.95H), #95%置信区间上限
           boxsize=0.2,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(0.8,1,1.2,1.4), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(0.8,"cm"), #固定行高
           graphwidth = unit(.5,"npc"), #图在表中的宽度比例
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines = list("3" = gpar(col = "black", lty = 1, lwd = 2)),
           # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
           #                 "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
           #                 "59" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
           # mar=unit(rep(0.01, times = 4), "cm"),#图形页边距
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=0.8),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Hazard Ratio",
           grid = T) # 垂直于x轴的网格线，对应每个刻度
dev.off()
png(filename = "univariate_cox_forest.png", height = 4000, width = 700)
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
           xticks = c(0, 0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(0.5,"cm"), #固定行高
           graphwidth = unit(.5,"npc"), #图在表中的宽度比例
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines = list("3" = gpar(col = "black", lty = 1, lwd = 2)),
           # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
           #                 "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
           #                 "59" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
           # mar=unit(rep(0.01, times = 4), "cm"),#图形页边距
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=0.8),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Hazard Ratio",
           grid = T) # 垂直于x轴的网格线，对应每个刻度

dev.off()

# 10 lasso----------
setwd("/data/nas1/luchunlin/project/BJTC-276")
if (! dir.exists("./09_lasso")){
  dir.create("./09_lasso")
}
setwd("./09_lasso")
library(glmnet)
x_all <- subset(train_data, select = -c(OS, OS.time))
x_all <- x_all[,rownames(res_results_0.05)]
y_all <- subset(train_data, select = c(OS, OS.time))

# 拟合模型
fit <- glmnet(as.matrix(x_all), Surv(y_all$OS.time,y_all$OS), 
              family = "cox") 

png(filename = "lasso_model.png", height = 400, width = 500)
plot(fit, xvar = "lambda",label = TRUE, las=1)
dev.off()
pdf(file = "lasso_model.pdf", height = 5)
plot(fit, xvar = "lambda",label = TRUE, las=1)
dev.off()
# 交叉验证拟合模型
##4  5
}
set.seed(5)
{
cvfit = cv.glmnet(as.matrix(x_all),
                  Surv(y_all$OS.time,y_all$OS),nfold=10,
                  family = "cox") 

png(filename = "lasso_verify.png", height = 400, width = 500)
plot(cvfit, las =1)
dev.off()
pdf(file = "lasso_verify.pdf", height = 5)
plot(cvfit, las =1)
dev.off()
cvfit
# 提取指定lambda时特征的系数
coef.min = coef(cvfit, s = "lambda.min")  ## lambda.min & lambda.1se 取一个
cvfit$lambda.min
coef.min
# [1] 0.02683339
# 找出那些回归系数没有被惩罚为0的
active.min = which(coef.min@i != 0)
# 提取基因名称
lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i+1]
}
lasso_geneids
## 8
##   "EPB41L3"  "KRT7"     "QKI"      "SERPINF1" "ENO2"     "PLXND1"   "CXCL12"   "IGFBP7"  
{
write(lasso_geneids, "lasso_genes.csv")
write.csv(x_all,file = "Lasso_x.csv",quote = F)
write.csv(y_all,file = "Lasso_y.csv",quote = F)
write.csv(y_all,file = "Lasso_y.csv",quote = F)

# 11 riskModel----------
setwd("/data/nas1/luchunlin/project/BJTC-276")
if (! dir.exists("./10_risk")){
  dir.create("./10_risk")
}
setwd("./10_risk")
## 11-1 riskscore-----

riskScore=predict(cvfit,newx = as.matrix(x_all),s=cvfit$lambda.min)
riskScore<-as.numeric(riskScore)
coxGene=lasso_geneids
outCol=c("OS","OS.time",coxGene)
mean(riskScore)
median(riskScore)
risk=as.vector(ifelse(riskScore>median(riskScore),0,1))
risk <- as.data.frame(c(cbind(id=rownames(cbind(train_data[,outCol],
                                                riskScore,
                                                risk)),
                              cbind(train_data[,outCol],
                                    riskScore,
                                    risk))))
table(risk$risk)
##  0  1 
## 133 134
risk_dis <- ggplot(risk, aes(x=reorder(id, riskScore), 
                             y=riskScore, 
                             color = factor(risk, 
                                            levels = c(0, 1), 
                                            labels = c("High Risk", "Low Risk")))) +
  geom_point() +
  scale_color_manual(values = c("#A73030FF", "#0073C2FF")) + 
  scale_x_discrete(breaks = risk[order(risk$riskScore),]$id[c(1,40,80,120,160,200,240,280)],
                   labels = c(1,40,80,120,160,200,240,280),
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
        legend.margin = margin(c(-5,4,4,3)),
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
  scale_x_discrete(breaks = risk[order(risk$riskScore),]$id[c(1,40,80,120,160,200,240,280)],
                   labels = c(1,40,80,120,160,200,240,280),
                   expand = c(0.02,0)) +
  ylim(c(0,25))+
  geom_vline(xintercept = nrow(risk[which(risk$risk==1),]) + 0.5,
             lty = 2) +
  labs(x = "Patients(increasing risk score)",
       y = "Survival time (years)",
       title = "Train survival state distribution") + 
  theme_base() +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.justification = c(0,1),
        legend.position = c(0,1),
        legend.margin = margin(c(-5,4,4,3)),
        legend.background = element_rect(color = "black", size = .3),
        plot.title = element_text(size = 15, hjust = 0.5))

surv_stat

## 11-2 KM曲线----------
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
}
## 11-3 ROC曲线-------- 
# BiocManager::install('survivalROC')
{
library(survivalROC)
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
for_multi_ROC <- multi_ROC(time_vector = c(365*seq(5,9,2)), 
                           risk_score_table = risk)
for_multi_ROC$Time <- factor(for_multi_ROC$Time)
#devtools::install_github('yikeshu0611/geomROC')
library(scales)
library(geomROC)
library(plotROC)
library(ggthemes)

auc_y9<- round(for_multi_ROC[which(for_multi_ROC$Time==3285),5][1],2)
auc_y5 <- round(for_multi_ROC[which(for_multi_ROC$Time==1825),5][1],2)
auc_y7 <- round(for_multi_ROC[which(for_multi_ROC$Time==2555),5][1],2)
ROC <- ggplot(for_multi_ROC, aes(x=False_positive,
                                 y=True_positive, 
                                 label=Cut_values, 
                                 color=Time)) + 
  scale_color_manual(breaks = c("1825","2555","3285"),
                     labels = c( "5 years", "7 years","9 years"),
                     values = c("#6495ED", "#FF4040", "#20B2AA")) +
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
           label = c(paste('AUC of 5 years =', format(auc_y5,nsmall=2)),
                     paste('AUC of 7 years =', format(auc_y7,nsmall=2)),
                     paste('AUC of 9 years =', format(auc_y9,nsmall=2))))
ROC
}
{
ggsave('Train ROC.png', ROC,width = 5, height = 4)
ggsave('Train ROC.pdf', ROC,width = 5, height = 4)


# 12 验证集----------
setwd("/data/nas1/luchunlin/project/BJTC-276")
if (! dir.exists("./11_External_va")){
  dir.create("./11_External_va")
}
setwd("./11_External_va")

library(GEOquery)
library(Biobase)
# 开始验证
##手动算
test_data2<-test_data[,lasso_geneids]
risk_out<-data.frame(test_data2)
risk_out$risk_out<-NA
risk_out$riskScore_out<-NA
cnt<-1

coef.min<-coef.min[lasso_geneids,]

while (cnt < 112) {
  risk_out$riskScore_out[cnt]<-sum(coef.min*test_data2[cnt,])
  cnt = cnt + 1
}

dim(test_data2)
cnt<-1
while (cnt < 112) {
  risk_out$risk_out[cnt]=as.vector(ifelse(risk_out$riskScore_out[cnt]>median(risk_out$riskScore_out),0,1))
  cnt = cnt + 1
}
riskScore_out<-as.numeric(risk_out$riskScore_out)

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
table(risk_out$risk_out)
## 55 56
risk_dis_out <- ggplot(risk_out, aes(x=reorder(id, riskScore_out), 
                                     y=riskScore_out, 
                                     color = factor(risk_out, 
                                                    levels = c(0, 1), 
                                                    labels = c("High Risk", "Low Risk")))) +
  geom_point() +
  scale_color_manual(values = c("#A73030FF", "#0073C2FF")) + 
  scale_x_discrete(breaks = risk_out[order(risk_out$riskScore_out),]$id[c(1,25,50,75,100,125,150,175,200)],
                   labels = c(1,25,50,75,100,125,150,175,200),
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
        legend.margin = margin(c(-5,4,4,3)),
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
  scale_x_discrete(breaks = risk_out[order(risk_out$riskScore_out),]$id[c(1,25,50,75,100,125,150,175,200)],
                   labels = c(1,25,50,75,100,125,150,175,200),
                   expand = c(0.02,0)) +
  ylim(x=c(0,25)) +
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
        legend.margin = margin(c(-5,4,4,3)),
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
}
verify_survival_median

{
  library(survivalROC)
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
  for_multi_ROC <- multi_ROC(time_vector = c(365*seq(5,9,2)), 
                             risk_score_table = risk_out)
  for_multi_ROC$Time <- factor(for_multi_ROC$Time)
  #devtools::install_github('yikeshu0611/geomROC')
  library(scales)
  library(geomROC)
  library(plotROC)
  library(ggthemes)
  
  auc_y9<- round(for_multi_ROC[which(for_multi_ROC$Time==3285),5][1],2)
  auc_y5 <- round(for_multi_ROC[which(for_multi_ROC$Time==1825),5][1],2)
  auc_y7 <- round(for_multi_ROC[which(for_multi_ROC$Time==2555),5][1],2)
  ROC <- ggplot(for_multi_ROC, aes(x=False_positive,
                                   y=True_positive, 
                                   label=Cut_values, 
                                   color=Time)) + 
    scale_color_manual(breaks = c("1825","2555","3285"),
                       labels = c( "5 years", "7 years","9 years"),
                       values = c("#6495ED", "#FF4040", "#20B2AA")) +
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
             label = c(paste('AUC of 5 years =', format(auc_y5,nsmall=2)),
                       paste('AUC of 7 years =', format(auc_y7,nsmall=2)),
                       paste('AUC of 9 years =', format(auc_y9,nsmall=2))))
  ROC
}
ggsave('Test ROC.png', ROC,width = 5, height = 4)
ggsave('Test ROC.pdf', ROC,width = 5, height = 4)



# 13 细胞表达量分析（风险评分相关基因在卵巢癌样本上皮细胞的表达量分析）----------
setwd("/data/nas1/luchunlin/project/BJTC-276")
if (! dir.exists("./12_scRNA")){
  dir.create("./12_scRNA")
}
setwd("./12_scRNA")
Biocductor_packages <- c("Seurat",
                         "scran",
                         "scater",
                         "monocle",
                         "DropletUtils",
                         "SingleR",
                         "monocle3"
)

options(BioC_mirror="https://mirrors.ustc.edu.cn/bioc/")
#options(BioC_mirror="https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
# use BiocManager to install
for (pkg in Biocductor_packages){
  if (! require(pkg, character.only = T)){
    BiocManager::install(pkg, ask = F, update = F)
    require(pkg)
  }
}

#最后再检查下成功与否
for (pkg in Biocductor_packages){
  require(pkg, character.only = T)
}
library(tidyverse)
library(Ipaper)
library(patchwork)

# 安装monocle3
bioc_pkgs <- c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
               'limma', 'S4Vectors', 'SingleCellExperiment',
               'SummarizedExperiment', 'batchelor', 'Matrix.utils')
for (bpkg in bioc_pkgs){
  if (! require(bpkg,character.only=T) ) {
    BiocManager::install(bpkg,ask = F,update = F)
    require(bpkg,character.only=T)
  }
}
library(monocle3)
## GSE165897   GSE173682 GSE147082 
## 导入数据，创建seurat对象
## GSE189843
dir='GSE189843_RAW/GSM5708485_RAW/'
list.files(dir)
GSM85<-Read10X(data.dir = 'GSE189843_RAW/GSM5708485_RAW/')
GSM86<-Read10X(data.dir = 'GSE189843_RAW/GSM5708486_RAW/')
GSM87<-Read10X(data.dir = 'GSE189843_RAW/GSM5708487_RAW/')
GSM88<-Read10X(data.dir = 'GSE189843_RAW/GSM5708488_RAW/')
GSM89<-Read10X(data.dir = 'GSE189843_RAW/GSM5708489_RAW/')
GSM90<-Read10X(data.dir = 'GSE189843_RAW/GSM5708490_RAW/')
GSM91<-Read10X(data.dir = 'GSE189843_RAW/GSM5708491_RAW/')
GSM92<-Read10X(data.dir = 'GSE189843_RAW/GSM5708492_RAW/')
GSM93<-Read10X(data.dir = 'GSE189843_RAW/GSM5708493_RAW/')
GSM94<-Read10X(data.dir = 'GSE189843_RAW/GSM5708494_RAW/')
GSM95<-Read10X(data.dir = 'GSE189843_RAW/GSM5708495_RAW/')
GSM96<-Read10X(data.dir = 'GSE189843_RAW/GSM5708496_RAW/')

count1<-CreateSeuratObject(counts = GSM85,min.cells = 3,min.features = 200)
count2<-CreateSeuratObject(counts = GSM86,min.cells = 3,min.features = 200)
count3<-CreateSeuratObject(counts = GSM87,min.cells = 3,min.features = 200)
count4<-CreateSeuratObject(counts = GSM88,min.cells = 3,min.features = 200)
count5<-CreateSeuratObject(counts = GSM89,min.cells = 3,min.features = 200)
count6<-CreateSeuratObject(counts = GSM90,min.cells = 3,min.features = 200)
count7<-CreateSeuratObject(counts = GSM91,min.cells = 3,min.features = 200)
count8<-CreateSeuratObject(counts = GSM92,min.cells = 3,min.features = 200)
count9<-CreateSeuratObject(counts = GSM93,min.cells = 3,min.features = 200)
count10<-CreateSeuratObject(counts = GSM94,min.cells = 3,min.features = 200)
count11<-CreateSeuratObject(counts = GSM95,min.cells = 3,min.features = 200)
count12<-CreateSeuratObject(counts = GSM96,min.cells = 3,min.features = 200)
scRNA<-merge(x=count1,y=c(count2,count3,count4,count5,count6,count7,count8,count9,count10,count11,count12),
             add.cell.ids=c('GSM5708485','GSM5708486','GSM5708487','GSM5708488','GSM5708489','GSM5708490','GSM5708491','GSM5708492','GSM5708493','GSM5708494','GSM5708495','GSM5708496'),
             merge.data=T)
#gunzip('GSE165897_UMIcounts_HGSOC.tsv.gz')
#counts<-read.table(file ='GSE165897_UMIcounts_HGSOC.tsv',header = T, row.names=1, sep="", as.is=T)
#head(counts[c(1:6)])
#scRNA <- CreateSeuratObject(counts = counts,min.cells=3,min.genes=200)
##初步过滤，>=3个细胞中表达的基因，>=200个基因的细胞
scRNA
## 13-1 质控------
scRNA[['percent.mt']]<-PercentageFeatureSet(scRNA,pattern = "^MT-")
head(scRNA@meta.data)
summary(scRNA@meta.data)
# 计算基因含量，MT为线粒体
VlnPlot(scRNA,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
## 可视化。线粒体基因含量占比5%以下的细胞保留。

## 可视化RNA-基因含量（RNA-feature）
plot1<-FeatureScatter(scRNA,feature1 = 'nCount_RNA',feature2 = 'percent.mt')
plot2<-FeatureScatter(scRNA,feature1 = 'nCount_RNA',feature2 = 'nFeature_RNA')
plot1+plot2


#去除线粒体基因表达比例过高的细胞，和一些极值细胞
scRNA <- subset(scRNA,
              subset = nFeature_RNA > 200 & nFeature_RNA < 45000 &
                percent.mt < 5)
VlnPlot(scRNA,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
#标准化
scRNA.norm<-NormalizeData(scRNA,normalization.method = "LogNormalize",scale.factor = 10000)
##scRNA.norm<-scRNA
## 13-1 PCA UMAP降维聚类-------
## 缩放数据
## 筛选1500个高变基因
all.genes<-rownames(scRNA.norm)
scRNA.norm<-FindVariableFeatures(scRNA.norm,selection.method = "vst", nfeatures = 1500)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(scRNA.norm), 10)
# plot variable features with and without labels
plot3 <- VariableFeaturePlot(scRNA.norm)
plot4 <- LabelPoints(plot = plot3, 
                     points = top10, 
                     repel = TRUE)
p3 <- plot3+plot4
p3 
plot4
library(Ipaper)
write_fig(plot4,
          file = "feature_selection.pdf",
          width = 5,
          height = 8,
          devices = NULL,
          res = 600,
          show = F)
write_fig(plot4,
          file = "feature_selection.png",
          width = 5,
          height = 8,
          devices = NULL,
          res = 300,
          show = F)

scRNA.nor.sca<-ScaleData(scRNA.norm)
## PCA降维
scRNA.norm.pca<-RunPCA(scRNA.nor.sca,features = VariableFeatures(object = scRNA.nor.sca))
print(scRNA.norm.pca[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(scRNA.norm.pca, dims = 1:2, reduction = "pca")
DimPlot(scRNA.norm.pca, reduction = "pca")
DimHeatmap(scRNA.norm.pca, dims = 1:5, cells = 500, balanced = TRUE)
## UMAP可视化
## 首先找到最佳聚类数
# 定义数据集的维度
#NOTE: This process can take a long time for big datasets, comment out for expediency. 
#More approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
scRNA.norm.pca <- JackStraw(scRNA.norm.pca, num.replicate = 100, dims = 50)
scRNA.norm.pca <- ScoreJackStraw(scRNA.norm.pca, dims = 1:50)
plot5 <- JackStrawPlot(scRNA.norm.pca, dims = 1:50)
plot5
### 27个
#T-SNE


write_fig(plot5,
          file = "pca_cluster.pdf",
          width = 15,
          height = 8,
          devices = NULL,
          res = 600,
          show = F)
write_fig(plot5,
          file = "pca_cluster.png",
          width = 15,
          height = 8,
          devices = NULL,
          res = 300,
          show = F)
plot6 <- ElbowPlot(scRNA.norm.pca, ndims = 50)
plot6
write_fig(plot6,
          file = "pca_sd.pdf",
          width = 15,
          height = 8,
          devices = NULL,
          res = 600,
          show = F)
write_fig(plot6,
          file = "pca_sd.png",
          width = 15,
          height = 8,
          devices = NULL,
          res = 300,
          show = F)


scRNA.norm.pca.c<-FindNeighbors(scRNA.norm.pca,dims = 1:27)
scRNA.norm.pca.c<-FindClusters(scRNA.norm.pca.c,resolution = 0.8)
UMAP<-RunUMAP(scRNA.norm.pca.c,dims = 1:27)
DimPlot(UMAP,reduction = 'umap')

#-tsne-----
set.seed(123)
TSNE <- RunTSNE(scRNA.norm.pca.c, dims = 1:27)
head(TSNE@reductions$tsne@cell.embeddings)
p5 <- DimPlot(TSNE, reduction = "tsne", label = T)
p5
saveRDS(TSNE,file = 'TSNE.rds')
## 13-2 细胞类型注释------
scRNA_nor_singleR <- GetAssayData(TSNE, slot = "data")
hpca.se <- HumanPrimaryCellAtlasData()
hpca.se
clusters=TSNE@meta.data$seurat_clusters
pred.hesc <- SingleR(test = scRNA_nor_singleR,
                     ref = hpca.se,
                     labels = hpca.se$label.main,
                     method = "cluster", 
                     clusters = clusters, 
                     quantile = 0.7,
                     assay.type.test = "logcounts", 
                     assay.type.ref = "logcounts")
table(pred.hesc$labels)
celltype = data.frame(ClusterID=rownames(pred.hesc), 
                      celltype=pred.hesc$labels, 
                      stringsAsFactors = F) 
TSNE@meta.data$singleR=celltype[match(clusters,celltype$ClusterID),'celltype']
p6 <- DimPlot(TSNE, reduction = "tsne", group.by = "singleR", label = T)
p6
phe=TSNE@meta.data
table(phe$singleR)
cell_type_stat <- as.data.frame(sort(table(phe$singleR)))
write.table(cell_type_stat,
            file = "cell_type_stat.txt",
            quote = F)

write_fig(p6,
          file = "cell_type_singleR.pdf",
          width = 10,
          height = 8,
          devices = NULL,
          res = 600,
          show = F)
write_fig(p6,
          file = "cell_type_singleR.png",
          width = 10,
          height = 8,
          devices = NULL,
          res = 300,
          show = F)

hub_gene <- lasso_geneids
TSNE$singleR
p8 <- FeaturePlot(TSNE, features = hub_gene, ncol = 3,split.by = "singleR")
p8
ggsave(filename = 'exp.png',p8,width = 14,height = 18)
ggsave(filename = 'exp.pdf',p8,width = 14,height = 18)

write_fig(p8,
          file = "hub_gene_tsne.pdf",
          width = 15,
          height = 10,
          devices = NULL,
          res = 600,
          show = F)
write_fig(p8,
          file = "hub_gene_tsne.png",
          width = 15,
          height = 10,
          devices = NULL,
          res = 300,
          show = F)

## 13-3 细胞轨迹分析------
library(future)
data <- as(as.matrix(TSNE@assays$RNA@counts), 'sparseMatrix')
# count矩阵
pd <- new('AnnotatedDataFrame', data = TSNE@meta.data)
# meta表转成特定格式
fData <- data.frame(gene_short_name = row.names(data), 
                    row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
# 基因名表转成特定格式
mycds <- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd,
                        expressionFamily = negbinomial.size())
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores=4, relative_expr = TRUE)
save(mycds, file = "mycds_raw.RData")
sce.markers <- FindAllMarkers(TSNE)
all.markers = sce.markers %>% dplyr::select(gene, everything()) %>% 
  subset(p_val<0.05 & abs(sce.markers$avg_log2FC) > 0.5)

write.table(all.markers, file = "All_Markers.xls",
            quote = F, row.names = F)
markers.gene <- all.markers$gene
mycds <- setOrderingFilter(mycds, markers.gene)
#降维
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
save(mycds, file = "mycds_reduce.RData")
#耗时，耗内存
#排序
load("mycds_reduce.RData")
mycds <- orderCells(mycds)
save(mycds, file = "mycds_order.RData")
load("mycds_order.RData")
table(mycds$singleR)
mycds$singleR <- factor(mycds$singleR, levels = c("Epithelial_cells" ,
                                                  "Fibroblasts",
                                                  "Macrophage",
                                                  "MSC",
                                                  "Tissue_stem_cells"))

plot11 <- plot_cell_trajectory(mycds, color_by = "seurat_clusters") +
  theme(legend.position = "right") +
  theme(legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 15))

plot11
write_fig(plot11,
          file = "sub_cell_trajectory_raw.pdf",
          width = 10,
          height = 8,
          devices = NULL,
          res = 600,
          show = F)
write_fig(plot11,
          file = "sub_cell_trajectory_raw.png",
          width = 10,
          height = 8,
          devices = NULL,
          res = 300,
          show = F)

plot12 <- plot_cell_trajectory(mycds, color_by = "singleR") +
  theme(legend.position = "right") +
  scale_color_discrete(breaks = c("Epithelial_cells" ,
                                  "Fibroblasts",
                                  "Macrophage",
                                  "MSC",
                                  "Tissue_stem_cells"))+
  theme(legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 15))
plot12


write_fig(plot12,
          file = "sub_cell_trajectory_singleR.pdf",
          width = 10,
          height = 8,
          devices = NULL,
          res = 600,
          show = F)
write_fig(plot12,
          file = "sub_cell_trajectory_singleR.png",
          width = 10,
          height = 8,
          devices = NULL,
          res = 300,
          show = F)

plot13 <- FeaturePlot(TSNE, features = hub_gene, ncol = 3)
plot13
write_fig(plot13,
          file = "subcell_hub_gene_tsne.pdf",
          width = 15,
          height = 10,
          devices = NULL,
          res = 600,
          show = F)
write_fig(plot13,
          file = "subcell_hub_gene_tsne.png",
          width = 15,
          height = 10,
          devices = NULL,
          res = 300,
          show = F)

plot14 <- plot_cell_trajectory(mycds, color_by = "State") +
  theme(legend.position = "right")+
  theme(legend.title = element_text(size = 15, face = "bold"),
        legend.text = element_text(size = 15))
plot14


write_fig(plot14,
          file = "sub_cell_trajectory_state.pdf",
          width = 10,
          height = 8,
          devices = NULL,
          res = 600,
          show = F)
write_fig(plot14,
          file = "sub_cell_trajectory_state.png",
          width = 10,
          height = 8,
          devices = NULL,
          res = 300,
          show = F)
# all_trajectory <- plot11 + plot12 + plot14 + plot_layout(ncol = 2)
# all_trajectory
# write_fig(all_trajectory,
#           file = "sub_cell_trajectory_all.pdf",
#           width = 12,
#           height = 8,
#           devices = NULL,
#           res = 600,
#           show = F)
# write_fig(all_trajectory,
#           file = "sub_cell_trajectory_all.png",
#           width = 12,
#           height = 8,
#           devices = NULL,
#           res = 300,
#           show = F)
# plot15 <- plot_cell_trajectory(mycds, color_by = "singleR") +
#   theme(legend.position = "right") +
#   facet_wrap(~singleR, nrow =1)+
#   scale_color_discrete(breaks = c("Epithelial_cells" ,
#                                   "Tissue_stem_cells",
#                                   "Endothelial_cells"))
# plot15
# 
# 
# write_fig(plot15,
#           file = "sub_cell_trajectory_facet.pdf",
#           width = 10,
#           height = 8,
#           devices = NULL,
#           res = 600,
#           show = F)
# write_fig(plot15,
#           file = "sub_cell_trajectory_facet.png",
#           width = 10,
#           height = 8,
#           devices = NULL,
#           res = 300,
#           show = F)
plot15 <- plot_genes_jitter(mycds[hub_gene,],
                            grouping = "State", 
                            ncol= 3) +
  theme(axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))
plot15

write_fig(plot15,
          file = "subcell_hub_gene_expr_state.pdf",
          width = 15,
          height = 10,
          devices = NULL,
          res = 600,
          show = F)
write_fig(plot15,
          file = "subcell_hub_gene_expr_state.png",
          width = 15,
          height = 10,
          devices = NULL,
          res = 300,
          show = F)

plot16 <- plot_genes_in_pseudotime(mycds[hub_gene,],
                                   color_by = "State",
                                   ncol = 3)+ 
  theme(axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))
plot16
write_fig(plot16,
          file = "subcell_hub_gene_expr_state2.pdf",
          width = 15,
          height = 10,
          devices = NULL,
          res = 600,
          show = F)
write_fig(plot16,
          file = "subcell_hub_gene_expr_state2.png",
          width = 15,
          height = 10,
          devices = NULL,
          res = 300,
          show = F)

plot17 <- plot_cells(mycds,
                     genes=hub_gene,
                     label_cell_groups=FALSE,
                     show_trajectory_graph=FALSE)
save(hub_gene, mycds, file = "plot17.RData")









## devtools::install_github("dynverse/dyno")
library(dyno)
library(tidyverse)
library(Matrix)
library(Seurat)
library(dynwrap)
library(hdf5r)
#读入seurat处理后的rds文件
sdata <-TSNE
#添加raw counts和normalised expression
#seurat的矩阵需要进行行列转换，以使行为细胞，列为基因
dataset <- wrap_expression(
  counts = t(sdata@assays$RNA@counts),
  expression = t(sdata@assays$RNA@data)
)

#添加先验信息，这里添加的是开始转换的“细胞id”，后期可视化可以根据具体的轨迹推断结果进行调整
dataset$cell_ids
dataset <- add_prior_information(
  dataset,
  start_id = "GSM5708485_AAACACCAATAACTGC-1"
)

#添加数据的cluster信息，这里我们直接用“seurat_clusters”即可
dataset <- add_grouping(
  dataset,
  sdata$seurat_clusters
)
guidelines <- guidelines_shiny(dataset)
methods_selected <- guidelines$methods_selected
model_slingshot <- infer_trajectory(dataset, 'slingshot')
#model_paga_tree <- infer_trajectory(dataset, methods_selected[2])
model_paga <- infer_trajectory(dataset, methods_selected[3])
model_mst <- infer_trajectory(dataset, methods_selected[4])

# 14 风险评分与临床相关性分析----------
setwd("/data/nas1/luchunlin/project/BJTC-276")
if (! dir.exists("./13_clinical")){
  dir.create("./13_clinical")
}
setwd("./13_clinical")
## 比较高低风险组在不同临床指标下（T、N、M分期）风险评分的差异，用小提琴图展
train_phenotype<-phenotype
train_phenotype$stage<-gsub('IV','4',train_phenotype$stage)
train_phenotype$stage<-gsub('III','3',train_phenotype$stage)
train_phenotype$stage<-gsub('II','1-2',train_phenotype$stage)
train_phenotype$stage<-gsub('I','1-2',train_phenotype$stage)
table(train_phenotype$age)
### 30-50  50-70  70-90  
#train_phenotype$age<-gsub('30-40','30-50',train_phenotype$age)
#train_phenotype$age<-gsub('40-50','30-50',train_phenotype$age)
#train_phenotype$age<-gsub('50-60','50-70',train_phenotype$age)
#train_phenotype$age<-gsub('60-70','50-70',train_phenotype$age)
#train_phenotype$age<-gsub('70-80','70-90',train_phenotype$age)
#train_phenotype$age<-gsub('80-90','70-90',train_phenotype$age)

train_phenotype<-merge(train_phenotype,survival,by='sample')
train_phenotype<-train_phenotype[,-5]
colnames(train_phenotype)<-c('id','TNM stage','Age','OS','OS.time')
sub_risk <- subset(risk, select = c(id, riskScore))
train_phenotype2 <- merge(train_phenotype,
                          sub_risk,
                          by = "id")
write.table(train_phenotype2,
            file = "clinical_risk.csv",
            row.names = T,
            sep = "\t",
            quote = F)
train_phenotype2$Group<-ifelse(train_phenotype2$riskScore>median(train_phenotype2$riskScore),'High risk','Low risk')
train_phenotype2<-na.omit(train_phenotype2)
library(ggpubr)
library(Ipaper)
library(ggthemes)
table(train_phenotype$`TNM stage`)
exp_plot <- ggplot(train_phenotype2,aes(x = `TNM stage`, y = riskScore, fill = Group)) +
  #geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#CD3700","#4682B4"), name = "Group")+
  labs(title="", x="", y = "",size=20) +
  stat_compare_means(data = train_phenotype2,
                     mapping = aes(group = Group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F) +
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
        panel.grid.minor = element_blank())+facet_wrap(~`TNM stage`,scales = "free",nrow = 1) 
exp_plot

my_comparisons <- list(c("Stage 1-2","Stage 3"),c("Stage 1-2","Stage 4"),c("Stage 3","Stage 4"))
stage_data <- data.frame(riskScore = train_phenotype2$riskScore,
                         stage = factor(train_phenotype2$`TNM stage`,
                                        levels = c("Stage 1-2","Stage 3","Stage 4")))
stage_data <- na.omit(stage_data)
stage<-ggplot(stage_data,aes(x = stage, y = riskScore, fill = stage)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "riskScore",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("TNM stage") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              map_signif_level = T,
              y_position = c(3,4,3.5))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 20))+
  guides(fill='none')
stage
table(train_phenotype$Age)
my_comparisons <- list(c("30-40","40-50"),c("40-50","50-60"),c("50-60","60-70"),c('60-70','70-80'),c('70-80','80-90'))
age_data <- data.frame(riskScore = train_phenotype2$riskScore,
                         age = factor(train_phenotype2$Age,
                                        levels = c("30-40","40-50","50-60",'60-70','70-80','80-90')))
age_data <- na.omit(age_data)
age<-ggplot(age_data,aes(x = age, y = riskScore, fill = age)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "riskScore",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("TNM age") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              map_signif_level = T,
              y_position = c(2,3,2.5,3,2.5,3,3,5))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 20))+
  guides(fill='none')
age

# 15 GSVA KEGG----------
setwd("/data/nas1/luchunlin/project/BJTC-276")
if (! dir.exists("./14_GSVA")){
  dir.create("./14_GSVA")
}
setwd("./14_GSVA")
## 分别对高低风险组样本进行GSVA富集分析，以KEGG为背景基因集。
library(GSVA)
library(GSEABase)
library(limma)
## 将表达矩阵提取出来
risk2 <- risk
risk2$risk_label <- ifelse(risk$risk == 0, "High", "Low")
gsva_exp<-dat_tpm[,risk2$id]
all(colnames(gsva_exp) == risk2$id)
dim(gsva_exp)
write.table(gsva_exp,
            file = "gsva_exp.xls",
            quote = F,
            sep = "\t",
            row.names = T)
# 分组
group_risk <- risk2$risk_label %>% as.factor()
design_risk <- model.matrix(~0 + group_risk)
rownames(design_risk) <- colnames(gsva_exp)
colnames(design_risk) <- levels(group_risk)
compare_risk <- makeContrasts("High-Low", levels = design_risk)

KEGG_ref <- getGmt("/data/nas1/luchunlin/pipeline/GSVA/c2.cp.kegg.v7.4.symbols.gmt")
es_KEGG <- gsva(as.matrix(gsva_exp), KEGG_ref,
                min.sz=10, max.sz=500, verbose=TRUE)
es_KEGG <- as.data.frame(es_KEGG)
fit <- lmFit(es_KEGG, design_risk)
fit2 <- contrasts.fit(fit ,compare_risk)
fit3 <- eBayes(fit2)
allGeneSets <- topTable(fit3, coef = 1, number = Inf)

logFCcutoff <- 0
allGeneSets$change = as.factor(
  ifelse(allGeneSets$adj.P.Val < 0.05 & abs(allGeneSets$logFC) > logFCcutoff,
         ifelse(allGeneSets$logFC > logFCcutoff,'UP','DOWN'),'NOT')
)

DEGeneSets <- subset(allGeneSets,
                     allGeneSets$adj.P.Val < 0.05 & abs(allGeneSets$logFC) > 0)

write.table(allGeneSets,
            file = "GSVA_KEGG.xls",
            quote = F,
            sep = "\t",
            row.names = T)
write.table(DEGeneSets,
            file = 'diff_KEGG.xls',
            quote = F,
            sep = "\t",
            row.names = T)
## 前20个画热图
group.gsva<-risk2[,c(1,14)]
TOP_20<-DEGeneSets[c(1:20),]
heat_20<-es_KEGG[rownames(es_KEGG)%in%rownames(TOP_20),]
## 按照分组排一下序
heat_20<-t(heat_20)
heat_20<-as.data.frame(heat_20)
group.gsva<-group.gsva[order(group.gsva$risk_label),]
heat_20<-heat_20[group.gsva$id,]
heat_20<-t(heat_20)
colnames(group.gsva)<-c('sample','group')
annotation_col<-as.data.frame(group.gsva$group)
colnames(annotation_col)='Group'
rownames(annotation_col)=colnames(heat_20)
color.key<-c("#3300CC", "#3399FF", "white", "#FF3333", "#CC0000")
ann_colors<-list(
  Group = c(Low="#00FFFF",High="#FFAEB9"))

pheatmap(heat_20,
         color = colorRampPalette(color.key)(50),
         border_color = NA,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         labels_row = NULL,
         clustering_method = 'ward.D2',
         show_rownames = T,
         show_colnames = F,
         fontsize_col = 5,
         cluster_cols = F,
         cluster_rows = T,
         width = 10,
         main = 'GSVA KEGG')

# 16 TMB评分---------
setwd("/data/nas1/luchunlin/project/BJTC-276")
if (! dir.exists("./15_TMB")){
  dir.create("./15_TMB")
}
setwd("./15_TMB")
## 高低风险组患者TMB评分箱线图
library(TCGAmutations)
tmp<-as.data.frame(tcga_available())
dt<-TCGAmutations::tcga_load(study = 'OV')
dt<-dt@data
dt1<-as.data.frame(table(dt$Tumor_Sample_Barcode))
names(dt1)<-c('Barcode','Freq')
dt1$tmb<-dt1$Freq/38
tmb<-dt1
tmb<-separate(tmb,'Barcode',into = 'sample',sep = '-01W')
tmb<-separate(tmb,'sample',into = 'sample',sep = '-01D')
tmb<-separate(tmb,'sample',into = 'sample',sep = '-02D')
tmb<-separate(tmb,'sample',into = 'sample',sep = '-03D')
tmb<-separate(tmb,'sample',into = 'sample',sep = '-11D')
tmb<-separate(tmb,'sample',into = 'sample',sep = '-12D')
tmb<-separate(tmb,'sample',into = 'sample',sep = '-02W')
tmb<-separate(tmb,'sample',into = 'sample',sep = '-03W')
tmb<-tmb[tmb$sample%in%group.gsva$sample,]
tmb<-merge(tmb,group.gsva,by='sample')
tmb<-tmb[-8,]
## 38是人类基因外显子的长度
my_comparisons <- list(c("High","Low"))
tmb.bar<-ggplot(tmb,aes(x = group, y = tmb, fill = group)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "TMB",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("TMB score") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              map_signif_level = T,
              y_position = c(10))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 20))+
  guides(fill='none')
tmb.bar

# 17 DNA 干性评分 (DNAs) 和 RNA 干性评分 (RNAs) ---------
setwd("/data/nas1/luchunlin/project/BJTC-276")
if (! dir.exists("./16_DNAss_RNAss")){
  dir.create("./16_DNAss_RNAss")
}
setwd("./16_DNAss_RNAss")
##分析DNA 干性评分 (DNAs) 和 RNA 干性评分 (RNAs) 与风险评分的相关性。
## 17-1 RNA干性评分-----
###样本fpkm表达量
exp.ss<-log2(dat_fpkm+1)
exp.ss<-exp.ss[,group.gsva$sample]
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

## 17-2 DNA干性评分-----

#gunzip("TCGA-OV.methylation27.tsv.gz", remove = FALSE)
met<-read_delim('TCGA-OV.methylation27.tsv',delim="\t",col_names = TRUE)  ##甲基化数据
met<-data.frame(met)
colnames(met)
met<-column_to_rownames(met,var = 'Composite.Element.REF')

##下载甲基化数据
# coad.samples <- matchedMetExp("TCGA-THCA", n = 10)
# query <- GDCquery(
#   project = "TCGA-THCA",
#   data.category = "DNA methylation",
#   platform = "Illumina Human Methylation 450",
#   legacy = TRUE,
#  # barcode = coad.samples
# )
# GDCdownload(query)  ##运行速度比较慢
# met <- GDCprepare(query, save = FALSE)

#save(met, file = "met.RData")

#load('met.RData')

data.met <- na.omit(met)

load('Stemness_index.rda')
names(mDNAsi$Weight)<-mDNAsi$`Probe ID`
common
common <- intersect(names(mDNAsi$Weight), rownames(data.met))
X <- data.met[common, ]
w <- mDNAsi$Weight[common]
w_t<-as.numeric(t(w))
score <-w_t %*% as.matrix(X)  ##%*%矩阵乘法
score <- score - min(score)
score <- score / max(score)
score <-data.frame(t(score ))
colnames(score)<-'mDNAsi'
write.table(score,'mDNAsi.txt',sep='\t',quote=F)
mDNA<-read.table('mDNAsi.txt',header=T)
mRNA<-read.table('mRNAsi.txt',header=T)
mDNA$sample<-NA
mDNA$sample<-rownames(mDNA)
mRNA$sample<-NA
mRNA$sample<-rownames(mRNA)
mDNA$sample<-gsub('.','-',mDNA$sample,fixed = T)
m<-merge(mDNA,mRNA,by='sample')
m<-na.omit(m)
rownames(m)<-m$sample
m<-m[,-1]
exp_m<-subset(exp.ss,select=rownames(m))

riskscore<-data.frame(sample = risk$id,
                      riskScore=risk$riskScore)


riskscore<-riskscore[riskscore$sample%in%rownames(m),]

rownames(riskscore)<-riskscore$sample
riskscore<-riskscore[rownames(m),]
riskscore<-data.frame(row.names = riskscore$sample,
                      riskScore=riskscore$riskScore)
riskscore<-t(riskscore)
group=as.vector(ifelse(as.numeric( riskscore)> median(as.numeric(riskscore)),'High','Low'))  
train<-data.frame(row.names=colnames(riskscore), riskscore=as.numeric(riskscore),mRNAsi=m$mRNAsi,mDNAsi=m$mDNAsi,Group=group)
train<-train[order(train$riskscore, decreasing = F),] 

##散点图线性拟合
p1<-ggplot(train,aes(x= riskscore,y=mRNAsi))+geom_point(aes( color = Group), size=5)+geom_smooth(method = 'lm', formula = y ~ x, se = T,color='red')+ stat_cor(data=train, method = "spearman")+
  scale_color_manual(values = c("#2774C4", "#BEBEBE"))+
  scale_fill_manual(values = c("#2774C4", "#BEBEBE")) +theme_bw()+
  theme(axis.title.x =element_text(size=20,family = "Times", face = "bold"),
        axis.text.x =element_text(size=16,family = "Times", face = "bold"),
        axis.title.y =element_text(size=20,family = "Times", face = "bold"),
        axis.text.y=element_text(size=16,family = "Times", face = "bold"),
        legend.title=element_text(size=15,family = "Times", face = "bold") , 
        legend.text=element_text(size=14,family = "Times", face = "bold"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(x="riskScore",y="mRNAsi")+
  theme(legend.position = "bottom")  ##修改图例位置
##箱线图显著性检验
my_comparisons = list( c("High", "Low"))
library(ggsci)
p2<-ggplot(train, aes(x=Group, y=mRNAsi,fill=Group)) + 
  geom_boxplot()+
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     method = "wilcox.test",cex=6)+
  scale_fill_npg()+
  scale_fill_manual(values=c("#2774C4", "#BEBEBE")) + 
  labs(x = "", y = "", title = "") + 
  #theme_bw() + 
  # geom_text(aes(label = B, vjust = 1.1, hjust = -0.5, angle = 45), show_guide = FALSE) + 
  theme(panel.grid =element_blank()) +   ## 删去网格线
  theme(axis.text = element_blank()) +   ## 删去所有刻度标签
  theme(axis.ticks = element_blank())+
  theme(panel.background = element_blank())+#去除背景
  guides(fill='none')  ##去除图例

##图片结合
library(patchwork)  
##ggarrange(p1,p2,ncol = 2,nrow =1,widths = c(2,0.3),heights = c(1,1))  ##y轴没对应
P1<-p1+p2+plot_layout(widths = c(2, 0.3))
P1
ggsave('mRNAsi.boxplot.pdf',width=8,height=7)
ggsave('mRNAsi.boxplot.png',width=8,height=7)


p1<-ggplot(train,aes(x= riskscore,y=mDNAsi))+geom_point(aes( color = Group), size=5)+geom_smooth(method = 'lm', formula = y ~ x, se = T,color='red')+ stat_cor(data=train, method = "spearman")+
  scale_color_manual(values = c("#2774C4", "#BEBEBE"))+
  scale_fill_manual(values = c("#2774C4", "#BEBEBE")) +theme_bw()+
  theme(axis.title.x =element_text(size=20,family = "Times", face = "bold"),
        axis.text.x =element_text(size=16,family = "Times", face = "bold"),
        axis.title.y =element_text(size=20,family = "Times", face = "bold"),
        axis.text.y=element_text(size=16,family = "Times", face = "bold"),
        legend.title=element_text(size=15,family = "Times", face = "bold") , 
        legend.text=element_text(size=14,family = "Times", face = "bold"))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+
  labs(x="riskScore",y="mDNAsi")+
  theme(legend.position = "bottom")  ##修改图例位置
##箱线图显著性检验
my_comparisons = list( c("High", "Low"))
library(ggsci)
p2<-ggplot(train, aes(x=Group, y=mDNAsi,fill=Group)) + 
  geom_boxplot()+
  stat_compare_means(comparisons = my_comparisons,
                     label = "p.signif",
                     method = "wilcox.test",cex=6)+
  scale_fill_npg()+
  scale_fill_manual(values=c("#2774C4", "#BEBEBE")) + 
  labs(x = "", y = "", title = "") + 
  #theme_bw() + 
  # geom_text(aes(label = B, vjust = 1.1, hjust = -0.5, angle = 45), show_guide = FALSE) + 
  theme(panel.grid =element_blank()) +   ## 删去网格线
  theme(axis.text = element_blank()) +   ## 删去所有刻度标签
  theme(axis.ticks = element_blank())+
  theme(panel.background = element_blank())+#去除背景
  guides(fill='none')  ##去除图例

##图片结合
library(patchwork)  
##ggarrange(p1,p2,ncol = 2,nrow =1,widths = c(2,0.3),heights = c(1,1))  ##y轴没对应
P2<-p1+p2+plot_layout(widths = c(2, 0.3))
P2
ggsave('mDNAsi.boxplot.pdf',width=8,height=7)
ggsave('mDNAsi.boxplot.png',width=8,height=7)


# 18 风险组和不同cluster相关性---------
setwd("/data/nas1/luchunlin/project/BJTC-276")
if (! dir.exists("./17_cor")){
  dir.create("./17_cor")
}
setwd("./17_cor")
## 风险组与不同簇与生存状态之间的相关性。（桑基图）
## 需要文件：1.不同簇和生存状态对应关系。2.高低风险组和生存状态对应关系。
cluster.os<-data.frame(sample=rownames(cluster),
                       cluster=cluster$cluster)
cluster.os<-merge(cluster.os,survival,by='sample')
cluster.os<-cluster.os[,c(1:3)]
cluster.os$OS<-ifelse(cluster.os$OS==1,'Dead','Alive')
cluster.os<-column_to_rownames(cluster.os,var = 'sample')
colnames(cluster.os)<-c('Cluster','Status')
risk.os<-risk[,c(1,2,13)]
risk.os$OS<-ifelse(risk.os$OS==1,'Dead','Alive')
risk.os$risk<-ifelse(risk.os$risk==1,'Low','High')
risk.os<-column_to_rownames(risk.os,var = 'id')
colnames(risk.os)<-c('Status','Risk group')
cluster.os.risk<-merge(cluster.os,risk.os,by='Status')
cluster.os.risk$link<-1
cluster.os.risk<-reshape::melt(cluster.os.risk,id='link')
variable<-summary(cluster.os.risk$variable)
cluster.os.risk$flow<-rep(1:variable[1],length(variable))
library(ggalluvial)
mycol <- c('#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462',
           '#B3DE69', '#FCCDE5', '#BC80BD', '#CCEBC5', '#FFED6F', '#E41A1C', '#377EB8',
           '#4DAF4A', '#984EA3', '#FF7F00', '#FFFF33', '#A65628', '#F781BF', '#66C2A5', 
           '#6181BD', '#F34800', '#64A10E', '#FF00FF', '#c7475b', '#049a0b', '#BEAED4', 
           '#FDC086', '#FFFF99', '#386CB0', '#F0027F', '#4253ff', '#ff4308', '#D8D155',
           '#64495D', '#7CC767')
p <- ggplot(cluster.os.risk, aes(x = variable, y = link,
                       stratum = value, alluvium = flow, fill = value)) +
  geom_stratum() +  #冲击图中的堆叠柱形图
  geom_flow(aes.flow = 'forward') +  #冲击图连线绘制
  scale_fill_manual(values = mycol) +  #颜色赋值
  geom_text(stat = 'stratum', infer.label = TRUE, size = 3) +  #添加 标签
  scale_x_discrete(limits = c('Cluster', 'Status', 'Risk group')) +  #定义 列的展示顺序
  labs(x = '', y = '') +  #去除 x 轴和 y 轴标题
  theme(legend.position = 'none', panel.background = element_blank(),
        line = element_blank(), axis.text.y = element_blank())  #去除背景和图例

p

# 19 ssGSEA---------
setwd("/data/nas1/luchunlin/project/BJTC-276")
if (! dir.exists("./18_ssGSEA")){
  dir.create("./18_ssGSEA")
}
setwd("./18_ssGSEA")
## 高低风险组的免疫浸润细胞进行比较，并分析风险评分与免疫浸润细胞的相关性。

## Spearman相关分析
##19-1 ssGSEA--------
library(GSVA)
gene_set <- read.table("/data/nas1/luchunlin/pipeline/ssGSEA/mmc3.txt",
                       header = T,
                       sep ="\t")
dat_final <- as.matrix(dat_tcga)
dat_final<-dat_final[,colnames(dat_final)%in%rownames(train_data)]
gene_list <- split(as.matrix(gene_set)[,1],
                   gene_set[,2])

ssgsea_score = gsva(dat_final, gene_list, 
                    method = "ssgsea", 
                    ssgsea.norm = TRUE, 
                    verbose = TRUE)
write.table(ssgsea_score,
            file = "ssgsea_result.xls",
            sep = "\t",
            quote = F)
## 19-2差异免疫细胞鉴定-------
tiics_result <- ssgsea_score
pvalue = padj = log2FoldChange <- matrix(0, nrow(tiics_result), 1)
group_High<-as.character(risk$id[which(risk$risk==0)])
group_Low<-as.character(risk$id[which(risk$risk==1)])

for (i in 1:nrow(tiics_result)){
  pvalue[i, 1] = p.value = wilcox.test(tiics_result[i, group_High],
                                       tiics_result[i, group_Low])$p.value
  log2FoldChange[i, 1] = mean(tiics_result[i, group_High]) - 
    mean(tiics_result[i, group_Low])
}
padj <- p.adjust(as.vector(pvalue), "fdr", n = length(pvalue))
rTable <- data.frame(log2FoldChange, 
                     pvalue, 
                     padj,
                     row.names = rownames(tiics_result))
Low <- signif(apply(tiics_result[rownames(rTable), group_Low], 
                        1,
                        mean), 4)
High <- signif(apply(tiics_result[rownames(rTable), group_High], 
                     1, 
                     mean), 4)
rTable <- data.frame(Low, 
                     High,
                     rTable[, c("padj", "pvalue", "log2FoldChange")])
rTable$immune_cell <- rownames(rTable)
rTable$sig <- ifelse(rTable$padj < 0.05,
                     ifelse(rTable$padj < 0.01, 
                            ifelse(rTable$padj < 0.001,
                                   ifelse(rTable$padj < 0.0001,
                                          paste(rTable$immune_cell, "****",  sep = ""),
                                          paste(rTable$immune_cell, "***", sep = "")),
                                   paste(rTable$immune_cell, "**", sep = "")),
                            paste(rTable$immune_cell, "*",  sep = "")), 
                     rTable$immune_cell)

diff_Table<-rTable[which(rTable$padj<0.05),]
##19
write.table(diff_Table,file = 'diff_tiics.xls',
            sep = '\t',
            row.names = T)
### 箱线图----
library(tidyr)
library(ggplot2)
library(ggpubr)
library(Ipaper)

xCell2 <- data.frame(Immune_Cell=rownames(tiics_result), 
                     tiics_result, 
                     pvalue=rTable$pvalue)
#xCell3 <- xCell2[which(xCell2$pvalue<0.05),]
xCell3<-xCell2
diff_tiics <- rownames(xCell3)
violin_dat <- gather(xCell3, key=Group, value=score, -c("Immune_Cell","pvalue"))

violin_dat$Group <- ifelse(gsub("\\.","-",violin_dat$Group) %in% group_Low,
                           "Low", "High") 
head(violin_dat)
boxplot_diff_TIICs <- ggplot(violin_dat, aes(x=Immune_Cell, 
                                             y=score,
                                             fill=Group)) +
  # geom_violin(trim=T,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)#"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar",
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  geom_point(aes(fill = Group),
             size = 0.05,
             position = position_dodge(0.9))+
  scale_fill_manual(values= c("#48D1CC", "#FF7256"))+ #设置填充的颜色
  labs(title="", x="", y = "Score",size=20) +
  stat_compare_means(data = violin_dat,
                     mapping = aes(group = Group),
                     label ="p.signif",
                     hide.ns = F,
                     method = 'wilcox') +
  theme_bw()+#把背景设置为白底
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18), # 将图表标题居中
        axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=10), #设置x轴刻度标签的字体显示倾斜角度为45度，并向下调整1(hjust = 1)，字体大小为14
        axis.text.y=element_text(hjust=0.5,colour="black",size=12), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.x=element_text(size=16,face="bold"),#设置x轴标题的字体属性
        axis.title.y=element_text(size=14,face="bold"), #设置y轴标题的字体属性
        legend.text=element_text(face="bold", hjust = 0.5,colour="black", size=11), #设置图例的子标题的字体属性
        legend.title=element_text(face="bold", colour="black", size=11),#设置图例的总标题的字体属性
        #legend.justification=c(-0.1,1.2), #可调整图例的位置。##(1,1)第一个1是调整图例在图的内外(左右移动)，第二个1是在图中上下移动。
        #legend.position=c(0, 1.04), #legend.position=c(0,1)左上角，(1,1)是在右上角。
        panel.grid.major = element_blank(), #不显示网格线
        panel.grid.minor = element_blank()) #不显示网格线
boxplot_diff_TIICs
## 风险评分与免疫浸润细胞相关性------
## 风险评分数据 和 免疫细胞浸润矩阵 (tiics_result)行名均为样本名称
risk.score<-data.frame(row.names = risk$id,riskScore=risk$riskScore)
risk.score<-as.matrix(risk.score)
tiics_exp<-t(tiics_result)
risk.score<-risk.score[rownames(tiics_exp),]
### 批量操作的具体实现过程：
### 1.设定容器,最终生成的数据放在什么地方？
correlation <- data.frame()
### 2.批量把数据导出到容器
for(i in 1:length(colnames(tiics_exp))){
  ## 1.指示
  print(i)
  ## 2.计算
  dd = cor.test(as.numeric(tiics_exp[,i]),risk.score,method="spearman")
  ## 3.填充
  correlation[i,1] = colnames(tiics_exp)[i]
  correlation[i,2] = dd$estimate
  correlation[i,3] = dd$p.value
}
### 修改名称
colnames(correlation) <- c("cell","cor","p.value")
#correlation_EPAS1<-correlation_EPAS1[order(correlation_EPAS1$cor,decreasing = T),]
correlation$correlation<-cut(abs(correlation$cor),breaks = c(0,0.1,0.2,0.3,0.4,0.5),labels = c("0.1","0.2","0.3","0.4","0.5"))
### 棒棒糖图
library(ggpubr)
#display.brewer.all()

ggdotchart(correlation,x='cell',y='cor',
           size = 'correlation',
           sorting='descending',
           rotate = T,
           color = 'correlation',
           palette =brewer.pal(5,"GnBu"), 
           add = 'segment',                     # 添加棒子
           add.params = list(color='black',size=0.5),
           ggtheme = theme_bw(base_size = 13)+
             theme(legend.text =element_text(size = 10)),# 改变主题
           xlab="",
           ylab = 'correlation',
           title = 'correlation')

# 20 预后模型--------
setwd("/data/nas1/luchunlin/project/BJTC-276")
if (! dir.exists("./19_progmod")){
  dir.create("./19_progmod")
}
setwd("./19_progmod")
## 20-1 单因素cox--------
train_phenotype3<-data.frame(sample=clinical$submitter_id.samples,
                             stage=clinical$clinical_stage,
                             age=clinical$age_at_index.demographic)
train_phenotype3<-merge(train_phenotype3,survival,by='sample')
train_phenotype3<-train_phenotype3[,-5]

train_phenotype3$stage<-gsub('A','',train_phenotype3$stage)
train_phenotype3$stage<-gsub('B','',train_phenotype3$stage)
train_phenotype3$stage<-gsub('C','',train_phenotype3$stage)

train_phenotype3$stage<-gsub('Stage IV','4',train_phenotype3$stage)
train_phenotype3$stage<-gsub('Stage III','3',train_phenotype3$stage)
train_phenotype3$stage<-gsub('Stage II','2',train_phenotype3$stage)
train_phenotype3$stage<-gsub('Stage I','1',train_phenotype3$stage)

colnames(train_phenotype3)<-c('id','stage','age','OS','OS.time')

train_risk_clinical <- merge(train_phenotype3,
                             sub_risk,
                             by = "id")
rownames(train_risk_clinical) <- train_risk_clinical$id
train_risk_clinical = subset(train_risk_clinical, select = -c(id))
dim(train_risk_clinical)
colnames_train <- colnames(train_risk_clinical)
covariates_train <- colnames_train[-which(colnames_train %in% c("OS", "OS.time"))]

train_risk_clinical$stage<-factor(train_risk_clinical$stage)
library(survival)
res.risk = coxph(Surv(time = OS.time, event = OS) ~ riskScore, data = train_risk_clinical) %>% summary
res.risk = c(res.risk$conf.int[-2], res.risk$coefficients[5])
res.age = coxph(Surv(time = OS.time, event = OS) ~ age, data = train_risk_clinical) %>% summary
res.age = c(res.age$conf.int[-2], res.age$coefficients[5])

res.stage = coxph(Surv(time = OS.time, event = OS) ~ stage, data = train_risk_clinical) %>% summary
res.stage = cbind(res.stage$conf.int[,-2], res.stage$coefficients[,5])

res.ref = c(1,1,1,NA)
res = rbind(res.risk, res.age,res.ref,res.stage) %>% as.data.frame()
rownames(res)
res$Indicators = c("riskScore","Age","Stage 1 (Reference)","Stage 2","Stage 3","Stage 4")
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
hz[3] <- ""

tabletext <- cbind(c(NA,rownames(res2)),
                   c("P value",ifelse(res2$p.value<0.0001,
                                      "< 0.0001",
                                      round(res2$p.value,4))),
                   c("Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1
library(forestplot)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,rep(FALSE, 7)),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,res2$HR),
           lower=c(NA,res2$HR.95L), #95%置信区间下限
           upper=c(NA,res2$HR.95H), #95%置信区间上限
           boxsize=0.1,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(-1,0,1,2,3,4,5,6), #横坐标刻度
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

## 20-2 多因素cox--------
res.mul = coxph(Surv(time = OS.time, event = OS) ~ riskScore + age+stage, data = train_risk_clinical)%>% summary
res.mul = cbind(res.mul$conf.int[,-2], res.mul$coefficients[,5]) %>% as.data.frame()
rownames(res.mul)
res.mul = rbind(res.mul[1:2,], res.ref, res.mul[3:5,])
res.mul$Indicators = c("riskScore","Age","Stage 1 (Reference)","Stage 2","Stage 3","Stage 4")
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
hz[3] <- ""

tabletext <- cbind(c(NA,rownames(multi_res)),
                   c("P value",ifelse(multi_res$p.value<0.0001,
                                      "< 0.0001",
                                      round(multi_res$p.value,4))),
                   c("Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,multi_res$HR),
           lower=c(NA,multi_res$HR.95L), #95%置信区间下限
           upper=c(NA,multi_res$HR.95H), #95%置信区间上限
           boxsize=0.1,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(-1,0,1,2,3,4,5,6), #横坐标刻度
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
## 20-3 构建COX模型，绘制列线图------
multi_cov<-c('riskScore','age','stage')
cox_data_prog <- as.formula(paste0('Surv(OS.time, OS)~',
                                   paste(multi_cov,
                                         sep = '',
                                         collapse = '+')))

# Nomogram
library(rms)
ddist <- datadist(train_risk_clinical)
options(datadist='ddist')

# 构建COX模型，绘制列线图

res.cox <- psm(cox_data_prog,
               data = train_risk_clinical, dist = 'lognormal')
surv <- Survival(res.cox) # 构建生存概率函数
function(x) surv(1825, x) # 5年事件发生概率
function(x) surv(2555, x) #7年事件发生概率
function(x) surv(3285, x) # 9年事件发生概率

nom.cox <- nomogram(res.cox,
                    fun = list(function(x) surv(1825, x),
                               function(x) surv(2555, x),
                               function(x) surv(3285, x)),
                    funlabel=c("5-year Survival Probability", "7-year Survival Probability", "9-year Survival Probability"),
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

## 20-4 构建校准曲线--------
##绘制1年生存期校准曲线
coxm_5 <- cph(cox_data_prog,
              data=train_risk_clinical,
              surv=T,
              x=T,
              y=T,
              time.inc = 5*365)
cal_5<-calibrate(coxm_5,u=5*365,cmethod='KM',m=100,B=1000)

##绘制3年生存期校曲线
##time.in 和 u 要是一样的，都是要评价的时间节点
coxm_7 <- cph(cox_data_prog,
              data=train_risk_clinical,
              surv=T,
              x=T,y=T,
              time.inc = 7*365)
cal_7 <-calibrate(coxm_7,u=7*365,cmethod='KM',m=100,B=1000)

coxm_9 <- cph(cox_data_prog,
              data=train_risk_clinical,
              surv=T,
              x=T,y=T,
              time.inc = 9*365)
cal_9 <-calibrate(coxm_9,u=9*365,cmethod='KM',m=100)

par(mar=c(7,4,4,3),cex=1.5)
plot(cal_5,
     subtitles = F,
     lwd=2,lty=1, ##设置线条形状和尺寸
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 5-9 year Progression-free Interval',#便签
     ylab='Actual 5-9 year Progression-free Interval(Proportion)',#标签
     col="#00468b",#设置一个颜色
     xlim = c(0,0.8),ylim = c(0,1)) ##x轴和y轴范围
plot(cal_7,
     add = T,
     subtitles = F,
     lwd=2,lty=1,  ##设置线条宽度和线条类型
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 5-9 year Progression-free Interval',#便签
     ylab='Actual 5-9 year Progression-free Interval(Proportion)',#标签
     col="#ed0000",#设置一个颜色
     xlim = c(0,0.8),ylim = c(0,1)) ##x轴和y轴范围
plot(cal_9,
     add = T,
     subtitles = F,
     lwd=2,lty=1, ##设置线条形状和尺寸
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 5-9 year Progression-free Interval',#便签
     ylab='Actual 5-9 year Progression-free Interval(Proportion)',#标签
     col="#42b540",#设置一个颜色
     xlim = c(0,0.8),ylim = c(0,1)) ##x轴和y轴范围



#加上图例
legend("bottomright", legend=c("5-year", "7-year", "9-year"), 
       col=c("#00468b", "#ed0000", "#42b540"), 
       lwd=2)
#调整对角线
abline(0,1,lty=5,lwd=2,col="grey")
## 20-5 KM曲线--------

res.mul = coxph(Surv(time = OS.time, event = OS)~ riskScore + age+stage, data = train_risk_clinical)
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
## 20-6 ROC曲线-------
# 开始验证
train_risk_clinical2 <- train_risk_clinical
library(survival)
library(survminer)

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

for_multi_ROC <- multi_ROC_out(time_vector = c(365*seq(5,9,2)), 
                               risk_score_table = train_risk_clinical2)
for_multi_ROC$Time <- factor(for_multi_ROC$Time)

# 画ROC曲线 
library(scales)
library(geomROC)
library(plotROC)
auc_y7 <- round(for_multi_ROC[which(for_multi_ROC$Time==2555),5][1],2)
auc_y9 <- round(for_multi_ROC[which(for_multi_ROC$Time==3285),5][1],2)
auc_y5 <- round(for_multi_ROC[which(for_multi_ROC$Time==1825),5][1],2)

ROC <- ggplot(for_multi_ROC, aes(x=False_positive,
                                 y=True_positive, 
                                 label=Cut_values, 
                                 color=Time)) + 
  scale_color_manual(breaks = c("1825", "2555", "3285"),
                     labels = c("5 years", "7 years", "9 years"),
                     values = c("#4682B4", "#FF4040", "#20B2AA")) +
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
           label = c(paste('AUC of 5 years =', format(auc_y5,nsmall=2)),
                     paste('AUC of 7 years =', format(auc_y7,nsmall=2)),
                     paste('AUC of 9 years =', format(auc_y9,nsmall=2))))
ROC
ggsave('ROC for Prognosis Model.png', ROC,width = 5, height = 4)
ggsave('ROC for Prognosis Model.pdf', ROC,width = 5, height = 4)
