## rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/BJTC-290")
if (! dir.exists("./00_rawdata")){
  dir.create("./00_rawdata")
}
setwd("./00_rawdata")
library(GEOquery)
library(Biobase)
library(tidyr)
library(GEOquery)
library(Biobase)
## 01-1 GSE93971--------
gset<-getGEO("GSE93971",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
a=gset[[1]]
pd<-pData(a)
library(AnnoProbe)
# GPL21827[Accession] - GEO DataSets Result - NCBI - NIH
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL21827
gpl='GPL21827'
probe2symbol=idmap(gpl,type = 'pipe')
head(probe2symbol)
colnames(probe2symbol)<-c('ID','symbol')
dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symbol$ID<-as.character(probe2symbol$ID)
library(tidyverse)
dat<-dat %>%
  inner_join(probe2symbol,by='ID')%>% 
  select(-ID)%>%     ## 去除多余信息
  select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
## 01-2 GSE93973--------
gset2<-getGEO("GSE93973",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr2<-as.data.frame(exprs(gset2[[1]]))
a2=gset2[[1]]
pd2<-pData(a2)
gpl2<-getGEO("GPL18058",destdir = '.')
gpl2<-Table(gpl2)
a2=gset2[[1]]
colnames(gpl2)
probe2symbol2<-gpl2 %>%
  select('ID','miRNA_ID')%>%
  filter('miRNA_ID'!='')%>%
  separate('miRNA_ID',c('symbol','drop'),sep = '/')%>%
  select(-drop)
probe2symbol2=probe2symbol2[probe2symbol2$symbol!='',]
dat2<-expr2
dat2$ID<-rownames(dat2)
dat2$ID<-as.character(dat2$ID)
probe2symbol2$ID<-as.character(probe2symbol2$ID)
dat2<-dat2 %>%
  inner_join(probe2symbol2,by='ID')%>% 
  select(-ID)%>%     ## 去除多余信息
  select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除

dim(dat2)
dat2[is.na(dat2)] <- 0
## 1963   14
# 02 差异分析---------
## limma 
setwd("/data/nas1/luchunlin/project/BJTC-290")
if (! dir.exists("./01_DEGs")){
  dir.create("./01_DEGs")
}
setwd("./01_DEGs")
## 02-1 GSE93971差异分析----

## 分组矩阵

group<-data.frame(sample=colnames(dat),
                  group=c(rep('β_thalassemia_HbF',7),rep('control',7)))
type<-group[,2]
design <- model.matrix(~ -1+factor(type,levels=c('control','β_thalassemia_HbF'))) 
colnames(design)<-c('control','β_thalassemia_HbF')
rownames(design)<-group$sample
library(limma)
# 对每一个基因进行线性模型构建
fit=lmFit(dat,design)
# 构建比较矩阵
contrast.matrix=makeContrasts(ControlVSMG=β_thalassemia_HbF-control,levels = design)
# 构建芯片数据线性模型，计算估计的相关系数和标准差
fit2=contrasts.fit(fit,contrast.matrix)
# 基于贝叶斯计算T值，F值和log-odds
fit2<-eBayes(fit2,0.01)
tempOutput = topTable(fit2, coef=1, n=Inf)
DEG= na.omit(tempOutput)

# 筛选差异基因
logFC_cutoff <- 0
DEG$change = as.factor(
  ifelse(DEG$adj.P.Val <0.05 & abs(DEG$logFC) > logFC_cutoff,
         ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
sig_diff <- subset(DEG,
                   DEG$adj.P.Val < 0.05 & abs(DEG$logFC) > logFC_cutoff)

dim(DEG)
dim(sig_diff)
## 622
summary(sig_diff$change)
# DOWN  NOT   UP 
# 149        473
write.table(group,file = "group.xls",
            quote = F,
            sep = "\t",
            row.names = T)
write.table(DEG,file = "DEG_all.xls",
            quote = F,
            sep = "\t",
            row.names = T)
write.table(sig_diff,file = "DEG_sig.xls",
            quote = F,
            sep = "\t",
            row.names = T)
sig_diff_gene<-as.data.frame(rownames(sig_diff))
## 火山图
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(Ipaper)
library(scales)
library(ggrepel)

dat_rep<-DEG[rownames(DEG)%in%
               c(head(rownames(subset(sig_diff,sig_diff$logFC>5.91)),10),
                 head(rownames(subset(sig_diff,sig_diff$logFC< -1.5)),13)),]
dat_rep<-dat_rep[-c(12,13),]
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
    size = 3,
    box.padding = unit(0.5, "lines"),
    min.segment.length = 0,
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )+
  labs(x = "log (Fold Change)",
       y = "-log10 (adj.P.Val)")
volcano_plot
ggsave('volcano.png', volcano_plot,width = 8, height = 7)
ggsave('volcano.pdf', volcano_plot,width = 8, height = 7)
## 热图
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(gplots)
group_rt<-group$group%>%as.data.frame()
rt<-dat
colnames(group_rt)<-'group'
rownames(group_rt)<-group$sample

heat<-rt[rownames(dat_rep),]
#x<-log2(heat+1)
x<-t(scale(t(heat)))
ann_colors<-list(
  group = c(control="#00CED1",β_thalassemia_HbF="#F08080"),
  Change=c(Up="#FF0000",Down="#436EEE"))
annotation_raw<-data.frame(row.names = rownames(heat),
                           Change=factor(rep(c('Up','Down'),c(10,10))))

pheatmap(mat=x,
         annotation_col = group_rt,
         annotation_row = annotation_raw,
         color=bluered(100),
         scale = "row",
         annotation_colors = ann_colors,
         fontsize = 10,
         show_colnames = F,
         cluster_cols = F,
         cluster_rows = T,
         show_rownames = T,
         annotation_names_row = F)
## 02-1 GSE93973差异分析----

## 分组矩阵

group2<-data.frame(sample=colnames(dat2),
                  group=c(rep('β_thalassemia_HbF',7),rep('control',7)))
type2<-group2[,2]
design2 <- model.matrix(~ -1+factor(type2,levels=c('control','β_thalassemia_HbF'))) 
colnames(design2)<-c('control','β_thalassemia_HbF')
rownames(design2)<-group2$sample
library(limma)

# 对每一个基因进行线性模型构建
fit2=lmFit(dat2,design2)
# 构建比较矩阵
contrast.matrix2=makeContrasts(ControlVSMG=β_thalassemia_HbF-control,levels = design2)
# 构建芯片数据线性模型，计算估计的相关系数和标准差
fit2_2=contrasts.fit(fit2,contrast.matrix2)
# 基于贝叶斯计算T值，F值和log-odds
fit2_2<-eBayes(fit2_2,0.01)
tempOutput2 = topTable(fit2_2, coef=1, n=Inf)
DEG_mir= na.omit(tempOutput2)
# 筛选差异基因
logFC_cutoff <- 0
DEG_mir$change = as.factor(
  ifelse(DEG_mir$P.Value & abs(DEG_mir$logFC) > logFC_cutoff,
         ifelse(DEG_mir$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
diff_mir <- subset(DEG_mir,
                   DEG_mir$P.Value < 0.05 & abs(DEG_mir$logFC) > logFC_cutoff)

dim(DEG_mir)
dim(diff_mir)
## 
summary(diff_mir$change)

# 03 WGCNA------
setwd("/data/nas1/luchunlin/project/BJTC-290")
if (! dir.exists("./02_WGCNA")){
  dir.create("./02_WGCNA")
}
setwd("./02_WGCNA")
library(WGCNA)
library(reshape2)
library(stringr)
options(stringsAsFactors = F)
enableWGCNAThreads()
dat_final<-dat
exprMat<-log2(dat_final)
dim(exprMat)
type = "unsigned"
corType = "pearson"
corFnc = ifelse(corType=="pearson", cor, bicor)
# 对二元变量，如样本性状信息计算相关性时，
# 或基因表达严重依赖于疾病状态时，需设置下面参数
maxPOutliers = ifelse(corType=="pearson",1,0.05)
# 关联样品性状的二元变量时，设置
robustY = ifelse(corType=="pearson",T,F)
dataExpr <- exprMat[rownames(dat_final),]
## 03-1 数据筛选-----
## 筛选中位绝对偏差（MAD）前75%的基因，至少MAD大于0.01
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
# [1] 14 26380

## 03-2 软阈值筛选----
## 样本聚类，查看是否有利群样本
sampleTree = hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
#abline(h = 200, col = "red")

# Determine cluster under the line
# 剪枝算法，cutHeight修剪树枝的高度，minSize集群最小数
#clust = cutreeStatic(sampleTree, cutHeight = 200, minSize = 10)
#table(clust)
# clust 1 contains the samples we want to keep.
#keepSamples = (clust==1)
## 符合要求的数据
#dataExpr = dataExpr[keepSamples,]
## 提取列
nGenes = ncol(dataExpr)
## 提取行
nSamples = nrow(dataExpr)
# 设置网络构建参数选择范围，计算无尺度分布拓扑矩阵
powers = c(c(1:10), seq(from = 12, to=30, by=2))
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
# 经验power (无满足条件的power时选用)
# 无向网络在power小于15或有向网络power小于30内，没有一个power值可以使
# 无标度网络图谱结构R^2达到0.8，平均连接度较高如在100以上，可能是由于
# 部分样品与其他样品差别太大。这可能由批次效应、样品异质性或实验条件对
# 表达影响太大等造成。可以通过绘制样品聚类查看分组信息和有无异常样品。
# 如果这确实是由有意义的生物变化引起的，也可以使用下面的经验power值。
 if (is.na(power)){
   power = ifelse(nSamples<20, ifelse(type == "unsigned", 9, 18),
                  ifelse(nSamples<30, ifelse(type == "unsigned", 8, 16),
                         ifelse(nSamples<40, ifelse(type == "unsigned", 7, 14),
                                ifelse(type == "unsigned", 6, 12))       
                  )
   )
 }
power
## 9
## 04-3一步法网络构建---------
## One-step network construction and module detection##
# power: 上一步计算的软阈值
net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                       TOMType = type, minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType,
                       maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                       saveTOMFileBase = "DiffGene_TOM",
                       verbose = 3)
# 根据模块中基因数目的多少，降序排列，依次编号为 `1-最大模块数`。
# 0 (grey)表示未分入任何模块的基因。
table(net$colors)
## 0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20   21   22   23  24   
##39 4571 4415 4400 3775 3020 2556  815  738  359  246  189  184  155  153  145  136  103   76   60   59   57   48   43   38 

## 04-4 层级聚类数展示各个模块
## 灰色的为**未分类**到模块的基因。
# Convert labels to colors for plotting
moduleLabels = net$colors
moduleColors = labels2colors(moduleLabels)
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
group_traits<-group
rownames(group_traits)<-group_traits$sample
group_traits<-group_traits[rownames(dataExpr),]
group_traits<-group_traits[,-1]
group_traits<-as.data.frame(group_traits)
colnames(group_traits)<-"Group"
rownames(group_traits)<-rownames(dataExpr)
datTraits=data.frame(samples=rownames(dataExpr),subtype=group_traits)
design_traits<-model.matrix(~0+datTraits$Group)
design_traits<-as.data.frame(design_traits)
colnames(design_traits)=levels(factor(datTraits$Group))
moduleColors <- labels2colors(net$colors)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(dataExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0); ##不同颜色的模块的ME值矩阵(样本vs模块)
moduleTraitCor = cor(MEs, design_traits , use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
sizeGrWindow(10,10)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(8, 14, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(design_traits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.8,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

## 04-7 筛选强相关性模块------------
## black
module='black'
probes=colnames(dataExpr)
inModule=(moduleColors==module)
modProbes=probes[inModule]
modGenes<-as.data.frame(modProbes)
colnames(modGenes)<-'modgene'
dim(modGenes)
## 815

## 加一个greenyellow 
module2='greenyellow'
probes=colnames(dataExpr)
inModule2=(moduleColors==module2)
modProbes2=probes[inModule2]
modGenes2<-as.data.frame(modProbes2)
colnames(modGenes2)<-'modgene'
dim(modGenes2)
## 189
modGenes<-rbind(modGenes,modGenes2)
dim(modGenes)
## 1004
## 加一个turquoise 
## module3='turquoise'
## probes=colnames(dataExpr)
## inModule3=(moduleColors==module3)
## modProbes3=probes[inModule3]
## modGenes3<-as.data.frame(modProbes3)
## colnames(modGenes3)<-'modgene'
## dim(modGenes3)
## 4571
## modGenes<-rbind(modGenes,modGenes3)
## dim(modGenes)
## 5575
## 加一个cyan 
## module4='cyan'
## probes=colnames(dataExpr)
## inModule4=(moduleColors==module4)
## modProbes4=probes[inModule4]
## modGenes4<-as.data.frame(modProbes4)
## colnames(modGenes4)<-'modgene'
## dim(modGenes4)
## 153
## modGenes<-rbind(modGenes,modGenes4)
## dim(modGenes)
## 5728
## 加一个grey60 
## module5='grey60'
## probes=colnames(dataExpr)
## inModule5=(moduleColors==module5)
## modProbes5=probes[inModule5]
## modGenes5<-as.data.frame(modProbes5)
## colnames(modGenes5)<-'modgene'
## dim(modGenes5)
## 103
## modGenes<-rbind(modGenes,modGenes5)
## dim(modGenes)
## 5831
## 'salmon'
## module6='salmon'
##probes=colnames(dataExpr)
## inModule6=(moduleColors==module6)
## modProbes6=probes[inModule6]
## modGenes6<-as.data.frame(modProbes6)
## colnames(modGenes6)<-'modgene'
## dim(modGenes6)
## 155
## modGenes<-rbind(modGenes,modGenes6)
## dim(modGenes)
## 5986

# 04 HbF-erythrocyte-mod-------
setwd("/data/nas1/luchunlin/project/BJTC-290")
if (! dir.exists("./03_DEGs-mod-Erythrocyte")){
  dir.create("./03_DEGs-mod-Erythrocyte")
}
setwd("./03_DEGs-mod-Erythrocyte")
## 04-1 Hbf-mod--------
Hbf_mod<-sig_diff[rownames(sig_diff)%in%modGenes$modgene,]
## 450
library(ggvenn)
mydata1<-list(DEGs=rownames(sig_diff),modGenes=modGenes$modgene)
ggvenn(mydata1,c('DEGs','modGenes'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = T,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_linetype="solid",
       set_name_color = c("#ffb2b2","#b2e7cb"),
       text_color = 'black')


## 04-2 Hbf-mod-erythrocyte-----
erythrocyte<-read_xlsx('/data/nas1/luchunlin/project/BJTC-290/03_HbF-erythrocyte-mod/erythrocyte.xlsx')
## 去重
erythrocyte<-erythrocyte[!duplicated(erythrocyte$Genset),]
## 706
Hbf_mod_erythrocyte<-Hbf_mod[rownames(Hbf_mod)%in%erythrocyte$Genset,]
## 11
mydata2<-list(DEGs_mod=rownames(Hbf_mod),Erythrocyte=erythrocyte$Genset)
ggvenn(mydata2,c('DEGs_mod','Erythrocyte'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = T,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_linetype="solid",
       set_name_color = c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
rownames(Hbf_mod_erythrocyte)
## "EPAS1" "HBD"   "GDF2"  "CSF3"  "CSF2"  "LARS2" "ITGB4" "GATA1" "DGKE"  "HBA2"  "HBA1" 

# 05 染色体位置及表达模式-----------
setwd("/data/nas1/luchunlin/project/BJTC-290")
if (! dir.exists("./04_chromosome")){
  dir.create("./04_chromosome")
}
setwd("./04_chromosome")
## 'OmicCircos'包 
# BiocManager::install('OmicCircos')
library(OmicCircos)
## 想要的结果：从外到内依次是基因名(SNP位点)、染色体定位、热图。
## 需要的数据：（1）一个gene list（chr po gene）(2)染色体图用自带的 （3）热图（chr po gene exp）
## 先将11个基因的表达矩阵提取出来
exp_chromosome<-dat_final[rownames(dat_final)%in%rownames(Hbf_mod_erythrocyte),]
## 查找并整理染色体位置信息
## 节段数据（segment data）chrom start ened 
## 映射数据（mapping data）前面两列固定为节段名和位置 后面是表达量和甲基化
location<-read_xlsx('/data/nas1/luchunlin/project/BJTC-290/04_chromosome/location.xlsx')

## 构建热图数据(分成两组画一下)
chr_heatmap<-cbind(location,exp_chromosome)
## 实验组样本7个
chr_heatmap1<-chr_heatmap[,c(1:10)]
## 对照组样本7个
chr_heatmap2<-chr_heatmap[,c(1:3,11:17)]
## gene fusion列表
gene_fus<-cbind(location,location)
colnames(gene_fus)<-c('chr1','po1','gene1','chr2','po2','gene2')
gene_fus<-gene_fus[,c(1:3,6)]
par(mar=c(2,2,2,2))
plot(c(1,800),c(1,800),type='n',axes=F,xlab='',ylab='',main='')
circos(R=300,cir = 'hg18',W=15,type = 'chr',print.chr.lab = T,scale = T)
circos(R=330,cir = 'hg18',W=10,mapping=gene_fus,type = 'label',side = 'out',col=c('black','blue','red'),cex = 0.6)
circos(R=190,cir = 'hg18',W=100,mapping = chr_heatmap1,col.v = 4,type = 'heatmap2',cluster = F,col.bar = F,lwd = 0.5,col = 'blue')
circos(R=100,cir = 'hg18',W=100,mapping = chr_heatmap2,col.v = 4,type = 'heatmap2',cluster = F,col.bar = F,lwd = 0,col = 'blue')

# 06 GO/KEGG富集----------
setwd("/data/nas1/luchunlin/project/BJTC-290")
if (! dir.exists("./05_GO_KEGG")){
  dir.create("./05_GO_KEGG")
}
setwd("./05_GO_KEGG")
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(GOplot)
gene_transform <- bitr(rownames(Hbf_mod_erythrocyte),
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

##05-2 KEGG富集分析（气泡图）
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "hsa",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.2)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk@result,file = "KEGG.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk, showCategory=15)
kk_dot

## KEGG 弦图
kegg_cir<-cnetplot(kk,circular = TRUE, colorEdge = TRUE,showCategory = 10,
                  color_category = "#E5C494",color_gene = "#B3B3B3",
)
kegg_cir


kk2<-data.frame(Category = "ALL",ID = kk$ID,Term = kk$Description, Genes = gsub("/", ", ", kk$geneID), adj_pval = kk$p.adjust)
## 读取logFC文件
genelist<-data.frame(ID=rownames(Hbf_mod_erythrocyte),
                     logFC=Hbf_mod_erythrocyte$logFC)
rownames(genelist)<-genelist$ID
genelist<-genelist[order(genelist$logFC,decreasing = T),]
circ<-circle_dat(kk2,genelist)
#circ<-circ[c(1:52),]
# GOBubble(circ, labels = 3,table.legend =F)
# GOCircle(circ,rad1=2.5,rad2=3.5,label.size=4,nsub=10) 
#termNum<-20
#geneNum<-nrow(genelist)
chord<-chord_dat(circ,genelist)
kegg_chord<-GOChord(chord,
                    gene.order = 'logFC',
                    gene.space = 0.25,
                    gene.size = 5,
                    space = 0.01,
                    lfc.col = c('red','white','blue'),
                    process.label = 10)
kegg_chord
# 07 PPI网络构建----------
setwd("/data/nas1/luchunlin/project/BJTC-290")
if (! dir.exists("./06_PPI")){
  dir.create("./06_PPI")
}
setwd("./06_PPI")
## "EPAS1" "HBD"   "CSF3"  "CSF2"   "GATA1" "HBA2"  "HBA1" 
## 7个
# 08 hub----------
setwd("/data/nas1/luchunlin/project/BJTC-290")
if (! dir.exists("./07_hub")){
  dir.create("./07_hub")
}
setwd("./07_hub")
## 08-1 expression--------

hubgene<-data.frame(hubgene=c("EPAS1","HBD","CSF3","CSF2","GATA1","HBA2","HBA1"))
hub_exp<-dat[hubgene$hubgene,]
hub_exp2<-hub_exp
hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))
group_exp<-data.frame(Group=c(rep("β_thalassemia_HbF",49),rep("control",49)))
hub_exp2<-cbind(group_exp,hub_exp2)
##分面图形
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
# 分面图形
exp_boxplot<-ggplot(hub_exp2,aes(x = Group, y = expr, fill = Group)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "Expression",expand = c(0.1,0.1))+
  scale_x_discrete(name = "group") +
  ggtitle("Expression of hub genes") +
  theme_bw() +
    geom_signif(comparisons = list(c("control","β_thalassemia_HbF")),
                test = t.test,
                map_signif_level = T)+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 18),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11))+
  guides(fill='none')
exp_boxplot
exp_boxplot2<-exp_boxplot+scale_fill_lancet()+facet_wrap(~Symbol,scales = "free",nrow = 2)
exp_boxplot2

## 08-2 roc--------
library(pROC)
library(ggplot2)
hub_exp3<-t(hub_exp)
hub_exp3<-cbind(design,hub_exp3)
hub_exp3<-cbind(group,hub_exp3)
## "EPAS1","HBD","CSF3","CSF2","LARS2","GATA1","HBA2","HBA1"
roc_EPAS1<-roc(hub_exp3$group,hub_exp3$EPAS1,
               levels=c("control","β_thalassemia_HbF"))
roc_HBD<-roc(hub_exp3$group,hub_exp3$HBD,
               levels=c("control","β_thalassemia_HbF"))
roc_CSF3<-roc(hub_exp3$group,hub_exp3$CSF3,
             levels=c("control","β_thalassemia_HbF"))
roc_CSF2<-roc(hub_exp3$group,hub_exp3$CSF2,
              levels=c("control","β_thalassemia_HbF"))
roc_GATA1<-roc(hub_exp3$group,hub_exp3$GATA1,
               levels=c("control","β_thalassemia_HbF"))
roc_HBA2<-roc(hub_exp3$group,hub_exp3$HBA2,
               levels=c("control","β_thalassemia_HbF"))
roc_HBA1<-roc(hub_exp3$group,hub_exp3$HBA1,
              levels=c("control","β_thalassemia_HbF"))
##EPAS1
plot(roc_EPAS1,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="EPAS1 ROC curve",
     col="#FF2E63",
     legacy.axes=T,)
##HBD
plot(roc_HBD,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="HBD ROC curve",
     col="#FF2E63",
     legacy.axes=T,)
##CSF3
plot(roc_CSF3,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="CSF3 ROC curve",
     col="#FF2E63",
     legacy.axes=T,)
##CSF2
plot(roc_CSF2,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="CSF2 ROC curve",
     col="#FF2E63",
     legacy.axes=T,)
##GATA1
plot(roc_GATA1,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="GATA1 ROC curve",
     col="#FF2E63",
     legacy.axes=T,)
##HBA2
plot(roc_HBA2,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="HBA2 ROC curve",
     col="#FF2E63",
     legacy.axes=T,)
##HBA1
plot(roc_HBA1,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="HBA1 ROC curve",
     col="#FF2E63",
     legacy.axes=T,)

# 09 correlation---------
setwd("/data/nas1/luchunlin/project/BJTC-290")
if (! dir.exists("./08_correlation")){
  dir.create("./08_correlation")
}
setwd("./08_correlation")
hub_corr<-round(cor(t(hub_exp)),3)
write.table(hub_corr,file = 'correlation.xls',
            sep = '\t',
            row.names = T)
## 检验基因之间的相关性p值
## 计算相关性系数并显示基因之间的相关性。相关性系数大于0为正相关，小于0为负相关。
## p小于0.05认为相关性显著
library(ggcorrplot)
library(corrplot)
hub_p.mat<-round(cor_pmat(t(hub_exp)),3)
col1 <- colorRampPalette(c("#4169E1","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF",
                           "#FDDBC7","#F4A582","#D6604D", "#B2182B","#CC0000","#990000"))
hub_corr_plot<-corrplot(hub_corr,
                        method = "circle",
                        is.corr = T,
                        type = "lower",
                        p.mat = hub_p.mat,
                        insig = "blank",
                        outline = "white",
                        addCoef.col ="black",
                        col = col1(200))

# 10 功能相似性分析------
setwd("/data/nas1/luchunlin/project/BJTC-290")
if (! dir.exists("./09_relationship")){
  dir.create("./09_relationship")
}
setwd("./09_relationship")
library(GOSemSim)
library(org.Hs.eg.db)
bp<-godata('org.Hs.eg.db',ont = "BP",computeIC = F)
cc<-godata('org.Hs.eg.db',ont = "CC",computeIC = F)
mf<-godata('org.Hs.eg.db',ont = "MF",computeIC = F)
geneid2<-hubgene$hubgene
gene_transform2<-bitr(geneid2,fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = org.Hs.eg.db)
genelist<-gene_transform2[-3,]
simbp<-mgeneSim(genelist$ENTREZID,
                semData = bp,
                measure = "Wang",
                drop = NULL,
                combine = "BMA")
simcc<-mgeneSim(genelist$ENTREZID,
                semData = cc,
                measure = "Wang",
                drop = NULL,
                combine = "BMA")
simmf<-mgeneSim(genelist$ENTREZID,
                semData = mf,
                measure = "Wang",
                drop = NULL,
                combine = "BMA")

## 基于相似度结果，进一步计算基因在BP、CC、MF层面的集合平均值，得到最终评分。
fsim<-(simmf*simcc*simbp)^(1/3)
View(fsim)
colnames(fsim)<-genelist$SYMBOL
rownames(fsim)<-genelist$SYMBOL

## 进一步去除基因与本身基因之间的相关性，使用melt函数将宽格式数据转化为长格式数据
for (i in 1:ncol(fsim)) {
  fsim[i,i] <- NA
}
dat_sim<-melt(fsim)
View(dat_sim)
dat_sim<-dat_sim[!is.na(dat_sim$value),]
dat_sim<-dat_sim[,c(1,3)]
head(dat_sim)
## 绘制boxplot图
dat.median<-aggregate(value~Var1,dat_sim,median)
dat.mean<-aggregate(value~Var1,dat_sim,mean)
View(dat.mean)
## 根据相似性评分的平均值或中位值（根据自身需要），对其进行排序，并根据评分的高低，将基因名设置为因子（factor）格式
m<-dat.mean$value
m2<-dat.median$value
names(m)<-dat.mean$Var1
names(m2)<-dat.median$Var1
dat_sim$Var1<-factor(dat_sim$Var1,
                     levels = names(sort(m2)))
str(dat_sim)

## 使用ggplot 对分析结果进行可视化。
ggplot(dat_sim,
       aes(x=Var1,y=value,fill=factor(Var1)))+
  scale_fill_brewer(palette = "Set3")+
  geom_boxplot()+
  coord_flip()+       # 坐标轴互换
  xlab("")+ylab("")+
  theme_bw()+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 18),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13))+
  guides(fill='none')
write.table(fsim,file = 'relationship.xls',
            sep = '\t',
            quote = F)
# 11 GSEA------
setwd("/data/nas1/luchunlin/project/BJTC-290")
if (! dir.exists("./10_GSEA")){
  dir.create("./10_GSEA")
}
setwd("./10_GSEA")
allGeneSets<-tempOutput
logFCcutoff <- 0
allGeneSets$change = as.factor(
  ifelse(allGeneSets$adj.P.Val < 0.05 & abs(allGeneSets$logFC) > logFCcutoff,
         ifelse(allGeneSets$logFC > logFCcutoff,'UP','DOWN'),'NOT')
)

DEGeneSets <- subset(allGeneSets,
                     allGeneSets$adj.P.Val < 0.05 & abs(allGeneSets$logFC) > 0)

genelist <- allGeneSets$logFC
names(genelist) <- rownames(allGeneSets)
geneList <- sort(genelist, decreasing = T)
DEGeneSets <- DEGeneSets[order(DEGeneSets$adj.P.Val),]
dim(DEGeneSets)
## GSEA KEGG-----
# BiocManager::install('clusterProfiler')
library(clusterProfiler)
library(enrichplot)
kegg_set<- read.gmt("c2.cp.kegg.v7.5.1.symbols.gmt")
set.seed(1)
kegg_gsea <- GSEA(geneList, TERM2GENE = kegg_set, pvalueCutoff = 0.05)
kegg_result <- kegg_gsea@result
dim(kegg_result)
gseaplot2(kegg_gsea,c(1:5),color = c('#7B68EE','#CD3333','#20B2AA','#FF8C00','#FF6666'),
          title = 'KEGG GSEA',
          base_size = 11,
          rel_heights = c(1.5, 0.3, 0.5))
write.table(kegg_result,file = 'KEGG_gsea.xls',
            sep = '\t',
            quote = F,
            row.names = T)
# 12 tiic_cor------
setwd("/data/nas1/luchunlin/project/BJTC-290")
if (! dir.exists("./13_tiic_cor")){
  dir.create("./13_tiic_cor")
}
setwd("./13_tiic_cor")
## 基因表达数据（hub_exp） 和 免疫细胞浸润矩阵 (tiics_result)行名均为样本名称
## Spearman相关分析
##12-1  ssGSEA--------
library(GSVA)
gene_set <- read.table("/data/nas1/luchunlin/pipeline/ssGSEA/mmc3.txt",
                       header = T,
                       sep ="\t")
dat_final <- as.matrix(dat)
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
tiics_result <- ssgsea_score
hub_exp<-t(hub_exp)
tiics_exp<-t(tiics_result)
## 12-2 
### EPAS1-----
gene_EPAS1<-'EPAS1'
EPAS1<-as.numeric(hub_exp[,gene_EPAS1])
### 批量操作的具体实现过程：
### 1.设定容器,最终生成的数据放在什么地方？
correlation_EPAS1 <- data.frame()
### 2.批量把数据导出到容器
for(i in 1:length(colnames(tiics_exp))){
  ## 1.指示
  print(i)
  ## 2.计算
  dd = cor.test(as.numeric(tiics_exp[,i]),EPAS1,method="spearman")
  ## 3.填充
  correlation_EPAS1[i,1] = colnames(tiics_exp)[i]
  correlation_EPAS1[i,2] = dd$estimate
  correlation_EPAS1[i,3] = dd$p.value
}
### 修改名称
colnames(correlation_EPAS1) <- c("cell","cor","p.value")
#correlation_EPAS1<-correlation_EPAS1[order(correlation_EPAS1$cor,decreasing = T),]
correlation_EPAS1$correlation<-cut(abs(correlation_EPAS1$cor),breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0),labels = c("0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"))

### 棒棒糖图
library(ggpubr)
#display.brewer.all()

ggdotchart(correlation_EPAS1,x='cell',y='cor',
           size = 'correlation',
           sorting='descending',
           rotate = T,
           color = 'correlation',
           palette =brewer.pal(7,"GnBu"), 
           add = 'segment',                     # 添加棒子
           add.params = list(color='black',size=0.5),
           ggtheme = theme_bw(base_size = 13)+
             theme(legend.text =element_text(size = 10)),# 改变主题
           xlab="",
           ylab = 'correlation',
           title = 'EPAS1')

### HBD-----
gene_HBD<-'HBD'
HBD<-as.numeric(hub_exp[,gene_HBD])
### 批量操作的具体实现过程：
### 1.设定容器,最终生成的数据放在什么地方？
correlation_HBD <- data.frame()
### 2.批量把数据导出到容器
for(i in 1:length(colnames(tiics_exp))){
  ## 1.指示
  print(i)
  ## 2.计算
  dd = cor.test(as.numeric(tiics_exp[,i]),HBD,method="spearman")
  ## 3.填充
  correlation_HBD[i,1] = colnames(tiics_exp)[i]
  correlation_HBD[i,2] = dd$estimate
  correlation_HBD[i,3] = dd$p.value
}
### 修改名称
colnames(correlation_HBD) <- c("cell","cor","p.value")
correlation_HBD<-correlation_HBD[order(correlation_HBD$cor,decreasing = T),]
correlation_HBD$correlation<-cut(abs(correlation_HBD$cor),breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0),labels = c("0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"))

### 棒棒糖图
library(ggpubr)
#display.brewer.all()

ggdotchart(correlation_HBD,x='cell',y='cor',
           size = 'correlation',
           sorting='descending',
           rotate = T,
           color = 'correlation',
           palette =brewer.pal(8,"GnBu"), 
           add = 'segment',                     # 添加棒子
           add.params = list(color='black',size=0.5),
           ggtheme = theme_bw(base_size = 13)+
             theme(legend.text =element_text(size = 10)),# 改变主题
           xlab="",
           ylab = 'correlation',
           title = 'HBD')
### CSF3-----
gene_CSF3<-'CSF3'
CSF3<-as.numeric(hub_exp[,gene_CSF3])
### 批量操作的具体实现过程：
### 1.设定容器,最终生成的数据放在什么地方？
correlation_CSF3 <- data.frame()
### 2.批量把数据导出到容器
for(i in 1:length(colnames(tiics_exp))){
  ## 1.指示
  print(i)
  ## 2.计算
  dd = cor.test(as.numeric(tiics_exp[,i]),CSF3,method="spearman")
  ## 3.填充
  correlation_CSF3[i,1] = colnames(tiics_exp)[i]
  correlation_CSF3[i,2] = dd$estimate
  correlation_CSF3[i,3] = dd$p.value
}
### 修改名称
colnames(correlation_CSF3) <- c("cell","cor","p.value")
correlation_CSF3<-correlation_CSF3[order(correlation_CSF3$cor,decreasing = T),]
correlation_CSF3$correlation<-cut(abs(correlation_CSF3$cor),breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0),labels = c("0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"))

### 棒棒糖图
library(ggpubr)
#display.brewer.all()

ggdotchart(correlation_CSF3,x='cell',y='cor',
           size = 'correlation',
           sorting='descending',
           rotate = T,
           color = 'correlation',
           palette =brewer.pal(9,"GnBu"), 
           add = 'segment',                     # 添加棒子
           add.params = list(color='black',size=0.5),
           ggtheme = theme_bw(base_size = 13)+
             theme(legend.text =element_text(size = 10)),# 改变主题
           xlab="",
           ylab = 'correlation',
           title = 'CSF3')
### CSF2-----
gene_CSF2<-'CSF2'
CSF2<-as.numeric(hub_exp[,gene_CSF2])
### 批量操作的具体实现过程：
### 1.设定容器,最终生成的数据放在什么地方？
correlation_CSF2 <- data.frame()
### 2.批量把数据导出到容器
for(i in 1:length(colnames(tiics_exp))){
  ## 1.指示
  print(i)
  ## 2.计算
  dd = cor.test(as.numeric(tiics_exp[,i]),CSF2,method="spearman")
  ## 3.填充
  correlation_CSF2[i,1] = colnames(tiics_exp)[i]
  correlation_CSF2[i,2] = dd$estimate
  correlation_CSF2[i,3] = dd$p.value
}
### 修改名称
colnames(correlation_CSF2) <- c("cell","cor","p.value")
correlation_CSF2<-correlation_CSF2[order(correlation_CSF2$cor,decreasing = T),]
correlation_CSF2$correlation<-cut(abs(correlation_CSF2$cor),breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0),labels = c("0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"))

### 棒棒糖图
library(ggpubr)
#display.brewer.all()

ggdotchart(correlation_CSF2,x='cell',y='cor',
           size = 'correlation',
           sorting='descending',
           rotate = T,
           color = 'correlation',
           palette =brewer.pal(7,"GnBu"), 
           add = 'segment',                     # 添加棒子
           add.params = list(color='black',size=0.5),
           ggtheme = theme_bw(base_size = 13)+
             theme(legend.text =element_text(size = 10)),# 改变主题
           xlab="",
           ylab = 'correlation',
           title = 'CSF2')
### GATA1-------
gene_GATA1<-'GATA1'
GATA1<-as.numeric(hub_exp[,gene_GATA1])
### 批量操作的具体实现过程：
### 1.设定容器,最终生成的数据放在什么地方？
correlation_GATA1 <- data.frame()
### 2.批量把数据导出到容器
for(i in 1:length(colnames(tiics_exp))){
  ## 1.指示
  print(i)
  ## 2.计算
  dd = cor.test(as.numeric(tiics_exp[,i]),GATA1,method="spearman")
  ## 3.填充
  correlation_GATA1[i,1] = colnames(tiics_exp)[i]
  correlation_GATA1[i,2] = dd$estimate
  correlation_GATA1[i,3] = dd$p.value
}
### 修改名称
colnames(correlation_GATA1) <- c("cell","cor","p.value")
correlation_GATA1<-correlation_GATA1[order(correlation_GATA1$cor,decreasing = T),]
correlation_GATA1$correlation<-cut(abs(correlation_GATA1$cor),breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0),labels = c("0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"))

### 棒棒糖图
library(ggpubr)
#display.brewer.all()

ggdotchart(correlation_GATA1,x='cell',y='cor',
           size = 'correlation',
           sorting='descending',
           rotate = T,
           color = 'correlation',
           palette =brewer.pal(7,"GnBu"), 
           add = 'segment',                     # 添加棒子
           add.params = list(color='black',size=0.5),
           ggtheme = theme_bw(base_size = 13)+
             theme(legend.text =element_text(size = 10)),# 改变主题
           xlab="",
           ylab = 'correlation',
           title = 'GATA1')
### HBA2----
gene_HBA2<-'HBA2'
HBA2<-as.numeric(hub_exp[,gene_HBA2])
### 批量操作的具体实现过程：
### 1.设定容器,最终生成的数据放在什么地方？
correlation_HBA2 <- data.frame()
### 2.批量把数据导出到容器
for(i in 1:length(colnames(tiics_exp))){
  ## 1.指示
  print(i)
  ## 2.计算
  dd = cor.test(as.numeric(tiics_exp[,i]),HBA2,method="spearman")
  ## 3.填充
  correlation_HBA2[i,1] = colnames(tiics_exp)[i]
  correlation_HBA2[i,2] = dd$estimate
  correlation_HBA2[i,3] = dd$p.value
}
### 修改名称
colnames(correlation_HBA2) <- c("cell","cor","p.value")
correlation_HBA2<-correlation_HBA2[order(abs(correlation_HBA2$cor)),]
correlation_HBA2$correlation<-cut(abs(correlation_HBA2$cor),breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0),labels = c("0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"))

### 棒棒糖图
library(ggpubr)
#display.brewer.all()

ggdotchart(correlation_HBA2,x='cell',y='cor',
           size = 'correlation',
           sorting='descending',
           rotate = T,
           color = 'correlation',
           palette =brewer.pal(7,"GnBu"), 
           add = 'segment',                     # 添加棒子
           add.params = list(color='black',size=0.5),
           ggtheme = theme_bw(base_size = 13)+
             theme(legend.text =element_text(size = 10)),# 改变主题
           xlab="",
           ylab = 'correlation',
           title = 'HBA2')
### HBA1-------
gene_HBA1<-'HBA1'
HBA1<-as.numeric(hub_exp[,gene_HBA1])
### 批量操作的具体实现过程：
### 1.设定容器,最终生成的数据放在什么地方？
correlation_HBA1 <- data.frame()
### 2.批量把数据导出到容器
for(i in 1:length(colnames(tiics_exp))){
  ## 1.指示
  print(i)
  ## 2.计算
  dd = cor.test(as.numeric(tiics_exp[,i]),HBA1,method="spearman")
  ## 3.填充
  correlation_HBA1[i,1] = colnames(tiics_exp)[i]
  correlation_HBA1[i,2] = dd$estimate
  correlation_HBA1[i,3] = dd$p.value
}
### 修改名称
colnames(correlation_HBA1) <- c("cell","cor","p.value")
correlation_HBA1<-correlation_HBA1[order(abs(correlation_HBA1$cor)),]
correlation_HBA1$correlation<-cut(abs(correlation_HBA1$cor),breaks = c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0),labels = c("0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"))

### 棒棒糖图
library(ggpubr)
#display.brewer.all()

ggdotchart(correlation_HBA1,x='cell',y='cor',
           size = 'correlation',
           sorting='descending',
           rotate = T,
           color = 'correlation',
           palette =brewer.pal(7,"GnBu"), 
           add = 'segment',                     # 添加棒子
           add.params = list(color='black',size=0.5),
           ggtheme = theme_bw(base_size = 13)+
             theme(legend.text =element_text(size = 10)),# 改变主题
           xlab="",
           ylab = 'correlation',
           title = 'HBA1')
correlation_HBA2<-correlation_HBA2[order(correlation_HBA2$cor,decreasing = T),]
# 13 网络构建------
setwd("/data/nas1/luchunlin/project/BJTC-290")
if (! dir.exists("./12_network")){
  dir.create("./12_network")
}
setwd("./12_network")

# 14 ANTXR1相关性------
setwd("/data/nas1/luchunlin/project/BJTC-290")
if (! dir.exists("./13_ANTXR1_cor")){
  dir.create("./13_ANTXR1_cor")
}
setwd("./13_ANTXR1_cor")
rownames(Hbf_mod_erythrocyte)
cor_gene<-data.frame(gene=c("EPAS1","HBD","DGKE","CSF3","CSF2","GATA1","GDF2","LARS2","HBA2","HBA1","ANTXR1"))
cor_exp<-dat[cor_gene$gene,]
corr<-round(cor(t(cor_exp)),3)
write.table(corr,file = 'correlation.xls',
            sep = '\t',
            row.names = T)
## 检验基因之间的相关性p值
## 计算相关性系数并显示基因之间的相关性。相关性系数大于0为正相关，小于0为负相关。
## p小于0.05认为相关性显著
library(ggcorrplot)
library(corrplot)
hub_p.mat<-round(cor_pmat(t(cor_exp)),3)
col1 <- colorRampPalette(c("#4169E1","#2166AC","#4393C3","#92C5DE","#D1E5F0","#FFFFFF",
                           "#FDDBC7","#F4A582","#D6604D", "#B2182B","#CC0000","#990000"))
corr_plot<-corrplot(corr,
                    method = "circle",
                    is.corr = T,
                    type = "lower",
                    p.mat = hub_p.mat,
                    insig = 'blank',
                    outline = "white",
                    addCoef.col ="black",
                    col = col1(200))
# 15 logistic-----
setwd("/data/nas1/luchunlin/project/BJTC-290")
if (! dir.exists("./14_logistic")){
  dir.create("./14_logistic")
}
setwd("./14_logistic")
## 15-1 模型构建-------
##构建逻辑回归模型，并评估逻辑回归模型诊断疾病的性能。这里不做lasso
## 混淆矩阵热图  ROC曲线  PR曲线  校准曲线
## 
type<-data.frame(sample=group$sample,group=group$group)
type$group<-ifelse(type$group=='β_thalassemia_HbF',1,0)
hub_exp<-as.data.frame(hub_exp)
Affairs<-data.frame(hub_exp,type)
#Affairs<-Affairs[,-8]
y<-as.matrix(type$group)
x<-as.matrix(hub_exp)
set.seed(15)
library(glmnet)
fit = glmnet(x, y, alpha=1,family="binomial") #这里alpha=1为LASSO回归，如果等于0就是岭回归
#参数 family 规定了回归模型的类型：
# family="gaussian" 适用于一维连续因变量（univariate）
# family="mgaussian" 适用于多维连续因变量（multivariate）
# family="poisson" 适用于非负次数因变量（count）
# family="binomial" 适用于二元离散因变量（binary）
# family="multinomial" 适用于多元离散因变量（category）
#我们这里结局指标是2分类变量，所以使用binomial
pdf("stromallambda.pdf",w=7,h=7)
png("stromallambda.png",w=700,h=700)
plot(fit, xvar="lambda", label=TRUE)
dev.off()
pdf("stromalcvfit.pdf",w=9,h=8)
png("stromalcvfit.png",w=900,h=800)
set.seed(2937)

fit_cv <- cv.glmnet(x, y, alpha=1, family = 'binomial')
plot(fit_cv)
dev.off()
fit_cv
coefficient<-coef(fit_cv,fit_cv$lambda.min)

coefficient<-as.matrix(coefficient)

write.table(coefficient,'coefficient1,xls',sep='\t',quote=F)

set.seed(355423)
fit.ful<-glm(group~EPAS1+HBA2,data = Affairs)
summary(fit.ful)

Affairs$prob <- predict(fit.ful,
                        newdata=Affairs,
                        type="response")

set.seed(355423)
fit<-glm(group~.,data = Affairs)
#fit<-lrm(group~EPAS1+HBA2,data = Affairs,x=T,y=T)
fit
## CSF3
summary(fit)
coefficient<-coef(fit,fit$coefficients)
coefficient<-as.matrix(coefficient)
fit2<-glm(group~CSF3,data = Affairs)
fit2
set.seed(12332)
Affairs$prob <- predict(fit2,
                        newdata=Affairs,
                        type="response")

## 15-2 混淆矩阵热图------
library(caret)
Actual<-factor(Affairs$group,levels = c(1,0),labels = c('β thalassemia HbF','control'))
Predict<-ifelse(Affairs$prob>median(Affairs$prob),1,0)
Predict<-factor(Predict,levels = c(1,0),labels = c('β thalassemia HbF','control'))
xtab<-table(Predict,Actual)
xtab
confusionMatrix(xtab)
pheatmap(mat = xtab,
         cluster_cols = F,
         cluster_rows = F,
         fontsize = 15,
         display_numbers = T,
         number_format = '%.0f',
         fontsize_number = 15)
## 15-3 ROC曲线--------
roc1 <- roc(Affairs$group, Affairs$prob)
plot(roc1, 
     print.auc=TRUE, #设置是否添加AUC值标签
     auc.polygon=TRUE, #设置是否添加AUC值面积多边形
     max.auc.polygon=TRUE, #设置是否添加最大AUC值面积多边形
     auc.polygon.col="skyblue", #设置AUC值面积多边形的填充色
     grid=c(0.1, 0.2), #添加网格线
     grid.col=c("green", "red"), #设置网格线颜色
     print.thres=TRUE)


## 15-4 DCA曲线 --------------------
##source("dca.R")#需要把这个代码放到加载目录里，文章最后会贴出，用的时候直接复查到.R后缀的文件并改名为dca.R就行
## install.packages("rmda")
library(rmda)
pdf("DCA.pdf",width=8,height=8)
png("DCA.png",width=800,height=800)
complex<-decision_curve(group~EPAS1+HBD+CSF3+CSF2+GATA1+HBA2+HBA1,data = Affairs,family = binomial(link ='logit'),
                        thresholds = seq(0,1, by = 0.01),
                        confidence.intervals= 0.95,
                        study.design = 'case-control',
                        population.prevalence= 0.3
)


plot_decision_curve(complex,
                    curve.names=c('logistic'),
                    cost.benefit.axis =FALSE,col= c('black'),
                    confidence.intervals=FALSE,
                    standardize = FALSE)
dev.off()

## 15-5 PCR曲线 ---------------------------------------------------------------------
#install.packages('modEvA')
library(modEvA)
png("pcr.png",width=800,height=700)
#pdf("pcr.pdf",width=8,height=7)
aupr=AUC(obs=Affairs$group,pred=Affairs$prob,curve = "PR", simplif=TRUE, main = "PR curve")
plot(aupr)
dev.off()

## 15-6 校准曲线 ---------------------------------------------------------------------
tibble_ex <- tibble(
  group = Affairs$group,
  pred = Affairs$prob
)

tibble_ex %>%
  arrange(tibble_ex) %>%
  ggplot(aes(x = pred, y = group)) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, by = 0.2)) +
  geom_smooth(formula = y~x,aes(x = pred, y = group), color = "darkgrey", se = F, method = "loess") + 
  # you can use stat_smooth in place of geom_smooth
  geom_abline()


library(rms)
Affairs<-Affairs[-c(8,10)]
Affairs$EPAS1<-as.factor(Affairs$EPAS1)
Affairs$HBD<-as.factor(Affairs$HBD)
Affairs$CSF3<-as.factor(Affairs$CSF3)
Affairs$CSF2<-as.factor(Affairs$CSF2)
Affairs$GATA1<-as.factor(Affairs$GATA1)
Affairs$HBA2<-as.factor(Affairs$HBA2)
formula<-as.formula(group~EPAS1+HBD+CSF3+CSF2+GATA1+HBA2+HBA1)
rms_fit <- lrm(formula,data = Affairs,x=TRUE,y=TRUE,'tol = 1e-09',)
rms_fit

cal1 <- calibrate(fit, cmethod='hare', method='boot', B=1000,data=vad)#建模组中绘制校准曲线
plot(cal1,xlim=c(0,1.0),ylim=c(0,1.0))#打印出校准曲线



