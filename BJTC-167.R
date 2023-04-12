# rm(list = ls())
# 01 获取数据集
setwd("/data/nas1/luchunlin/project/BJTC-167")
if (! dir.exists("./00_rawdata")){
dir.create("./00_rawdata")
}
setwd("./00_rawdata")
library(GEOquery)
library(Biobase)
library(limma)
## 两个训练集，111006为a ，111016为b
gset_a<-getGEO("GSE111006",
              destdir = '.',
              GSEMatrix = T,
              getGPL = F)
a=gset_a[[1]]
pd_a<-pData(a)
gset_b<-getGEO("GSE111016",
              destdir = '.',
              GSEMatrix = T,
              getGPL = F)
b=gset_b[[1]]
pd_b<-pData(b)
exp_a<-read_csv(file = "GSE111006_lcpmy1_hss_forGEO.csv",
                col_names = T)
exp_a<-as.data.frame(exp_a)
rownames(exp_a)<-exp_a[,1]
exp_a<-exp_a[,-1]
exp_b<-read_csv(file = "GSE111016_lcpmy1_sss_forGEO.csv",
                col_names = T)
exp_b<-as.data.frame(exp_b)
rownames(exp_b)<-exp_b[,1]
exp_b<-exp_b[,-1]
## 验证集111010
gset_va<-getGEO("GSE111010",
                destdir = '.',
                GSEMatrix = T,
                getGPL = F)
va=gset_va[[1]]
pd_va<-pData(va)
exp_va<-read.csv(file = "GSE111010_lcpmy1_jss_forGEO.csv")
exp_va<-as.data.frame(exp_va)
rownames(exp_va)<-exp_va[,1]
exp_va<-exp_va[,-1]
## 将Ensembl换成Symbol
library(org.Hs.eg.db)
library(clusterProfiler)
id_a<-rownames(exp_a)
#gene_transform_a<-mapIds(org.Hs.eg.db,keys = id_a,keytype = "ENSEMBL",column = "SYMBOL")
gene_transform_a<-bitr(id_a,fromType = "ENSEMBL",
                       toType = "SYMBOL",
                       OrgDb = org.Hs.eg.db)
gene_transform_a<-as.data.frame(gene_transform_a)
length(unique(gene_transform_a$SYMBOL))
gene_transform_a<-gene_transform_a[!duplicated(gene_transform_a$SYMBOL),]
gene_transform_a<-gene_transform_a[!duplicated(gene_transform_a$ENSEMBL),]
exp_a<-exp_a[rownames(exp_a)%in%gene_transform_a$ENSEMBL,]
gene_transform_a<-gene_transform_a[gene_transform_a$ENSEMBL%in%rownames(exp_a),]
rownames(exp_a)<-gene_transform_a$SYMBOL
## 然后换b数据集
id_b<-rownames(exp_b)
#gene_transform_a<-mapIds(org.Hs.eg.db,keys = id_a,keytype = "ENSEMBL",column = "SYMBOL")
gene_transform_b<-bitr(id_b,fromType = "ENSEMBL",
                       toType = "SYMBOL",
                       OrgDb = org.Hs.eg.db)
gene_transform_b<-as.data.frame(gene_transform_b)
length(unique(gene_transform_b$SYMBOL))
gene_transform_b<-gene_transform_b[!duplicated(gene_transform_b$SYMBOL),]
gene_transform_b<-gene_transform_b[!duplicated(gene_transform_b$ENSEMBL),]
exp_b<-exp_b[rownames(exp_b)%in%gene_transform_b$ENSEMBL,]
gene_transform_b<-gene_transform_b[gene_transform_b$ENSEMBL%in%rownames(exp_b),]
rownames(exp_b)<-gene_transform_b$SYMBOL
## 然后换验证集
id_va<-rownames(exp_va)
#gene_transform_a<-mapIds(org.Hs.eg.db,keys = id_a,keytype = "ENSEMBL",column = "SYMBOL")
gene_transform_va<-bitr(id_va,fromType = "ENSEMBL",
                        toType = "SYMBOL",
                        OrgDb = org.Hs.eg.db)
gene_transform_va<-as.data.frame(gene_transform_va)
length(unique(gene_transform_va$SYMBOL))
gene_transform_va<-gene_transform_va[!duplicated(gene_transform_va$SYMBOL),]
gene_transform_va<-gene_transform_b[!duplicated(gene_transform_va$ENSEMBL),]
exp_va<-exp_va[rownames(exp_va)%in%gene_transform_va$ENSEMBL,]
gene_transform_va<-gene_transform_va[gene_transform_va$ENSEMBL%in%rownames(exp_va),]
rownames(exp_va)<-gene_transform_va$SYMBOL
exp_a<-exp_a[,-c(1,9:11,13,24,30,36)]
exp_b<-exp_b[,-3]
write.table(exp_a,file = "GSE111006_EXP.xls",
            quote = F,
            sep = '\t',
            row.names = T)
write.table(exp_b,file = "GSE111016_EXP.xls",
            quote = F,
            sep = '\t',
            row.names = T)
write.table(exp_va,file = "GSE111010_EXP.xls",
            quote = F,
            sep = '\t',
            row.names = T)
# 02 DEGs
setwd("/data/nas1/luchunlin/project/BJTC-167")
if (! dir.exists("./01_DEGs")){
dir.create("./01_DEGs")
}
setwd("./01_DEGs")
# 01-1差异分析
## limma包pvalue<0.05,|logFC|>1  火山图，热图，2个数据集
group_a<-read_xlsx("/data/nas1/luchunlin/project/BJTC-167/01_DEGs/Group_a.xlsx")
design_a<-as.matrix(group_a)
type_a<-design_a[,2]
design_a <- model.matrix(~ -1+factor(type_a,levels=c('Health','Sarcopenia')))
colnames(design_a)<-c('Health','Sarcopenia')
#构建差异表达矩阵
library(limma)
#DGElist_a<-DGEList(counts=exp_a,group=group_a$Group)
#keep_gene_a<-rowSums(cpm(DGElist_a)>1)>=2  #自定义
#table(keep_gene_a)
#DGElist_a<-DGElist_a[keep_gene_a,,keep.lib.sizes=F]
#DGElist_a<-calcNormFactors(DGElist_a)
#v_a<-voom(DGElist_a,design_a,plot = T,normalize.method = "quantile")
# 对每一个基因进行线性模型构建
fit_a=lmFit(exp_a,design_a)
# 构建比较矩阵
contrast.matrix_a=makeContrasts(ControlVSMG=Health-Sarcopenia,levels = design_a)
#构建芯片数据线性模型，计算估计的相关系数和标准差
fit2_a=contrasts.fit(fit_a,contrast.matrix_a)
# 基于贝叶斯计算T值，F值和log-odds
fit2_a<-eBayes(fit2_a,0.01)
DEGs_a=topTable(fit2_a,coef =1,n=Inf)
DEGs_a=na.omit(DEGs_a)
head(DEGs_a)
write.table(DEGs_a,file = "DEG_a.xls",
            quote = F,
            sep = "\t",
            row.names = T)
# 筛选差异基因
logFC_cutoff <- 0
DEGs_a$change = as.factor(
ifelse(DEGs_a$P.Value < 0.05 & abs(DEGs_a$logFC) > logFC_cutoff,
       ifelse(DEGs_a$logFC > logFC_cutoff ,'UP','DOWN'),'NOT'))
sig_diff_a <- subset(DEGs_a,
                     DEGs_a$P.Value < 0.05& abs(DEGs_a$logFC) >= 0)
dim(DEGs_a)
dim(sig_diff_a)
summary(sig_diff_a$change)
# DOWN  NOT   UP
# 779    0  740
group_b<-read_xlsx("/data/nas1/luchunlin/project/BJTC-167/01_DEGs/Group_b.xlsx")
design_b<-as.matrix(group_b)
type_b<-design_b[,2]
design_b <- model.matrix(~ -1+factor(type_b,levels=c('Health','Sarcopenia')))
colnames(design_b)<-c('Health','Sarcopenia')
#构建差异表达矩阵
library(limma)
# 对每一个基因进行线性模型构建
fit_b=lmFit(exp_b,design_b)
# 构建比较矩阵
contrast.matrix_b=makeContrasts(ControlVSMG=Health-Sarcopenia,levels = design_b)
#构建芯片数据线性模型，计算估计的相关系数和标准差
fit2_b=contrasts.fit(fit_b,contrast.matrix_b)
# 基于贝叶斯计算T值，F值和log-odds
fit2_b<-eBayes(fit2_b,0.01)
DEGs_b=topTable(fit2_b,coef =1,n=Inf)
DEGs_b=na.omit(DEGs_b)
head(DEGs_b)
write.table(DEGs_b,file = "DEG_b.xls",
            quote = F,
            sep = "\t",
            row.names = T)
# 筛选差异基因
logFC_cutoff <- 0
DEGs_b$change = as.factor(
ifelse(DEGs_b$P.Value < 0.05 & abs(DEGs_b$logFC) > logFC_cutoff,
       ifelse(DEGs_b$logFC > logFC_cutoff ,'UP','DOWN'),'NOT'))
sig_diff_b <- subset(DEGs_b,
            DEGs_b$P.Value < 0.05& abs(DEGs_b$logFC) >= 0)
dim(DEGs_b)
dim(sig_diff_b)
summary(sig_diff_b$change)
#DOWN  NOT   UP
#1012    0  1092
# 01-2绘制火山图
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(Ipaper)
library(scales)
library(ggrepel)
volcano_plot_a<- ggplot(data = DEGs_a,
                        aes(x = logFC,
                            y = -log10(P.Value),
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
                                size=8),
      plot.title = element_text(hjust = 0.5),
      axis.text.x = element_text(face = "bold",
                                 color = "black",
                                 size = 12),
      axis.text.y = element_text(face = "bold",
                                 color = "black",
                                 size = 12),
      axis.title.x = element_text(face = "bold",
                                 color = "black",
                                 size = 12),
      axis.title.y = element_text(face = "bold",
                                 color = "black",
                                 size = 12)) +
labs(x = "log (Fold Change)",
     y = "-log10 (P.Val)")
volcano_plot_a
volcano_plot_b<- ggplot(data = DEGs_b,
                        aes(x = logFC,
                            y = -log10(P.Value),
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
                                   size=8),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face = "bold",
                                   color = "black",
                                   size = 12),
        axis.text.y = element_text(face = "bold",
                                   color = "black",
                                   size = 12),
        axis.title.x = element_text(face = "bold",
                                    color = "black",
                                    size = 12),
        axis.title.y = element_text(face = "bold",
                                    color = "black",
                                    size = 12)) +
  labs(x = "log (Fold Change)",
       y = "-log10 (P.Val)")
volcano_plot_b
# 01-3 绘制热图
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(gplots)
group<-group_a
group<-as.data.frame(group)
rt<-exp_a
rownames(group)<-colnames(rt)
## 排序，让rt和group顺序保持一致
colnames(rt)<-group$Group
group<-group[order(group$Group),]
class(group)
rt<-rt[,order(colnames(rt))]
colnames(rt)<-group$Sample
group<-group[,2]
group<-as.data.frame(group)
rownames(group)<-colnames(rt)
heat<-rt[rownames(rt)%in%
         c(head(rownames(subset(sig_diff_a,sig_diff_a$logFC>0)),10),head(rownames(subset(sig_diff_a,sig_diff_a$logFC<0)),10)),]
x<-t(scale(t(heat)))
ann_colors<-list(
  Group = c(Health="lightblue",Sarcopenia="darkorange"))
pheatmap(mat=x,
         annotation_col = group,
         color=bluered(100),
         scale = "row",
         annotation_colors = ann_colors,
         fontsize = 10,
         show_colnames = F,
         cluster_cols = F,
         cluster_rows = T)
## b数据集
group2<-group_b
group2<-as.data.frame(group2)
rt2<-exp_b
rownames(group2)<-colnames(rt2)
## 排序，让rt和group顺序保持一致
colnames(rt2)<-group2$Group
group2<-group2[order(group2$Group),]
class(group2)
rt2<-rt2[,order(colnames(rt2))]
colnames(rt2)<-group2$Sample
group2<-group2[,2]
group2<-as.data.frame(group2)
rownames(group2)<-colnames(rt2)
heat2<-rt2[rownames(rt2)%in%
           c(head(rownames(subset(sig_diff_b,sig_diff_b$logFC>0)),10),head(rownames(subset(sig_diff_b,sig_diff_b$logFC<0)),10)),]
x2<-t(scale(t(heat2)))
ann_colors2<-list(
 Group = c(Health="lightblue",Sarcopenia="darkorange"))
#绘制基因的热图，并调整参数
pheatmap(mat=x2,
         annotation_col = group2,
         color=bluered(100),
         scale = "row",
         annotation_colors = ann_colors2,
         fontsize = 10,
         show_colnames = F,
         cluster_cols = F,
         cluster_rows = T)
# 03 WGCNA
setwd("/data/nas1/luchunlin/project/BJTC-167")
if (! dir.exists("./02_WGCNA")){
dir.create("./02_WGCNA")
}
setwd("./02_WGCNA")
library(WGCNA)
library(reshape2)
library(stringr)
###数据集b
options(stringsAsFactors = F)
enableWGCNAThreads()
exprMat_b<-exp_b
dim(exprMat_b)
type = "unsigned"
corType = "pearson"
corFnc = ifelse(corType=="pearson", cor, bicor)
# 对二元变量，如样本性状信息计算相关性时，
# 或基因表达严重依赖于疾病状态时，需设置下面参数
maxPOutliers = ifelse(corType=="pearson",1,0.05)
# 关联样品性状的二元变量时，设置
robustY = ifelse(corType=="pearson",T,F)
dataExpr_b <- exprMat_b[rownames(exp_b),]
# 03-1 数据筛选
## 筛选中位绝对偏差（MAD）前75%的基因，至少MAD大于0.01
## 筛选后会降低运算量，也会失去部分信息
## 也可不做筛选，使MAD大于0即可
m.mad <- apply(dataExpr_b,1,mad)
dataExprVar <- dataExpr_b[which(m.mad >
                                max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
## 转换为样品在行，基因在列的矩阵
dataExpr_b <- as.data.frame(t(dataExprVar))
## 检测缺失值
gsg = goodSamplesGenes(dataExpr_b, verbose = 3)
gsg$allOK
if (!gsg$allOK){
# Optionally, print the gene and sample names that were removed:
if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:",
                      paste(names(dataExpr_b)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:",
    paste(rownames(dataExpr_b)[!gsg$goodSamples], collapse = ",")));
# Remove the offending genes and samples from the data:
dataExpr_b = dataExpr_b[gsg$goodSamples, gsg$goodGenes]
}
nGenes = ncol(dataExpr_b)
nSamples = nrow(dataExpr_b)
dim(dataExpr_b)
# [1]   39 10881
# 03-2 软阈值筛选
## 样本聚类，查看是否有利群样本
sampleTree = hclust(dist(dataExpr_b), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
#abline(h = 15, col = "red")
# Determine cluster under the line
# 剪枝算法，cutHeight修剪树枝的高度，minSize集群最小数
#clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
#table(clust)
# clust 1 contains the samples we want to keep.
#keepSamples = (clust==1)
## 符合要求的数据
#dataExpr_b = dataExpr_b[keepSamples,]
## 提取列
nGenes = ncol(dataExpr_b)
## 提取行
nSamples = nrow(dataExpr_b)
# 设置网络构建参数选择范围，计算无尺度分布拓扑矩阵
powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(dataExpr_b, powerVector=powers,
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
##4
##一步法网络构建：One-step network construction and module detection##
# power: 上一步计算的软阈值
net = blockwiseModules(dataExpr_b, power = power, maxBlockSize = nGenes,
                       TOMType = type, minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs=TRUE, corType = corType,
                       maxPOutliers=maxPOutliers, loadTOMs=TRUE,
                       saveTOMFileBase = "DiffGene_TOM",
                       verbose = 3)
# 根据模块中基因数目的多少，降序排列，依次编号为 `1-最大模块数`。
# **0 (grey)**表示**未**分入任何模块的基因。
table(net$colors)
#   0    1    2    3    4    5    6    7    8    9   10   11   12   13   4   15   16   17   18   19   20   21   22
#1093 2016 1658 424 1294  754  330  286  273  262  253  161  159  134  133  127  120   90   87   73   56   50   48
# 03-4 层级聚类树展示各个模块
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
group_traits<-group_b
sample<-group_traits[,1]
sample<-as.data.frame(sample)
group_traits<-group_traits[,-1]
group_traits<-as.data.frame(group_traits)
rownames(group_traits)<-sample$sample
colnames(group_traits)<-"Group"
datTraits=data.frame(samples=rownames(dataExpr_b),subtype=group_traits)
design_traits<-model.matrix(~0+datTraits$Group)
design_traits<-as.data.frame(design_traits)
colnames(design_traits)=levels(factor(datTraits$Group))
moduleColors <- labels2colors(net$colors)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(dataExpr_b, moduleColors)$eigengenes
MEs = orderMEs(MEs0); ##不同颜色的模块的ME值矩阵(样本vs模块)
moduleTraitCor = cor(MEs, design_traits , use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
sizeGrWindow(10,10)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
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
# 03-6 绘制模块间的相关性热图
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
sizeGrWindow(9,5)
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
#03-7 导出模块基因
### 计算模块与基因的相关性矩阵
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(dataExpr_b, MEs, use = "p"));
## 算出每个模块跟基因的皮尔森相关系数矩阵
## MEs是每个模块在每个样本里面的值
## datExpr是每个基因在每个样本的表达量
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
### 计算性状与基因的相关性矩阵
## 只有连续型性状才能只有计算
## 这里把是否属于  表型这个变量用0,1进行数值化。
Sarcopenia = as.data.frame(design_traits[,2]);
names(Sarcopenia) = "Sarcopenia"
geneTraitSignificance = as.data.frame(cor(dataExpr_b, Sarcopenia, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(Sarcopenia), sep="");
names(GSPvalue) = paste("p.GS.", names(Sarcopenia), sep="")
### 提取模块的基因
module="brown"
probes=colnames(dataExpr_b)
inModule=(moduleColors==module)
modProbes=probes[inModule]
modProbes<-as.data.frame(modProbes)

## 424
## 数据集b交集
sig_diff_b2<-sig_diff_b[rownames(sig_diff_b)%in%modProbes$modProbes,]

##728
## 绘制韦恩图
library(VennDiagram)
grid.newpage();
venn_plot<-draw.pairwise.venn(
  area1 = 2104,
  area2 = 424,
  cross.area = 728,
  category = c('DEGs','MOD gene'),
  fill = c("#FFC0CB","#ADD8E6"),
  cex = 3,
  cat.cex = 2,
  cat.pos = c(185,150),
  cat.dist = 0.05,
  cat.just = list(c(0,1),c(1,1)),
  ext.pos = 30,
  ext.dist = -0.05,
  ext.length = 0.85,
  ext.line.lwd = 2,
  ext.line.lty = "dashed"
  )
# 04 肌肉减少症中衰老相关差异表达基因（DEARG）的鉴定
setwd("/data/nas1/luchunlin/project/BJTC-167")
if (! dir.exists("./03_DEARG")){
  dir.create("./03_DEARG")
  }
setwd("./03_DEARG")
## HAGR数据库得到衰老相关基因（aging-related gene），将衰老相关基因与GSE111016数据集DEGs取交集。得到DEARG
ARGs<-read.csv(file = "genage_human.csv")
DEARGs<-sig_diff_b2[rownames(sig_diff_b2)%in%ARGs$symbol,]  # 16个
write.table(DEARGs,file = "DEARGs_list.xls",
            quote = F,
            sep = "\t")
## 绘制韦恩图
library(VennDiagram)
grid.newpage();
DEARGs_venn<-draw.pairwise.venn(
  area1 = 728,
  area2 = 307,
  cross.area = 16,
  category = c('Sarcopenia DEGs','ARGs'),
  fill = c("#FFC0CB","#ADD8E6"),
  cex = 3,
  cat.cex = 2,
  cat.pos = c(200,160),
  cat.dist = 0.05,
  cat.just = list(c(0,1),c(1,1)),
  ext.pos = 30,
  ext.dist = -0.05,
  ext.length = 0.85,
  ext.line.lwd = 2,
  ext.line.lty = "dashed"
)
# 05 DEARG GO/KEGG富集
setwd("/data/nas1/luchunlin/project/BJTC-167")
if (! dir.exists("./04_GO_KEGG")){
dir.create("./04_GO_KEGG")
}
setwd("./04_GO_KEGG")
##05-1 GO 富集分析（条形图）
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(GOplot)
diff_ARGs<-DEARGs
diff_gene_names <- rownames(diff_ARGs)
gene_transform <- bitr(diff_gene_names,
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
go_bar <- barplot(ego, showCategory=5, split="ONTOLOGY") +
  facet_grid(ONTOLOGY ~ ., scales = "free")
go_bar
##05-2 KEGG富集分析（气泡图）
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "hsa",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 2)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk@result,file = "KEGG.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk, showCategory=25)
kk_dot
# 06 PPI网络筛选Hub基因
setwd("/data/nas1/luchunlin/project/BJTC-167")
if (! dir.exists("./05_PPI")){
dir.create("./05_PPI")
}
setwd("./05_PPI")
# 07 诊断基因筛选
setwd("/data/nas1/luchunlin/project/BJTC-167")
if (! dir.exists("./06_hub")){
  dir.create("./06_hub")
  }
setwd("./06_hub")
## 07-1 ROC曲线验证hub基因对肌肉减少症的诊断价值
##将5个hub基因的表达矩阵提取出来
#hubgene<-rownames(DEARGs)
hubgene<-c("VEGFA","FOXO3","FOXO1","STAT5A","PPARGC1A")
hubgene<-as.data.frame(hubgene)
exp_hub<-exp_b[rownames(exp_b)%in%hubgene$hubgene,]
library(pROC)
library(ggplot2)
exp_hub2<-t(exp_hub)
exp_hub2<-cbind(design_b,exp_hub2)
exp_hub2<-as.data.frame(exp_hub2)
exp_hub2<-cbind(group_traits,exp_hub2)
exp_hub2<-as.data.frame(exp_hub2)
## 绘制ROC曲线
roc_HSPA8<-roc(exp_hub2$Group,exp_hub2$HSPA8,
               levels=c("Health","Sarcopenia"))
roc_DDIT3<-roc(exp_hub2$Group,exp_hub2$DDIT3,
               levels=c("Health","Sarcopenia"))
roc_AIFM1<-roc(exp_hub2$Group,exp_hub2$AIFM1,
               levels=c("Health","Sarcopenia"))
roc_PLAU<-roc(exp_hub2$Group,exp_hub2$PLAU,
               levels=c("Health","Sarcopenia"))
roc_RGN<-roc(exp_hub2$Group,exp_hub2$RGN,
               levels=c("Health","Sarcopenia"))
roc_ERCC8<-roc(exp_hub2$Group,exp_hub2$ERCC8,
               levels=c("Health","Sarcopenia"))
roc_VEGFA<-roc(exp_hub2$Group,exp_hub2$VEGFA,
              levels=c("Health","Sarcopenia"))
roc_FOXO3<-roc(exp_hub2$Group,exp_hub2$FOXO3,
              levels=c("Health","Sarcopenia"))
roc_FOXO1<-roc(exp_hub2$Group,exp_hub2$FOXO1,
              levels=c("Health","Sarcopenia"))
roc_STAT5A<-roc(exp_hub2$Group,exp_hub2$STAT5A,
              levels=c("Health","Sarcopenia"))
roc_PPARGC1A<-roc(exp_hub2$Group,exp_hub2$PPARGC1A,
              levels=c("Health","Sarcopenia"))
roc_RPA1<-roc(exp_hub2$Group,exp_hub2$RPA1,
               levels=c("Health","Sarcopenia"))
roc_TAF1<-roc(exp_hub2$Group,exp_hub2$TAF1,
               levels=c("Health","Sarcopenia"))
roc_SIN3A<-roc(exp_hub2$Group,exp_hub2$SIN3A,
               levels=c("Health","Sarcopenia"))
roc_HSPD1<-roc(exp_hub2$Group,exp_hub2$HSPD1,
                levels=c("Health","Sarcopenia"))
roc_HSPA9<-roc(exp_hub2$Group,exp_hub2$HSPA9,
                  levels=c("Health","Sarcopenia"))
##SIN3A
plot(roc_HSPA8,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="HSPA8 ROC curve",
     col="#FF2E63",
     legacy.axes=T,)
##DDIT3
plot(roc_DDIT3,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="DDIT3 ROC curve",
     col="#FF2E63",
     legacy.axes=T,)
##AIFM1
plot(roc_AIFM1,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="AIFM1 ROC curve",
     col="#FF2E63",
     legacy.axes=T,)
##PLAU
plot(roc_PLAU,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="PLAU ROC curve",
     col="#FF2E63",
     legacy.axes=T,)
##RGN
plot(roc_RGN,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="RGN ROC curve",
     col="#FF2E63",
     legacy.axes=T,)
##ERCC8
plot(roc_ERCC8,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="ERCC8 ROC curve",
     col="#FF2E63",
     legacy.axes=T,)
## VEGFA
plot(roc_VEGFA,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="VEGFA ROC curve",
     col="#FF2E63",
     legacy.axes=T,)
## FOXO3
plot(roc_FOXO3,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="FOXO3 ROC curve",
     col="#FF2E63",
     legacy.axes=T,)
## FOXO1
plot(roc_FOXO1,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="FOXO1 ROC curve",
     col="#FF2E63",
     legacy.axes=T,)
## STAT5A
plot(roc_STAT5A,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="STAT5A ROC curve",
     col="#FF2E63",
     legacy.axes=T,)
## PPARGC1A
plot(roc_PPARGC1A,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="PPARGC1A ROC curve",
     col="#FF2E63",
     legacy.axes=T,)
#RPA1
plot(roc_RPA1,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="RPA1 ROC curve",
     col="#FF2E63",
     legacy.axes=T,)
#TAF1
plot(roc_TAF1,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="TAF1 ROC curve",
     col="#FF2E63",
     legacy.axes=T,)
#SIN3A
plot(roc_SIN3A,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="SIN3A ROC curve",
     col="#FF2E63",
     legacy.axes=T,)
#HSPD1
plot(roc_HSPD1,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="HSPD1 ROC curve",
     col="#FF2E63",
     legacy.axes=T,)
#HSPA9
plot(roc_HSPA9,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="HSPA9 ROC curve",
     col="#FF2E63",
     legacy.axes=T,)
## AUC>0.7的hub基因在肌肉减少症和对照组中的差异表达（箱线图）
# 大于0.7的共4个 VEGFA,FOXO3,FOXO1,STAT5A
## 把表达矩阵提取出来
exp_hub3<-exp_hub[-1,]
## 数据清洗
class(exp_hub3)
exp_hub4<-exp_hub3
exp_hub4$Symbol<-rownames(exp_hub4)
exp_hub4<-gather(exp_hub4,
                 key = sample,
                 value = expr,
                 -c("Symbol"))
## 样本分组
group_exp<-data.frame(group=c(rep('Health',8),rep('Sarcopenia',4),rep('Health',8),rep('Sarcopenia',8),
                              rep('Health',8),rep('Sarcopenia',8),rep('Health',8),rep('Sarcopenia',8),rep('Health',8),rep('Sarcopenia',4),
                              rep('Health',8),rep('Sarcopenia',4),rep('Health',8),rep('Sarcopenia',4),rep('Health',8),rep('Sarcopenia',4),
                              rep('Health',4),rep('Sarcopenia',8),rep('Health',4),rep('Sarcopenia',8),rep('Health',4),rep('Sarcopenia',8),
                              rep('Health',4),rep('Sarcopenia',8)))
colnames(group_exp)<-"Group"
exp_hub4<-cbind(group_exp,exp_hub4)
##分面图形
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
# 分面图形
exp_boxplot<-ggplot(exp_hub4,aes(x = Group, y = expr, fill = Group)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "Expression",expand = c(0.1,0.1))+
  scale_x_discrete(name = "group") +
  ggtitle("Expression of hub genes") +
  theme_bw() +
#  geom_signif(comparisons = list(c("Health","Sarcopenia")),
#              test = t.test,
#              map_signif_level = T)+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 18),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13))
exp_boxplot
exp_boxplot2<-exp_boxplot+scale_fill_lancet()+facet_wrap(~Symbol,scales = "free")
exp_boxplot2

# 08 诊断基因的外部数据集验证
setwd("/data/nas1/luchunlin/project/BJTC-167")
if (! dir.exists("./07_hub_va")){
dir.create("./07_hub_va")
}
setwd("./07_hub_va")
#将hub基因表达矩阵提取出来
# "VEGFA" "SIN3A""TAF1" "DDIT3"   "AIFM1"  "HSPD1""PLAU" "FOXO3" "FOXO1" "RPA1"  "HSPA9" "RGN"  "ERCC8"  "STAT5A"  
hubgene<-hubgene[-5,]
hubgene<-as.data.frame(hubgene)
exp_va<-exp_a
#exp_va<-read_xlsx("/data/nas1/luchunlin/project/BJTC-167/00_rawdata/GSE111010_EXP.xlsx")
#exp_va<-as.data.frame(exp_va)
#rownames(exp_va)<-exp_va[,1]
exp_hub_va<-exp_va[rownames(exp_va)%in%hubgene$hubgene,]
exp_hub_va$Symbol<-rownames(exp_hub_va)
## 数据清洗
exp_hub_va2<-gather(exp_hub_va,
                    key = sample,
                    value = expr,
                    -c("Symbol"))
## 分组
#colnames(pd_va)
#group_va<-pd_va[,c(1,50)]
#group_va<-as.data.frame(group_va)
group_va<-data.frame(group=c(rep('Health',4),rep('Sarcopenia',4),rep('Health',16),rep('Sarcopenia',4),
                             rep('Health',32),rep('Sarcopenia',4),rep('Health',52),rep('Sarcopenia',4),rep('Health',8)))
colnames(group_va)<-'Group'
exp_hub_va2<-cbind(exp_hub_va2,group_va)
# 分面图形
exp_boxplot_va<-ggplot(exp_hub_va2,aes(x = Group, y = expr, fill = Group)) +
geom_boxplot(alpha=0.7) +
scale_y_continuous(name = "Expression",expand = c(0.1,0.1))+
scale_x_discrete(name = "group") +
  ggtitle("Expression of hub genes") +
  theme_bw() +
#  geom_signif(comparisons = list(c("Health","Sarcopenia")),
#                            test = t.test,
#                            map_signif_level = T)+
  theme(plot.title = element_text(size = 16, face =  "bold"),
                      text = element_text(size = 18),
                      axis.title = element_text(face="bold"),
                      axis.text.x=element_text(size = 13))
exp_boxplot_va
exp_boxplot_va2<-exp_boxplot_va+scale_fill_lancet()+facet_wrap(~Symbol,scales = "free")
exp_boxplot_va2
## 08-2 ROC曲线
exp_hub_va2<-exp_hub_va[,-1]
exp_hub_va2<-t(exp_hub_va2)
exp_hub_va2<-cbind(design_a,exp_hub_va2)
exp_hub_va2<-as.data.frame(exp_hub_va2)
group_va2<-group_a[,2]
exp_hub_va2<-cbind(group_va2,exp_hub_va2)
exp_hub_va2<-as.data.frame(exp_hub_va2)
## 绘制ROC曲线
roc_SIN3A_va<-roc(exp_hub_va2$Group,exp_hub_va2$SIN3A,
               levels=c("Health","Sarcopenia"))
roc_DDIT3_va<-roc(exp_hub_va2$Group,exp_hub_va2$DDIT3,
                    levels=c("Health","Sarcopenia"))
roc_AIFM1_va<-roc(exp_hub_va2$Group,exp_hub_va2$AIFM1,
               levels=c("Health","Sarcopenia"))
roc_PLAU_va<-roc(exp_hub_va2$Group,exp_hub_va2$PLAU,
              levels=c("Health","Sarcopenia"))
roc_RGN_va<-roc(exp_hub_va2$Group,exp_hub_va2$RGN,
             levels=c("Health","Sarcopenia"))
roc_ERCC8_va<-roc(exp_hub_va2$Group,exp_hub_va2$ERCC8,
               levels=c("Health","Sarcopenia"))
roc_RPA1_va<-roc(exp_hub_va2$Group,exp_hub_va2$RPA1,
              levels=c("Health","Sarcopenia"))
roc_TAF1_va<-roc(exp_hub_va2$Group,exp_hub_va2$TAF1,
              levels=c("Health","Sarcopenia"),direction='<')
roc_HSPD1_va<-roc(exp_hub_va2$Group,exp_hub_va2$HSPD1,
               levels=c("Health","Sarcopenia"),direction='>')
roc_HSPA9_va<-roc(exp_hub_va2$Group,exp_hub_va2$HSPA9,
               levels=c("Health","Sarcopenia"),direction='>')
roc_VEGFA_va<-roc(exp_hub_va2$Group,exp_hub_va2$VEGFA,
                                levels=c("Health","Sarcopenia"))
roc_FOXO3_va<-roc(exp_hub_va2$Group,exp_hub_va2$FOXO3,
                                levels=c("Health","Sarcopenia"))
roc_FOXO1_va<-roc(exp_hub_va2$Group,exp_hub_va2$FOXO1,
                               levels=c("Health","Sarcopenia"))
roc_STAT5A_va<-roc(exp_hub_va2$Group,exp_hub_va2$STAT5A,
                                 levels=c("Health","Sarcopenia"))
##SIN3A
plot(roc_SIN3A_va,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="SIN3A ROC curve",
     col="#FF2E63",
     legacy.axes=T,)
##DDIT3
plot(roc_DDIT3_va,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="DDIT3 ROC curve",
     col="#FF2E63",
     legacy.axes=T,)
##AIFM1
plot(roc_AIFM1_va,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="AIFM1 ROC curve",
     col="#FF2E63",
     legacy.axes=T,)
##PLAU
plot(roc_PLAU_va,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="PLAU ROC curve",
     col="#FF2E63",
     legacy.axes=T,)
##RGN
plot(roc_RGN_va,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="RGN ROC curve",
     col="#FF2E63",
     legacy.axes=T,)
##ERCC8
plot(roc_ERCC8_va,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="ERCC8 ROC curve",
     col="#FF2E63",
     legacy.axes=T,)
## VEGFA
plot(roc_VEGFA_va,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="VEGFA ROC curve",
     col="#FF2E63",
      legacy.axes=T,)
## FOXO3
plot(roc_FOXO3_va,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="FOXO3 ROC curve",
     col="#FF2E63",
     legacy.axes=T,) 
## FOXO1
plot(roc_FOXO1_va,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="FOXO1 ROC curve",
     col="#FF2E63",
     legacy.axes=T,) 
## STAT5A
plot(roc_STAT5A_va,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="STAT5A ROC curve",
     col="#FF2E63",
     legacy.axes=T,) 

#RPA1
plot(roc_RPA1_va,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="RPA1 ROC curve",
     col="#FF2E63",
     legacy.axes=T,) 
#TAF1
plot(roc_TAF1_va,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="TAF1 ROC curve",
     col="#FF2E63",
     legacy.axes=T,) 
#SIN3A
plot(roc_SIN3A_va,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="SIN3A ROC curve",
     col="#FF2E63",
     legacy.axes=T,) 
#HSPD1
plot(roc_HSPD1_va,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="HSPD1 ROC curve",
     col="#FF2E63",
     legacy.axes=T,) 
#HSPA9
plot(roc_HSPA9_va,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="HSPA9 ROC curve",
     col="#FF2E63",
     legacy.axes=T,) 
# 09 hub基因的表达相关性分析
setwd("/data/nas1/luchunlin/project/BJTC-167")
if (! dir.exists("./08_correlation")){
dir.create("./08_correlation")
}
setwd("./08_correlation")
# 筛选到的5个hub基因：FOXO1,STAT5A,VEGFA,FOXO3
library(ggcorrplot)
library(corrplot)
hub_gene_final<-c("FOXO3","VEGFA","FOXO1","STAT5A","PPARGC1A")
hub_gene_final<-as.data.frame(hub_gene_final)
#提取5个诊断基因的表达矩阵
exp_final<-exp_b[rownames(exp_b)%in%hub_gene_final$hub_gene_final,]
hub_corr<-round(cor(t(exp_final)),3)
## 检验基因之间的相关性p值
## 计算相关性系数并显示基因之间的相关性。相关性系数大于0为正相关，小于0为负相关。
## p小于0.05认为相关性显著
hub_p.mat<-round(cor_pmat(t(exp_final)),3)
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
## 10 GSEA
setwd("/data/nas1/luchunlin/project/BJTC-167")
if (! dir.exists("./09_GSEA")){
  dir.create("./09_GSEA")
}
setwd("./09_GSEA")
## 根据诊断基因表达水平的中位值，将肌肉减少症样本高低分组，鉴定差异基因后，GSEA，分析各个诊断基因与哪些信号通路有关。
## （参考基因集为MSigDB数据库中的C2：KEGG gene sets）
GSEA_exp<-exp_b
# 将肌肉减少症样本提取出来
GSEA_exp<-t(GSEA_exp)
GSEA_exp<-cbind(group_b,GSEA_exp)
GSEA_exp<-GSEA_exp[order(GSEA_exp$Group),]
GSEA_exp<-as.data.frame(GSEA_exp)
rownames(GSEA_exp)<-GSEA_exp[,1]
GSEA_exp<-GSEA_exp[-c(1:20),]
GSEA_exp<-as.data.frame(GSEA_exp)
GSEA_exp<-GSEA_exp[,-c(1,2)]
GSEA_exp<-t(GSEA_exp)
## 按照不同的诊断基因高低分组（中位值），鉴定差异基因，GSEA    FOXO1,STAT5A,RGN,SIN3A
##10-1 PPARGC1A
## 分组
PPARGC1A_exp<-t(GSEA_exp)
PPARGC1A_exp<-as.data.frame(PPARGC1A_exp)
PPARGC1A_exp<-PPARGC1A_exp[order(PPARGC1A_exp$PPARGC1A),]
PPARGC1A_exp2<-t(PPARGC1A_exp)
median(PPARGC1A_exp$PPARGC1A)
PPARGC1A_exp$PPARGC1A
group_PPARGC1A<-data.frame(sample=rownames(PPARGC1A_exp),
                        group=c(rep("Low",9),rep("High",10)))
rownames(group_PPARGC1A)<-rownames(PPARGC1A_exp)
design_PPARGC1A<-as.matrix(group_PPARGC1A)
type_PPARGC1A<-design_PPARGC1A[,2]
design_PPARGC1A<-model.matrix(~-1+factor(type_PPARGC1A,levels = c('Low','High')))
colnames(design_PPARGC1A)<-c('Low','High')
## 差异分析
library(limma)
fit_PPARGC1A=lmFit(PPARGC1A_exp2,design_PPARGC1A)
# 构建比较矩阵
contrast.matrix_PPARGC1A=makeContrasts(ControlVSMG=Low-High,levels = design_PPARGC1A)
#构建芯片数据线性模型，计算估计的相关系数和标准差
fit2_PPARGC1A=contrasts.fit(fit_PPARGC1A,contrast.matrix_PPARGC1A)
# 基于贝叶斯计算T值，F值和log-odds
fit2_PPARGC1A<-eBayes(fit2_PPARGC1A,0.01)
DEGs_PPARGC1A=topTable(fit2_PPARGC1A,coef =1,n=Inf)
DEGs_PPARGC1A=na.omit(DEGs_PPARGC1A)
head(DEGs_PPARGC1A)
write.table(DEGs_PPARGC1A,file = "DEG_PPARGC1A.xls",
            quote = F,
            sep = "\t",
            row.names = T)
# 筛选差异基因
logFC_cutoff <- 0
DEGs_PPARGC1A$change = as.factor(
  ifelse(DEGs_PPARGC1A$P.Value < 0.05 & abs(DEGs_PPARGC1A$logFC) > logFC_cutoff,
         ifelse(DEGs_PPARGC1A$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
sig_diff_PPARGC1A <- subset(DEGs_PPARGC1A,
                         DEGs_PPARGC1A$P.Value < 0.05& abs(DEGs_PPARGC1A$logFC) >= 0)
dim(DEGs_PPARGC1A)
dim(sig_diff_PPARGC1A)
summary(sig_diff_PPARGC1A$change)
## GSEA
##需要包含基因名和logFC的数据文件。并把基因名转化为ENTREZID
GSEA_PPARGC1A<-data.frame(sig_diff_PPARGC1A[,1])
GSEA_PPARGC1A$Symbol<-rownames(sig_diff_PPARGC1A)
rownames(GSEA_PPARGC1A)<-GSEA_PPARGC1A$Symbol
colnames(GSEA_PPARGC1A)<-c('logFC','Symbol')
PPARGC1A_gene<-rownames(GSEA_PPARGC1A)
gene_transform_PPARGC1A<-bitr(PPARGC1A_gene,
                           fromType = "SYMBOL",
                           toType = c("ENTREZID"),
                           OrgDb = "org.Hs.eg.db")
GSEA_PPARGC1A$Symbol<-gene_transform_PPARGC1A$ENTREZID
colnames(GSEA_PPARGC1A)<-c('logFC','Entrezid')
## 按照logFC降序排列
GSEA_PPARGC1A<-GSEA_PPARGC1A[order(GSEA_PPARGC1A$logFC,decreasing = T),]
## 整理成GSEA分析格式
genelist_PPARGC1A<-GSEA_PPARGC1A[,1]
names(genelist_PPARGC1A)=as.character(GSEA_PPARGC1A[,2])
genelist_PPARGC1A
kegg_set<-read.gmt("c2.cp.kegg.v7.4.entrez.gmt")
kegg_PPARGC1A<-GSEA(genelist_PPARGC1A,TERM2GENE = kegg_set,pvalueCutoff = 1)
kegg_result_PPARGC1A<-kegg_PPARGC1A@result
gseaplot2(kegg_PPARGC1A,1:5,color = 'red')
write.table(kegg_result_PPARGC1A,file = "PPARGC1A_GSEA.xls",sep = "\t",quote = F,row.names = F)
##10-2 FOXO3
## 分组
FOXO3_exp<-t(GSEA_exp)
FOXO3_exp<-as.data.frame(FOXO3_exp)
FOXO3_exp<-FOXO3_exp[order(FOXO3_exp$FOXO3),]
FOXO3_exp2<-t(FOXO3_exp)
median(FOXO3_exp$FOXO3)
FOXO3_exp$FOXO3
group_FOXO3<-data.frame(sample=rownames(FOXO3_exp),
                      group=c(rep("Low",9),rep("High",10)))
rownames(group_FOXO3)<-rownames(FOXO3_exp)
design_FOXO3<-as.matrix(group_FOXO3)
type_FOXO3<-design_FOXO3[,2]
design_FOXO3<-model.matrix(~-1+factor(type_FOXO3,levels = c('Low','High')))
colnames(design_FOXO3)<-c('Low','High')
## 差异分析
library(limma)
fit_FOXO3=lmFit(FOXO3_exp2,design_FOXO3)
# 构建比较矩阵
contrast.matrix_FOXO3=makeContrasts(ControlVSMG=Low-High,levels = design_FOXO3)
#构建芯片数据线性模型，计算估计的相关系数和标准差
fit2_FOXO3=contrasts.fit(fit_FOXO3,contrast.matrix_FOXO3)
# 基于贝叶斯计算T值，F值和log-odds
fit2_FOXO3<-eBayes(fit2_FOXO3,0.01)
DEGs_FOXO3=topTable(fit2_FOXO3,coef =1,n=Inf)
DEGs_FOXO3=na.omit(DEGs_FOXO3)
head(DEGs_FOXO3)
write.table(DEGs_FOXO3,file = "DEG_FOXO3.xls",
            quote = F,
            sep = "\t",
            row.names = T)
# 筛选差异基因
logFC_cutoff <- 0
DEGs_FOXO3$change = as.factor(
  ifelse(DEGs_FOXO3$P.Value < 0.05 & abs(DEGs_FOXO3$logFC) > logFC_cutoff,
         ifelse(DEGs_FOXO3$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
sig_diff_FOXO3 <- subset(DEGs_FOXO3,
                       DEGs_FOXO3$P.Value < 0.05& abs(DEGs_FOXO3$logFC) >= 0)
dim(DEGs_FOXO3)
dim(sig_diff_FOXO3)
summary(sig_diff_FOXO3$change)
## GSEA
##需要包含基因名和logFC的数据文件。并把基因名转化为ENTREZID
GSEA_FOXO3<-data.frame(sig_diff_FOXO3[,1])
GSEA_FOXO3$Symbol<-rownames(sig_diff_FOXO3)
rownames(GSEA_FOXO3)<-GSEA_FOXO3$Symbol
colnames(GSEA_FOXO3)<-c('logFC','Symbol')
FOXO3_gene<-rownames(GSEA_FOXO3)
gene_transform_FOXO3<-bitr(FOXO3_gene,
                         fromType = "SYMBOL",
                         toType = c("ENTREZID"),
                         OrgDb = "org.Hs.eg.db")
GSEA_FOXO3$Symbol<-gene_transform_FOXO3$ENTREZID
colnames(GSEA_FOXO3)<-c('logFC','Entrezid')
## 按照logFC降序排列
GSEA_FOXO3<-GSEA_FOXO3[order(GSEA_FOXO3$logFC,decreasing = T),]
## 整理成GSEA分析格式
genelist_FOXO3<-GSEA_FOXO3[,1]
names(genelist_FOXO3)=as.character(GSEA_FOXO3[,2])
genelist_FOXO3
kegg_set<-read.gmt("c2.cp.kegg.v7.4.entrez.gmt")
kegg_FOXO3<-GSEA(genelist_FOXO3,TERM2GENE = kegg_set,pvalueCutoff = 1)
kegg_result_FOXO3<-kegg_FOXO3@result
gseaplot2(kegg_FOXO3,1:5,color = 'red')
write.table(kegg_result_FOXO3,file = 'FOXO3_GSEA.xls',
            quote = F,
            sep = '\t')

##10-3 FOXO1
FOXO1_exp<-t(GSEA_exp)
FOXO1_exp<-as.data.frame(FOXO1_exp)
FOXO1_exp<-FOXO1_exp[order(FOXO1_exp$FOXO1),]
FOXO1_exp2<-t(FOXO1_exp)
median(FOXO1_exp$FOXO1)
FOXO1_exp$FOXO1
group_FOXO1<-data.frame(sample=rownames(FOXO1_exp),
                        group=c(rep("Low",9),rep("High",10)))
rownames(group_FOXO1)<-rownames(FOXO1_exp)
design_FOXO1<-as.matrix(group_FOXO1)
type_FOXO1<-design_FOXO1[,2]
design_FOXO1<-model.matrix(~-1+factor(type_FOXO1,levels = c('Low','High')))
colnames(design_FOXO1)<-c('Low','High')
## 差异分析
library(limma)
fit_FOXO1=lmFit(FOXO1_exp2,design_FOXO1)
# 构建比较矩阵
contrast.matrix_FOXO1=makeContrasts(ControlVSMG=Low-High,levels = design_FOXO1)
#构建芯片数据线性模型，计算估计的相关系数和标准差
fit2_FOXO1=contrasts.fit(fit_FOXO1,contrast.matrix_FOXO1)
# 基于贝叶斯计算T值，F值和log-odds
fit2_FOXO1<-eBayes(fit2_FOXO1,0.01)
DEGs_FOXO1=topTable(fit2_FOXO1,coef =1,n=Inf)
DEGs_FOXO1=na.omit(DEGs_FOXO1)
head(DEGs_FOXO1)
write.table(DEGs_FOXO1,file = "DEG_FOXO1.xls",
            quote = F,
            sep = "\t",
            row.names = T)
# 筛选差异基因
logFC_cutoff <- 0
DEGs_FOXO1$change = as.factor(
  ifelse(DEGs_FOXO1$P.Value < 0.05 & abs(DEGs_FOXO1$logFC) > logFC_cutoff,
         ifelse(DEGs_FOXO1$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
sig_diff_FOXO1 <- subset(DEGs_FOXO1,
                         DEGs_FOXO1$P.Value < 0.05& abs(DEGs_FOXO1$logFC) >= 0)
dim(DEGs_FOXO1)
dim(sig_diff_FOXO1)
summary(sig_diff_FOXO1$change)

## GSEA
##需要包含基因名和logFC的数据文件。并把基因名转化为ENTREZID
GSEA_FOXO1<-data.frame(sig_diff_FOXO1[,1])
GSEA_FOXO1$Symbol<-rownames(sig_diff_FOXO1)
rownames(GSEA_FOXO1)<-GSEA_FOXO1$Symbol
colnames(GSEA_FOXO1)<-c('logFC','Symbol')
FOXO1_gene<-rownames(GSEA_FOXO1)
gene_transform_FOXO1<-bitr(FOXO1_gene,
                           fromType = "SYMBOL",
                           toType = c("ENTREZID"),
                           OrgDb = "org.Hs.eg.db")
GSEA_FOXO1$Symbol<-gene_transform_FOXO1$ENTREZID
colnames(GSEA_FOXO1)<-c('logFC','Entrezid')
## 按照logFC降序排列
GSEA_FOXO1<-GSEA_FOXO1[order(GSEA_FOXO1$logFC,decreasing = T),]
## 整理成GSEA分析格式
genelist_FOXO1<-GSEA_FOXO1[,1]
names(genelist_FOXO1)=as.character(GSEA_FOXO1[,2])
genelist_FOXO1
kegg_FOXO1<-GSEA(genelist_FOXO1,TERM2GENE = kegg_set,pvalueCutoff = 0.5)
kegg_result_FOXO1<-kegg_FOXO1@result
gseaplot2(kegg_FOXO1,1:5,color = 'red')
write.table(kegg_result_FOXO1,file = 'FOXO1_GSEA.xls',
            quote = F,
            sep = '\t')
##10-4 STAT5A
STAT5A_exp<-t(GSEA_exp)
STAT5A_exp<-as.data.frame(STAT5A_exp)
STAT5A_exp<-STAT5A_exp[order(STAT5A_exp$STAT5A),]
STAT5A_exp2<-t(STAT5A_exp)
median(STAT5A_exp$STAT5A)
STAT5A_exp$STAT5A
group_STAT5A<-data.frame(sample=rownames(STAT5A_exp),
                         group=c(rep("Low",9),rep("High",10)))
rownames(group_STAT5A)<-rownames(STAT5A_exp)
design_STAT5A<-as.matrix(group_STAT5A)
type_STAT5A<-design_STAT5A[,2]
design_STAT5A<-model.matrix(~-1+factor(type_STAT5A,levels = c('Low','High')))
colnames(design_STAT5A)<-c('Low','High')
## 差异分析
library(limma)
fit_STAT5A=lmFit(STAT5A_exp2,design_STAT5A)
# 构建比较矩阵
contrast.matrix_STAT5A=makeContrasts(ControlVSMG=Low-High,levels = design_STAT5A)
#构建芯片数据线性模型，计算估计的相关系数和标准差
fit2_STAT5A=contrasts.fit(fit_STAT5A,contrast.matrix_STAT5A)
# 基于贝叶斯计算T值，F值和log-odds
fit2_STAT5A<-eBayes(fit2_STAT5A,0.01)
DEGs_STAT5A=topTable(fit2_STAT5A,coef =1,n=Inf)
DEGs_STAT5A=na.omit(DEGs_STAT5A)
head(DEGs_STAT5A)
write.table(DEGs_STAT5A,file = "DEG_STAT5A.xls",
            quote = F,
            sep = "\t",
            row.names = T)
# 筛选差异基因
logFC_cutoff <- 0
DEGs_STAT5A$change = as.factor(
  ifelse(DEGs_STAT5A$P.Value < 0.05 & abs(DEGs_STAT5A$logFC) > logFC_cutoff,
         ifelse(DEGs_STAT5A$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
sig_diff_STAT5A <- subset(DEGs_STAT5A,
                          DEGs_STAT5A$P.Value < 0.05& abs(DEGs_STAT5A$logFC) >= 0)
dim(DEGs_STAT5A)
dim(sig_diff_STAT5A)
summary(sig_diff_STAT5A$change)
# DOWN  NOT   UP    ## 2466
#  1380    0 1086
## GSEA
##需要包含基因名和logFC的数据文件。并把基因名转化为ENTREZID
GSEA_STAT5A<-data.frame(sig_diff_STAT5A[,1])
GSEA_STAT5A$Symbol<-rownames(sig_diff_STAT5A)
rownames(GSEA_STAT5A)<-GSEA_STAT5A$Symbol
colnames(GSEA_STAT5A)<-c('logFC','Symbol')
STAT5A_gene<-rownames(GSEA_STAT5A)
gene_transform_STAT5A<-bitr(STAT5A_gene,
                            fromType = "SYMBOL",
                            toType = c("ENTREZID"),
                            OrgDb = "org.Hs.eg.db")
GSEA_STAT5A$Symbol<-gene_transform_STAT5A$ENTREZID
colnames(GSEA_STAT5A)<-c('logFC','Entrezid')
## 按照logFC降序排列
GSEA_STAT5A<-GSEA_STAT5A[order(GSEA_STAT5A$logFC,decreasing = T),]
## 整理成GSEA分析格式
genelist_STAT5A<-GSEA_STAT5A[,1]
names(genelist_STAT5A)=as.character(GSEA_STAT5A[,2])
genelist_STAT5A
kegg_STAT5A<-GSEA(genelist_STAT5A,TERM2GENE = kegg_set,pvalueCutoff = 1)
kegg_result_STAT5A<-kegg_STAT5A@result
gseaplot2(kegg_STAT5A,1:5,color = 'red')
write.table(kegg_result_STAT5A,file = 'STAT5A_GSEA.xls',
            quote = F,
            sep = '\t')
##10-5 VEGFA
VEGFA_exp<-t(GSEA_exp)
VEGFA_exp<-as.data.frame(VEGFA_exp)
VEGFA_exp<-VEGFA_exp[order(VEGFA_exp$VEGFA),]
VEGFA_exp2<-t(VEGFA_exp)
median(VEGFA_exp$VEGFA)
VEGFA_exp$VEGFA
group_VEGFA<-data.frame(sample=rownames(VEGFA_exp),
                         group=c(rep("Low",9),rep("High",10)))
rownames(group_VEGFA)<-rownames(VEGFA_exp)
design_VEGFA<-as.matrix(group_VEGFA)
type_VEGFA<-design_VEGFA[,2]
design_VEGFA<-model.matrix(~-1+factor(type_VEGFA,levels = c('Low','High')))
colnames(design_VEGFA)<-c('Low','High')
## 差异分析
library(limma)
fit_VEGFA=lmFit(VEGFA_exp2,design_VEGFA)
# 构建比较矩阵
contrast.matrix_VEGFA=makeContrasts(ControlVSMG=Low-High,levels = design_VEGFA)
#构建芯片数据线性模型，计算估计的相关系数和标准差
fit2_VEGFA=contrasts.fit(fit_VEGFA,contrast.matrix_VEGFA)
# 基于贝叶斯计算T值，F值和log-odds
fit2_VEGFA<-eBayes(fit2_VEGFA,0.01)
DEGs_VEGFA=topTable(fit2_VEGFA,coef =1,n=Inf)
DEGs_VEGFA=na.omit(DEGs_VEGFA)
head(DEGs_VEGFA)
write.table(DEGs_VEGFA,file = "DEG_VEGFA.xls",
            quote = F,
            sep = "\t",
            row.names = T)
# 筛选差异基因
logFC_cutoff <- 0
DEGs_VEGFA$change = as.factor(
  ifelse(DEGs_VEGFA$P.Value < 0.05 & abs(DEGs_VEGFA$logFC) > logFC_cutoff,
         ifelse(DEGs_VEGFA$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
sig_diff_VEGFA <- subset(DEGs_VEGFA,
                          DEGs_VEGFA$P.Value < 0.05& abs(DEGs_VEGFA$logFC) >= 0)
dim(DEGs_VEGFA)
dim(sig_diff_VEGFA)
summary(sig_diff_VEGFA$change)
## GSEA
##需要包含基因名和logFC的数据文件。并把基因名转化为ENTREZID
GSEA_VEGFA<-data.frame(sig_diff_VEGFA[,1])
GSEA_VEGFA$Symbol<-rownames(sig_diff_VEGFA)
rownames(GSEA_VEGFA)<-GSEA_VEGFA$Symbol
colnames(GSEA_VEGFA)<-c('logFC','Symbol')
VEGFA_gene<-rownames(GSEA_VEGFA)
gene_transform_VEGFA<-bitr(VEGFA_gene,
                            fromType = "SYMBOL",
                            toType = c("ENTREZID"),
                            OrgDb = "org.Hs.eg.db")
GSEA_VEGFA$Symbol<-gene_transform_VEGFA$ENTREZID
colnames(GSEA_VEGFA)<-c('logFC','Entrezid')
## 按照logFC降序排列
GSEA_VEGFA<-GSEA_VEGFA[order(GSEA_VEGFA$logFC,decreasing = T),]
## 整理成GSEA分析格式
genelist_VEGFA<-GSEA_VEGFA[,1]
names(genelist_VEGFA)=as.character(GSEA_VEGFA[,2])
genelist_VEGFA
kegg_VEGFA<-GSEA(genelist_VEGFA,TERM2GENE = kegg_set,pvalueCutoff = 1)
kegg_result_VEGFA<-kegg_VEGFA@result
gseaplot2(kegg_VEGFA,1:4,color = 'red')
write.table(kegg_result_VEGFA,file = 'VEGFA_GSEA.xls',
            quote = F,
            sep = '\t')
## 11 诊断基因功能相似性分析（GOSemSim包）
setwd("/data/nas1/luchunlin/project/BJTC-167")
if (! dir.exists("./10_relationship")){
dir.create("./10_relationship")
}
setwd("./10_relationship")
library(GOSemSim)
library(org.Hs.eg.db)
BiocManager::install("org.Hs.eg.db")
bp<-godata('org.Hs.eg.db',ont = "BP",computeIC = F)
cc<-godata('org.Hs.eg.db',ont = "CC",computeIC = F)
mf<-godata('org.Hs.eg.db',ont = "MF",computeIC = F)
geneid2<-hub_gene_final$hub_gene_final
gene_transform2<-bitr(geneid2,fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = org.Hs.eg.db)
genelist<-gene_transform2
genelist<-hub_gene_final
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
colnames(fsim)<-genelist$SYMBOL
rownames(fsim)<-genelist$SYMBOL
## 进一步去除基因与本身基因之间的相关性，使用melt函数将宽格式数据转化为长格式数据
for (i in 1:ncol(fsim)) {
fsim[i,i] <- NA
}
dat<-melt(fsim)
View(dat)
dat<-dat[!is.na(dat$value),]
dat<-dat[,c(1,3)]
head(dat)
## 绘制boxplot图
dat.mean<-aggregate(value~Var1,dat,mean)
View(dat.mean)
## 根据相似性评分的平均值，对其进行排序，并根据评分的高低，将基因名设置为因子（factor）格式
m<-dat.mean$value
names(m)<-dat.mean$Var1
dat$Var1<-factor(dat$Var1,
                 levels = names(sort(m)))
str(dat)

## 使用ggplot 对分析结果进行可视化。
ggplot(dat,
       aes(x=Var1,y=value,fill=factor(Var1)))+
  scale_fill_brewer(palette = "Set2")+
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
