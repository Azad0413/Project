rm(list = ls())
# 01 获取数据集--------------
setwd("/data/nas1/luchunlin/project/SJZZK-431-10//")
if (! dir.exists("./02_WGCNA")){
  dir.create("./02_WGCNA")
}
setwd("./02_WGCNA")
library(tidyverse)
library(lance)
library(GSVA)
library(readxl)
anikios <- read_xlsx('anikios.xlsx')
anikios$anikios <- c(rep('anikios',35))

dat<-read.delim2('../00_rawdata/dat.xls',row.names = 1)%>%lc.tableToNum()
group<-read_xlsx('../00_rawdata/group.xlsx')
table(group$group)

gene_list <- split(as.matrix(anikios)[,1],
                   anikios[,2])
dat2 <- as.matrix(dat)
anikios_score = gsva(dat2, gene_list,
                    method = "gsva",
                    ssgsea.norm = TRUE,
                    verbose = TRUE)
write.table(anikios_score,
            file = "anikios_score.xls",
            sep = "\t",
            quote = F)
anikios_score<-t(anikios_score)%>%as.data.frame()



library(WGCNA)
library(reshape2)
library(stringr)
options(stringsAsFactors = F)
enableWGCNAThreads()
exprMat<-dat
dim(exprMat)
type = "unsigned"
corType = "pearson"
corFnc = ifelse(corType=="pearson", WGCNA::cor, bicor)
# 对二元变量，如样本性状信息计算相关性时，
# 或基因表达严重依赖于疾病状态时，需设置下面参数
maxPOutliers = ifelse(corType=="pearson",1,0.05)
# 关联样品性状的二元变量时，设置
robustY = ifelse(corType=="pearson",T,F)
dataExpr <- exprMat[rownames(dat),]
## 04-1 数据筛选-----
## 筛选中位绝对偏差（MAD）前75%的基因，至少MAD大于0.01
## 筛选后会降低运算量，也会失去部分信息
## 也可不做筛选，使MAD大于0即可
m.mad <- apply(dataExpr,1,mad)
dataExprVar <- dataExpr[which(m.mad >
                                max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
#dataExprVar <- dataExpr
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
# [1] 109 18879

## 04-2 软阈值筛选----

## 样本聚类，查看是否有利群样本

sampleTree = hclust(dist(dataExpr,method = 'euclidean'), method = "complete")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",label=F)
pdf('01.sample_tree.pdf',w=8,h=5)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",label=F)
#abline(h =70, col = "red")
dev.off()
png('01.sample_tree.png',w=700,h=500)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="",label=F)
#abline(h =70, col = "red")
dev.off()
# Determine cluster under the line
# 剪枝算法，cutHeight修剪树枝的高度，minSize集群最小数
#clust = cutreeStatic(sampleTree, cutHeight = 70, minSize = 10)
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
powers = c(c(1:20))
sft = pickSoftThreshold(dataExpr, powerVector=powers,
                        networkType=type, verbose=5)
par(mfrow = c(1,2))
cex1 = 0.9
# 横轴是Soft threshold (power)，纵轴是无标度网络的评估参数，数值越高，
# 网络越符合无标度特征 (non-scale)
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",
     ylab="Scale Free Topology Model Fitsigned R^2",type="n",
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
## 3
## 04-3一步法网络构建---------
## One-step network construction and module detection##
# power: 上一步计算的软阈值
# cor <- WGCNA::cor
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
# 0    1    2    3    4    5    6    7   
# 3048 5043 4585 3821 1268  620  282  212 
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
png(filename = "03.cluster_dendrogram.png", height = 600, width = 800)
plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
pdf(file = "03.cluster_dendrogram.pdf", height = 6, width = 8)
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
# sizeGrWindow(6,6)
# pdf('04.Eigengene_dendrogram_and_heatmap.pdf',w=5,h=5)
# plotEigengeneNetworks(MEs_col, 
#                       setLabels = "Eigengene dendrogram and heatmap", 
#                       marDendro = c(1,3,4,4),
#                       marHeatmap = c(3,3,0,2), 
#                       plotDendrograms = T, 
#                       xLabelsAngle = 90)
# 
# dev.off()
# png('04.Eigengene_dendrogram_and_heatmap.png',w=500,h=500)
# plotEigengeneNetworks(MEs_col, 
#                       setLabels = "Eigengene dendrogram and heatmap", 
#                       marDendro = c(1,3,4,4),
#                       marHeatmap = c(3,3,0,2), 
#                       plotDendrograms = T, 
#                       xLabelsAngle = 90)
# 
# dev.off()
# Plot the dendrogram
plotEigengeneNetworks(MEs_col, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)


## 04-6 关联表型数据------
group_traits<-anikios_score

datTraits=data.frame(samples=rownames(dataExpr),subtype=group_traits)

moduleColors <- labels2colors(net$colors)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(dataExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0); ##不同颜色的模块的ME值矩阵(样本vs模块)
moduleTraitCor = cor(MEs, group_traits , use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
# moduleTraitCor <- moduleTraitCor[-8,]
# moduleTraitPvalue <- moduleTraitPvalue[-8,]
# MEs <- MEs[,-8]
sizeGrWindow(10,10)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))

labeledHeatmap(Matrix = as.matrix(moduleTraitCor[-8,]),
               xLabels = names(group_traits),
               yLabels = names(MEs[,-8]),
               ySymbols = names(MEs[,-8]),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix[-8,],
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = paste("Module-Anikios relationships"))

### 模块与基因相关性
if (corType=="pearson") {
  geneModuleMembership = as.data.frame(cor(dataExpr, MEs_col, use = "p"))
  MMPvalue = as.data.frame(corPvalueStudent(
    as.matrix(geneModuleMembership), nSamples))
} else {
  geneModuleMembershipA = bicorAndPvalue(dataExpr, MEs_col, robustY=robustY)
  geneModuleMembership = geneModuleMembershipA$bicor
  MMPvalue   = geneModuleMembershipA$p
}
# 计算性状与基因的相关性矩阵

## 只有连续型性状才能进行计算，如果是离散变量，在构建样品表时就转为0-1矩阵。

if (corType=="pearson") {
  geneTraitCor = as.data.frame(cor(dataExpr, group_traits, use = "p"))
  geneTraitP = as.data.frame(corPvalueStudent(
    as.matrix(geneTraitCor), nSamples))
} else {
  geneTraitCorA = bicorAndPvalue(dataExpr, group_traits, robustY=robustY)
  geneTraitCor = as.data.frame(geneTraitCorA$bicor)
  geneTraitP   = as.data.frame(geneTraitCorA$p)
}

### black-----
module = "black"
probes=colnames(dataExpr)
inModule=(moduleColors==module)
modProbes=probes[inModule]
modGenes<-as.data.frame(modProbes)
colnames(modGenes)<-'modgene'


pheno = colnames(group_traits)
modNames = substring(colnames(MEs_col), 3)
# 获取关注的列
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(group_traits))
# 获取模块内的基因
moduleGenes = moduleColors == module
MM<-abs(geneModuleMembership[moduleGenes,module_column])
GS<-abs(geneTraitCor[moduleGenes, 1])
c<-as.data.frame(cbind(MM,GS))
rownames(c)=modGenes$modgene
black_hub <- subset(c,c$MM > 0.5 & c$GS > 0.2)
##490
sizeGrWindow(7, 7)
par(mfrow = c(1,1),
    mai = c(1,1,1,0.5))
# 与性状高度相关的基因，也是与性状相关的模型的关键基因
pdf('05.black_cor.pdf',w=7,h=6)
verboseScatterplot(abs(geneModuleMembership[moduleGenes,module_column]),
                   abs(geneTraitCor[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for DEimmune Cell",
                   main = paste("Module membership vs. gene significance"), 
                   pch = 20,col=module,
                   cex.lab = 1.2,cex.axis = 1.2,cex.main = 1.2
)
abline(h=0.2,v=0.5,col="red",lwd=1.5)
dev.off()
png('05.black_cor.png',w=500,h=400)
verboseScatterplot(abs(geneModuleMembership[moduleGenes,module_column]),
                   abs(geneTraitCor[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for DEimmune Cell",
                   main = paste("Module membership vs. gene significance"), 
                   pch = 20,col=module,
                   cex.lab = 1.2,cex.axis = 1.2,cex.main = 1.2
)
abline(h=0.2,v=0.5,col="red",lwd=1.5)
dev.off()

## MElightgreen---------
module='turquoise'
probes=colnames(dataExpr)
inModule=(moduleColors==module)
modProbes=probes[inModule]
modGenes<-as.data.frame(modProbes)
colnames(modGenes)<-'modgene'
##66
pheno = colnames(group_traits)
modNames = substring(colnames(MEs_col), 3)
# 获取关注的列
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(group_traits))
# 获取模块内的基因
moduleGenes = moduleColors == module
MM<-abs(geneModuleMembership[moduleGenes,module_column])
GS<-abs(geneTraitCor[moduleGenes, 1])
c<-as.data.frame(cbind(MM,GS))
rownames(c)=modGenes$modgene
turquoise_hub <- subset(c,c$MM > 0.5 & c$GS > 0.2)
##3
sizeGrWindow(7, 7)
par(mfrow = c(1,1),
    mai = c(1,1,1,0.5))
# 与性状高度相关的基因，也是与性状相关的模型的关键基因
pdf('06.turquoise_cor.pdf',w=7,h=6)
verboseScatterplot(abs(geneModuleMembership[moduleGenes,module_column]),
                   abs(geneTraitCor[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for DEimmune Cell",
                   main = paste("Module membership vs. gene significance"), 
                   pch = 20,col=module,
                   cex.lab = 1.2,cex.axis = 1.2,cex.main = 1.2
)
abline(h=0.2,v=0.5,col="red",lwd=1.5)
dev.off()
png('06.turquoise_cor.png',w=500,h=400)
verboseScatterplot(abs(geneModuleMembership[moduleGenes,module_column]),
                   abs(geneTraitCor[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for DEimmune Cell",
                   main = paste("Module membership vs. gene significance"), 
                   pch = 20,col=module,
                   cex.lab = 1.2,cex.axis = 1.2,cex.main = 1.2
)
abline(h=0.2,v=0.5,col="red",lwd=1.5)
dev.off()

modgene <- data.frame(rownames(rbind(black_hub,turquoise_hub)))
colnames(modgene) <- 'modgene'
## 2678
write.table(modgene,file = 'modgene.xls',sep = '\t',row.names = F,quote = F)


