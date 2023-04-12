# rm(list = ls())
# 01 获取数据集--------------
setwd("/data/nas1/luchunlin/project/BJTC-251")
if (! dir.exists("./00_rawdata")){
  dir.create("./00_rawdata")
}
setwd("./00_rawdata")
library(GEOquery)
library(Biobase)
library(limma)
## GSE27276
gset<-getGEO("GSE27276",
              destdir = '.',
              GSEMatrix = T,
              getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL2507",destdir = '.')
a=gset[[1]]
pd<-pData(a)
View(pd)
gpl1<-Table(gpl)                  
colnames(Table(gpl))         
View(gpl1)  
## 没有Symbol名，先换成GB_ACC     
probe2GB<-gpl1[,c(1,3)]           
colnames(probe2GB)=c('probe_id','GB_ACC')
length(unique(probe2GB$GB_ACC))
# [1] 25360
table(sort(table(probe2GB$GB_ACC)))
# 去掉没有注释的探针
ids=probe2GB[probe2GB$GB_ACC!='',]
#判断是否匹配
ids=probe2GB[probe2GB$probe_id%in%rownames(expr),]
dat=expr[ids$probe_id,]
ids$mean=apply(dat,1,mean)
ids=ids[order(ids$GB_ACC,ids$mean,decreasing = T),]
ids<-separate(ids,GB_ACC,into = c('GB_ACC'),sep = '\\.')
ids=ids[!duplicated(ids$GB_ACC),]
dat=dat[ids$probe_id,]
rownames(dat)=ids$GB_ACC
dat=dat[-25354,]

## 然后将GB_ACC换成Symbol
library(org.Hs.eg.db)
library(clusterProfiler)
geneid<-rownames(dat)
gene_transform<-bitr(geneid,fromType = "REFSEQ",
                     toType = "SYMBOL",
                     OrgDb = org.Hs.eg.db)
dat2<-dat
colnames(gene_transform)=c('REFSEQ','Symbol')
gene_transform=gene_transform[!duplicated(gene_transform$Symbol),]
dat2=dat2[gene_transform$REFSEQ,]
rownames(dat2)=gene_transform$Symbol
## 16351
write.table(dat2,file = "EXP.xls",
            quote = F,
            sep = '\t',
            row.names = T)

## 验证集GSE138125
gset_va<-getGEO("GSE138125",
               destdir = '.',
               GSEMatrix = T,
               getGPL = F)
va=gset_va[[1]]
pd_va<-pData(va)
exp_va<-exprs(va)
##4个疾病 4个对照
## gpl_va<-getGEO("GPL21827",destdir = '.')
## gpl_va1<-Table(gpl_va)
### 发现GEO上芯片信息只提供了探针序列，无法进行转化。
### 使用idmap3包转化
#library(devtools)
#install_github("jmzeng1314/idmap3")
#library(idmap3)
#ids2=idmap3::get_pipe_IDs('GPL21827')
#head(ids2) 
#length(unique(ids2$probe_id))
#ids2=ids2[ids2$probe_id!='',]
#ids2=ids2[ids2$probe_id%in%rownames(exp_va),]
#dat_va=exp_va[ids2$probe_id,]
#ids2$mean=apply(dat_va,1,mean)
#ids2=ids2[order(ids2$probe_id,ids2$mean,decreasing = T),]
#ids2=ids2[!duplicated(ids2$symbol),]
#dat_va=dat_va[ids2$probe_id,]
#rownames(dat_va)=ids2$symbol
#write.table(ids2,file = 'ids2.xls',
#            sep='\t')
## idmap3注释后，hub基因少了两个，试一下annoprobe
library(AnnoProbe)
# GPL21827[Accession] - GEO DataSets Result - NCBI - NIH
# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL21827
gpl_va='GPL21827'
probe2gene=idmap(gpl_va,type = 'pipe')
head(probe2gene)
length(unique(probe2gene$symbol))
table(sort(table(probe2gene$symbol)))
ids_va=probe2gene[probe2gene$symbol!='',]
ids_va=probe2gene[probe2gene$probe_id%in%rownames(exp_va),]
exp_va=exp_va[ids_va$probe_id,]
ids_va$mean=apply(exp_va,1,mean)
ids_va=ids_va[order(ids_va$symbol,ids_va$mean,decreasing = T),]
ids_va=ids_va[!duplicated(ids_va$symbol),]
exp_va=exp_va[ids_va$probe_id,]
rownames(exp_va)=ids_va$symbol

# 02 差异基因鉴定----------
setwd("/data/nas1/luchunlin/project/BJTC-251")
if (! dir.exists("./01_DEGs")){
  dir.create("./01_DEGs")
}
setwd("./01_DEGs")
## GSE27276 limma p<0.05 logfc>1 火山图热图
## 02-1样本分组 ---------
group<-data.frame(sample=pd$geo_accession,
                  group=c(rep('POAG',4),rep('control',19),rep('POAG',13)))
write.table(group,file = 'group.xls',
            sep = '\t',
            quote = F,
            row.names = T)
type<-group[,2]
design <- model.matrix(~ -1+factor(type,levels=c('control','POAG'))) 
colnames(design)<-c('control','POAG')
rownames(design)<-group$sample

## 02-2差异分析-----------
# 计算logFC和FDR
library(limma)
# 对每一个基因进行线性模型构建
fit=lmFit(dat2,design)
# 构建比较矩阵
contrast.matrix=makeContrasts(ControlVSMG=POAG-control,levels = design)
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
logFC_cutoff <- 1
DEG$change = as.factor(
  ifelse(DEG$adj.P.Val < 0.05 & abs(DEG$logFC) > logFC_cutoff,
         ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
sig_diff <- subset(DEG,
                   DEG$adj.P.Val < 0.05 & abs(DEG$logFC) > 1)
dim(DEG)
dim(sig_diff)
summary(sig_diff$change)
##171
# DOWN  NOT   UP 
# 110    0   61 
# 输出结果
write.table(DEG,file = "DEG_all.xls",
            quote = F,
            sep = "\t",
            row.names = T)
write.table(sig_diff,file = "DEG_sig.xls",
            quote = F,
            sep = "\t",
            row.names = T)
# 02-3 绘制火山图-----------
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(Ipaper)
library(scales)
library(ggrepel)
dat_rep<-DEG[c('CEACAM5','S100A8','S100A11','DEFB1','S100A9','S100A11'),]
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
    size = 4,
    box.padding = unit(0.5, "lines"),
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )+
  labs(x = "log (Fold Change)",
       y = "-log10 (adj.P.Val)")
volcano_plot
-
# 02-4 绘制热图------------
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(gplots)
group_rt<-group
group_rt<-as.data.frame(group_rt)
rt<-dat2
rownames(group_rt)<-colnames(rt)
## 排序，让rt和group顺序保持一致
colnames(rt)<-group_rt$group
group_rt<-group_rt[order(group_rt$group),]
class(group_rt)
rt<-rt[,order(colnames(rt))]
colnames(rt)<-group_rt$sample
group_rt<-group_rt[,2]
group_rt<-as.data.frame(group_rt)
colnames(group_rt)<-'group'
rownames(group_rt)<-colnames(rt)
heat<-rt[rownames(rt)%in%
           c(head(rownames(subset(sig_diff,sig_diff$logFC>1)),10),head(rownames(subset(sig_diff,sig_diff$logFC<1)),10)),]
x<-t(scale(t(heat)))
ann_colors<-list(
  Group = c(control="lightblue",POAG="darkorange"))
pheatmap(mat=x,
         annotation_col = group_rt,
         color=bluered(100),
         scale = "row",
         annotation_colors = ann_colors,
         fontsize = 10,
         show_colnames = F,
         cluster_cols = F,
         cluster_rows = T)
# 03 WGACN----------
setwd("/data/nas1/luchunlin/project/BJTC-251")
if (! dir.exists("./02_WGACN")){
  dir.create("./02_WGACN")
}
setwd("./02_WGACN")
## 03-1 WGCNA 鉴定疾病相关差异基因，热图展示模块特征基因与青光眼的相关性。强相关基因模块（|Cor| > 0.5, P < 0.05）
library(WGCNA)
library(reshape2)
library(stringr)
options(stringsAsFactors = F)
enableWGCNAThreads()
exprMat<-dat2
dim(exprMat)
type = "unsigned"
corType = "pearson"
corFnc = ifelse(corType=="pearson", cor, bicor)
# 对二元变量，如样本性状信息计算相关性时，
# 或基因表达严重依赖于疾病状态时，需设置下面参数
maxPOutliers = ifelse(corType=="pearson",1,0.05)
# 关联样品性状的二元变量时，设置
robustY = ifelse(corType=="pearson",T,F)
dataExpr <- exprMat[rownames(dat2),]
# 03-1 数据筛选
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
# [1]   36 12263
# 03-2 软阈值筛选
## 样本聚类，查看是否有利群样本
sampleTree = hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
#abline(h = 15, col = "red")
# Determine cluster under the line
# 剪枝算法，cutHeight修剪树枝的高度，minSize集群最小数
#clust = cutreeStatic(sampleTree, cutHeight = 15, minSize = 10)
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
## 7
##一步法网络构建：One-step network construction and module detection##
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
# **0 (grey)**表示**未**分入任何模块的基因。
table(net$colors)
# 0    1    2     3    4    5    6    7    8    9    10   11   12   13   14   15   16   17   18 
#2706 2834 1308  961  768  636  537  338  314  286  275  267  238  172  155  146  126  120   76 
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
group_traits<-group
sample<-group_traits[,1]
sample<-as.data.frame(sample)
group_traits<-group_traits[,-1]
group_traits<-as.data.frame(group_traits)
rownames(group_traits)<-sample$sample
colnames(group_traits)<-"Group"
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
#names (colors) of the modules
#modNames = substring(names(MEs), 3)
#geneModuleMembership = as.data.frame(cor(dataExpr, MEs, use = "p"));
## 算出每个模块跟基因的皮尔森相关系数矩阵
## MEs是每个模块在每个样本里面的值
## datExpr是每个基因在每个样本的表达量
#MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
#names(geneModuleMembership) = paste("MM", modNames, sep="");
#names(MMPvalue) = paste("p.MM", modNames, sep="");
### 计算性状与基因的相关性矩阵
## 只有连续型性状才能只有计算
## 这里把是否属于  表型这个变量用0,1进行数值化。
#POAG = as.data.frame(design_traits[,2]);
#names(POAG) = "POAG"
#geneTraitSignificance = as.data.frame(cor(dataExpr, POAG, use = "p"));
#GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
#names(geneTraitSignificance) = paste("GS.", names(POAG), sep="");
#names(GSPvalue) = paste("p.GS.", names(POAG), sep="")
### 提取模块的基因
## 提取强相关模块基因（|Cor| > 0.5, P < 0.05） cyan  red  greenyellow  midnightblue   brown   lightcyan  magenta  
module1="red"
probes1=colnames(dataExpr)
inModule1=(moduleColors==module1)
modProbes1=probes1[inModule1]
modProbes1<-as.data.frame(modProbes1)
colnames(modProbes1)<-'modProbes'
module2="brown"
probes2=colnames(dataExpr)
inModule2=(moduleColors==module2)
modProbes2=probes1[inModule2]
modProbes2<-as.data.frame(modProbes2)
colnames(modProbes2)<-'modProbes'
modProbes<-rbind(modProbes1,modProbes2)
## 1498个强相关模块基因
## 03-2 GSE27276和WGCNA模块基因取交集获得青光眼差异表达基因#
sig_diff_mod<-sig_diff[rownames(sig_diff)%in%modProbes$modProbes,]
## 119
## 绘制韦恩图
library(VennDiagram)
grid.newpage();
venn_plot<-draw.pairwise.venn(
  area1 = 171,
  area2 = 1498,
  cross.area = 119,
  category = c('DEGs','MOD gene'),
  fill = c("#FFC0CB","#ADD8E6"),
  ext.text = F,
  cex = 3,
  cat.cex = 2,
  cat.pos = c(0,-15),
  cat.dist = 0.03,
  cat.just = list(c(0,1),c(1,1)),
  ext.pos = 30,
  ext.dist = -0.05,
  ext.length = 0.85,
  ext.line.lwd = 2,
  ext.line.lty = "dashed"
)

# 04 青光眼中免疫相关差异表达基因（DEIRG）的鉴定----------
setwd("/data/nas1/luchunlin/project/BJTC-251")
if (! dir.exists("./03_DEIRG")){
  dir.create("./03_DEIRG")
}
setwd("./03_DEIRG")
## Immport、TISIDB和InnateDB取交集（或并集，根据基因数量调整）获得免疫相关基因（Immune-related genes），将Immune-related genes与上述得到的青光眼差异基因做交集，得到青光眼中的免疫相关差异基因（DEIRG）。
Immport<-read_xlsx("/data/nas1/luchunlin/project/BJTC-251/03_DEIRG/Immport_IRGs.xlsx")
TISIDB<-read_xlsx("/data/nas1/luchunlin/project/BJTC-251/03_DEIRG/TISIDB.xlsx")
InnateDB<-read_xls("/data/nas1/luchunlin/project/BJTC-251/03_DEIRG/innatedb.xls")
## 取并集
Immport<-Immport[,1]
TISIDB<-TISIDB[,3]
InnateDB<-InnateDB[,2]
colnames(Immport)<-'symbol'
colnames(TISIDB)<-'symbol'
colnames(InnateDB)<-'symbol'
IRGs<-rbind(Immport,TISIDB,InnateDB)
## 去重
IRGs<-IRGs[!duplicated(IRGs$symbol),]
## 2991
## 与119个疾病相关差异基因取交集
DEIRGs<-sig_diff_mod[rownames(sig_diff_mod)%in%IRGs$symbol,]
## 29 个
## 绘制韦恩图
library(VennDiagram)
grid.newpage();
venn_plot<-draw.pairwise.venn(
  area1 = 119,
  area2 = 2991,
  cross.area = 29,
  category = c('POAG DEGs','IRGs'),
  fill = c("#FFC0CB","#ADD8E6"),
  cex = 3,
  cat.cex = 2,
  cat.pos = c(20,-5),
  cat.dist = 0.05,
  cat.just = list(c(0,1),c(1,1)),
  ext.pos = 30,
  ext.dist = -0.03,
  ext.length = 0.85,
  ext.line.lwd = 1,
  ext.line.lty = "dashed"
)
write.table(DEIRGs,file = "DEIRGs.xls",
            quote = F,
            sep = "\t",
            row.names = T)
# 05 GO/KEGG富集（条形图，气泡图）----------
setwd("/data/nas1/luchunlin/project/BJTC-251")
if (! dir.exists("./04_GO_KEGG")){
  dir.create("./04_GO_KEGG")
}
setwd("./04_GO_KEGG")
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(GOplot)
diff_gene_names <- rownames(DEIRGs)
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
go_result<-as.data.frame(ego)
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

# 06 PPI网络分析筛选hub---------
setwd("/data/nas1/luchunlin/project/BJTC-251")
if (! dir.exists("./05_PPI")){
  dir.create("./05_PPI")
}
setwd("./05_PPI")

# 使用Cytohubba degree 10  degree>2的7个. 分别为：SLPI,S100A8,CEACAM5,S100A9,S100A11,DEFB1,LCN2

# 07 诊断基因的筛选--------------
setwd("/data/nas1/luchunlin/project/BJTC-251")
if (! dir.exists("./06_hub")){
  dir.create("./06_hub")
}
setwd("./06_hub")
## 07-1 ROC曲线验证--------
## 将7个基因的矩阵提取出来
hubgene<-c('SLPI','S100A8','CEACAM5','S100A9','S100A11','DEFB1','LCN2')
hubgene<-as.data.frame(hubgene)
hub_exp<-dat2[rownames(dat2)%in%hubgene$hubgene,]
library(pROC)
library(ggplot2)
hub_exp2<-t(hub_exp)
hub_exp2<-cbind(design,hub_exp2)
hub_exp2<-cbind(group_traits,hub_exp2)
## 绘制ROC曲线
## SLPI
roc_SLPI<-roc(hub_exp2$Group,hub_exp2$SLPI,
              levels=c("control","POAG"))
print(ci.thresholds(roc_SLPI))
## GDF15
#roc_GDF15<-roc(hub_exp2$Group,hub_exp2$GDF15,
#              levels=c("control","POAG"))
## S100A8
roc_S100A8<-roc(hub_exp2$Group,hub_exp2$S100A8,
                levels=c("control","POAG"))
print(ci.thresholds(roc_S100A8))

## CEACAM5
roc_CEACAM5<-roc(hub_exp2$Group,hub_exp2$CEACAM5,
                levels=c("control","POAG"))
print(ci.thresholds(roc_CEACAM5))
## GRP
#roc_GRP<-roc(hub_exp2$Group,hub_exp2$GRP,
#                 levels=c("control","POAG"))
## S100A9
roc_S100A9<-roc(hub_exp2$Group,hub_exp2$S100A9,
             levels=c("control","POAG"))
print(ci.thresholds(roc_S100A9))

## S100A11
roc_S100A11<-roc(hub_exp2$Group,hub_exp2$S100A11,
             levels=c("control","POAG"))
print(ci.thresholds(roc_S100A11))
## S100A14
#roc_S100A14<-roc(hub_exp2$Group,hub_exp2$S100A14,
#             levels=c("control","POAG"))
## DEFB1
roc_DEFB1<-roc(hub_exp2$Group,hub_exp2$DEFB1,
             levels=c("control","POAG"))
print(ci.thresholds(roc_DEFB1))
## LCN2
roc_LCN2<-roc(hub_exp2$Group,hub_exp2$LCN2,
             levels=c("control","POAG"))
print(ci.thresholds(roc_LCN2))
## SLPI
plot(roc_SLPI,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="SLPI ROC curve",
     col="#FF2E63",
     legacy.axes=T)
## GDF15
plot(roc_GDF15,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="GDF15 ROC curve",
     col="#FF2E63",
     legacy.axes=T)
## S100A8
plot(roc_S100A8,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="S100A8 ROC curve",
     col="#FF2E63",
     legacy.axes=T)
## CEACAM5
plot(roc_CEACAM5,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="CEACAM5 ROC curve",
     col="#FF2E63",
     legacy.axes=T)
## GRP
plot(roc_GRP,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="GRP ROC curve",
     col="#FF2E63",
     legacy.axes=T)
## S100A9
plot(roc_S100A9,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="S100A9 ROC curve",
     col="#FF2E63",
     legacy.axes=T)
## S100A11
plot(roc_S100A11,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="S100A11 ROC curve",
     col="#FF2E63",
     legacy.axes=T)
## S100A14
plot(roc_S100A14,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="S100A14 ROC curve",
     col="#FF2E63",
     legacy.axes=T)
## DEFB1
plot(roc_DEFB1,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="DEFB1 ROC curve",
     col="#FF2E63",
     legacy.axes=T)
## LCN2
plot(roc_LCN2,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="LCN2 ROC curve",
     col="#FF2E63",
     legacy.axes=T)
## 7 个基因全部AUC>0.7 SLPI,S100A8,CEACAM5,S100A9,S100A14,DEFB1,LCN2

## 07-2 诊断基因的差异表达-----------

## 数据清洗
hub_exp3<-hub_exp
hub_exp3$Symbol<-rownames(hub_exp3)
hub_exp3<-hub_exp3[,c(37,1:36)]
hub_exp3<-gather(hub_exp3,
                 key = sample,
                 value = expr,
                 -c("Symbol"))
library(rstatix)
stat.test<-hub_exp3%>%
  group_by(Symbol)%>%
  t_test(expr ~ Group)%>%
  adjust_pvalue(method = 'fdr')
write.table(stat.test,file = 'hubexp_result.xls',
            sep = '\t',
            row.names = F)
## 样本分组
group_exp<-data.frame(group=c(rep('POAG',28),rep('control',133),rep('POAG',91)))
colnames(group_exp)<-"Group"
hub_exp3<-cbind(group_exp,hub_exp3)

##分面图形
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
exp_plot <- ggplot(hub_exp3,aes(x = Symbol, y = expr, fill = Group)) +
  #geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#4682B4","#CD3700"), name = "Group")+
  labs(title="", x="", y = "",size=20) +
  stat_compare_means(data = hub_exp3,
                     mapping = aes(group = Group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F) +
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
        panel.grid.minor = element_blank())+facet_wrap(~Symbol,scales = "free",nrow = 2) 
exp_plot


# 分面图形
exp_boxplot<-ggplot(hub_exp3,aes(x = Symbol, y = expr, fill = Group)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "Expression",expand = c(0.1,0.1))+
  scale_x_discrete(name = "group") +
  ggtitle("Expression of hub genes") +
  stat_compare_means(data = hub_exp3,
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
exp_boxplot2<-exp_boxplot+scale_fill_lancet()+facet_wrap(~Symbol)
exp_boxplot2

# 08 外部数据集验证------------
setwd("/data/nas1/luchunlin/project/BJTC-251")
if (! dir.exists("./07_hub_va")){
  dir.create("./07_hub_va")
}
setwd("./07_hub_va")
## 验证集换成GSE138125，分析时先用这个验证
## 08-1 表达情况--------
## 将表达矩阵提取出来
hub_exp_va<-exp_va[rownames(exp_va)%in%hubgene$hubgene,]
hub_exp_va3<-hub_exp_va
hub_exp_va3<-as.data.frame(hub_exp_va3)
hub_exp_va3$Symbol<-rownames(hub_exp_va3)
hub_exp_va3<-gather(hub_exp_va3,
                    key = sample,
                    value = expr,
                    -c("Symbol"))
stat.test.va<-hub_exp_va3%>%
  group_by(Symbol)%>%
  t_test(expr ~ Group)%>%
  adjust_pvalue(method = 'fdr')
write.table(stat.test.va,file = 'hubexp_result_va.xls',
            sep = '\t',
            row.names = F)
## 样本分组
group_exp_va<-data.frame(group=c(rep('control',28),rep('POAG',28)))
colnames(group_exp_va)<-"Group"
hub_exp_va3<-cbind(group_exp_va,hub_exp_va3)
# 分面图形
exp_boxplot_va<-ggplot(hub_exp_va3,aes(x = Symbol, y = expr, fill = Group)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "Expression",expand = c(0.1,0.1))+
  scale_x_discrete(name = "group") +
  ggtitle("Expression of hub genes") +
  stat_compare_means(data = hub_exp_va3,
                     mapping = aes(group = Group),
                     label ="p.signif",
                     method = 't.test',
                     paired = F) +
  theme_bw() +
  #  geom_signif(comparisons = list(c("control","POAG")),
  #              test = t.test,
  #              map_signif_level = T)+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 18),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13))
exp_boxplot_va
exp_boxplot_va2<-exp_boxplot_va+scale_fill_lancet()+facet_wrap(~Symbol,scales = "free")
exp_boxplot_va2


## 08-2 ROC曲线验证----------
group_va<-data.frame(Group=c(rep('control',4),rep('POAG',4)))
type_va<-group_va$Group
design_va <- model.matrix(~ -1+factor(type_va,levels=c('control','POAG'))) 
colnames(design_va)<-c('control','POAG')
rownames(design_va)<-rownames(pd_va)
hub_exp_va2<-t(hub_exp_va)
hub_exp_va2<-cbind(design_va,hub_exp_va2)
hub_exp_va2<-cbind(group_va,hub_exp_va2)
## SLPI
roc_SLPI_va<-roc(hub_exp_va2$Group,hub_exp_va2$SLPI,
               levels=c("control","POAG"))
roc_SLPI$predictor
print(ci.thresholds(roc_S100A9_va))

## GDF15
roc_GDF15_va<-roc(hub_exp_va2$Group,hub_exp_va2$GDF15,
                 levels=c("control","POAG"))
## S100A8
roc_S100A8_va<-roc(hub_exp_va2$Group,hub_exp_va2$S100A8,
                 levels=c("control","POAG"))
## CEACAM5
roc_CEACAM5_va<-roc(hub_exp_va2$Group,hub_exp_va2$CEACAM5,
                 levels=c("control","POAG"))
## GRP
roc_GRP_va<-roc(hub_exp_va2$Group,hub_exp_va2$GRP,
                 levels=c("control","POAG"))
## S100A9
roc_S100A9_va<-roc(hub_exp_va2$Group,hub_exp_va2$S100A9,
                 levels=c("control","POAG"))
## S100A11
roc_S100A11_va<-roc(hub_exp_va2$Group,hub_exp_va2$S100A11,
                 levels=c("control","POAG"))
## S100A14
roc_S100A14_va<-roc(hub_exp_va2$Group,hub_exp_va2$S100A14,
                 levels=c("control","POAG"))
## DEFB1
roc_DEFB1_va<-roc(hub_exp_va2$Group,hub_exp_va2$DEFB1,
                 levels=c("control","POAG"))
## LCN2
roc_LCN2_va<-roc(hub_exp_va2$Group,hub_exp_va2$LCN2,
                 levels=c("control","POAG"),direction='>')

## SLPI
plot(roc_SLPI_va,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="SLPI ROC curve",
     col="#FF2E63",
     legacy.axes=T)
## GDF15
plot(roc_GDF15_va,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="GDF15 ROC curve",
     col="#FF2E63",
     legacy.axes=T)
## S100A8
plot(roc_S100A8_va,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="S100A8 ROC curve",
     col="#FF2E63",
     legacy.axes=T)
## CEACAM5
plot(roc_CEACAM5_va,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="CEACAM5 ROC curve",
     col="#FF2E63",
     legacy.axes=T)
## GRP
plot(roc_GRP_va,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="GRP ROC curve",
     col="#FF2E63",
     legacy.axes=T)
## S100A9
plot(roc_S100A9_va,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="S100A9 ROC curve",
     col="#FF2E63",
     legacy.axes=T)
## S100A11
plot(roc_S100A11_va,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="S100A11 ROC curve",
     col="#FF2E63",
     legacy.axes=T)
## S100A14
plot(roc_S100A14_va,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="S100A14 ROC curve",
     col="#FF2E63",
     legacy.axes=T)
## DEFB1
plot(roc_DEFB1_va,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="DEFB1 ROC curve",
     col="#FF2E63",
     legacy.axes=T)
## LCN2
plot(roc_LCN2_va,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,
     #auc.polygon=T,
     #auc.polygon.con="#fff7f7",
     grid=c(0.5,0.2),
     grid.col=c("black","black"),
     #print.thres=T,
     main="LCN2 ROC curve",
     col="#FF2E63",
     legacy.axes=T)
## 6个基因通过了验证SLPI,S100A8,CEACAM5,S100A9,S100A11,DEFB1
# 09 诊断基因表达相关性分析-----------
setwd("/data/nas1/luchunlin/project/BJTC-251")
if (! dir.exists("./08_correlation")){
  dir.create("./08_correlation")
}
setwd("./08_correlation")
library(ggcorrplot)
library(corrplot)
hub_gene_final<-hubgene[-7,]
hub_gene_final<-as.data.frame(hub_gene_final)
## 提取诊断基因表达矩阵
exp_final<-dat2[rownames(dat2)%in%hub_gene_final$hub_gene_final,]
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
write.table(hub_corr,file = 'hub_corr.xls',
            sep = '\t',
            row.names = T)

#10 诊断基因功能相似性分析（GOSemSim包）----------
setwd("/data/nas1/luchunlin/project/BJTC-251")
if (! dir.exists("./09_relationship")){
  dir.create("./09_relationship")
}
setwd("./09_relationship")
library(GOSemSim)
library(org.Hs.eg.db)
bp<-godata('org.Hs.eg.db',ont = "BP",computeIC = F)
cc<-godata('org.Hs.eg.db',ont = "CC",computeIC = F)
mf<-godata('org.Hs.eg.db',ont = "MF",computeIC = F)
geneid2<-hub_gene_final$hub_gene_final
gene_transform2<-bitr(geneid2,fromType = "SYMBOL",
                     toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db)
genelist<-gene_transform2
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
write.table(fsim,file = 'relationship.xls',
            sep = '\t',
            row.names = T)
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
# 10 诊断基因调控网络预测
setwd("/data/nas1/luchunlin/project/BJTC-251")
if (! dir.exists("./10_network")){
  dir.create("./10_network")
}
setwd("./10_network")
## 
## 使用mirTarBase  JASPAR数据库对青光眼诊断基因的miRNA和TF进行预测，将诊断基因及其miRNA和TF的数据整合到一个调控网络中（根据实际情况，数据库可进行替换）。
## SLPI,GDF15,S100A8,CEACAM5,GRP,S100A9,S100A11,S100A14,DEFB1











