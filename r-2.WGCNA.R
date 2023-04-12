rm(list = ls())
# 03 WGCNA----------
setwd("/data/nas1/luchunlin/project/BJTC-321")
if (! dir.exists("./02_WGCNA")){
  dir.create("./02_WGCNA")
}
setwd("./02_WGCNA")

## 03-1 GSVA-------
# library(GSVA)
# inflammatory<-read_xlsx('inflammatory.xlsx')
# inflammatory$inflammator<-c(rep('inflammatory',200))
library(tidyverse)
library(lance)
dat<-read.delim2('/data/nas1/luchunlin/project/BJTC-321/00_rawdata/dat.final.xls',row.names = 1)%>%lc.tableToNum()
group<-read.delim2('/data/nas1/luchunlin/project/BJTC-321/00_rawdata/group.xls')
table(group$group)
# asthma.sample<-group$sample[which(group$group=='asthma')]
# dat.asthma<-dat[,asthma.sample]%>%as.matrix()
# gene_set <- read.table("/data/nas1/luchunlin/pipeline/ssGSEA/mmc3.txt",
#                        header = T,
#                        sep ="\t")
# gene_list <- split(as.matrix(inflammatory)[,1],
#                    inflammatory[,2])
# 
# inflammatory_score = gsva(dat.asthma, gene_list, 
#                     method = "gsva", 
#                     ssgsea.norm = TRUE, 
#                     verbose = TRUE)
# write.table(inflammatory_score,
#             file = "inflammatory_result.xls",
#             sep = "\t",
#             quote = F)
# inflammatory_score<-t(inflammatory_score)%>%as.data.frame()
# group2<-data.frame(sample=rownames(inflammatory_score),group=ifelse(inflammatory_score$inflammatory>median(inflammatory_score$inflammatory),'High_inflammatory_score','Low_inflammatory_score'))
# # 03-2 WGCNA------
library(WGCNA)
library(reshape2)
library(stringr)
options(stringsAsFactors = F)
enableWGCNAThreads()
exprMat<-dat
dim(exprMat)
type = "unsigned"
corType = "pearson"
corFnc = ifelse(corType=="pearson", cor, bicor)
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
# [1]105 13797

## 04-2 软阈值筛选----
## 样本聚类，查看是否有利群样本
sampleTree = hclust(dist(dataExpr,method = 'euclidean'), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")

abline(h =65, col = "red")
# Determine cluster under the line
# 剪枝算法，cutHeight修剪树枝的高度，minSize集群最小数
clust = cutreeStatic(sampleTree, cutHeight = 65, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
## 符合要求的数据
dataExpr = dataExpr[keepSamples,]
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
## 7
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
# **0 (grey)**表示**未**分入任何模块的基因。
table(net$colors)  
#    0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20  21 22
## 3526 1838 1364 1209 1050  651  603  595  535  392  339  301  277  225  147  138  133  121  118   81   65   59   30 
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
par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor[-23,],
               xLabels = names(design_traits),
               yLabels = names(MEs[,-23]),
               ySymbols = names(MEs[,-23]),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix[-23,],
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
sizeGrWindow(8,7)
par(mfrow = c(1,1))
GeneSignificance=abs(moduleTraitCor[,1])
moduleTraitCor
color<-rownames(moduleTraitCor)
color<-gsub('ME','',color)
library(ggthemes)
plotModuleSignificance(GeneSignificance[-23],color[-23])

## 04-7 筛选强相关性模块------------
## blue
module='lightgreen'
probes=colnames(dataExpr)
inModule=(moduleColors==module)
modProbes=probes[inModule]
modGenes<-as.data.frame(modProbes)
colnames(modGenes)<-'modgene'
dim(modGenes)
#1413

## 加1个magenta
module2='magenta'
probes=colnames(dataExpr)
inModule2=(moduleColors==module2)
modProbes2=probes[inModule2]
modGenes2<-as.data.frame(modProbes2)
colnames(modGenes2)<-'modgene'
modGenes<-rbind(modGenes,modGenes2)


module3='blue'
probes=colnames(dataExpr)
inModule3=(moduleColors==module3)
modProbes3=probes[inModule3]
modGenes3<-as.data.frame(modProbes3)
colnames(modGenes3)<-'modgene'
modGenes<-rbind(modGenes,modGenes3)
#1874
write.table(modGenes,file = 'modGene.xls',sep = '\t',row.names = F,quote = F)

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
  geneTraitCor = as.data.frame(cor(dataExpr, design_traits, use = "p"))
  geneTraitP = as.data.frame(corPvalueStudent(
    as.matrix(geneTraitCor), nSamples))
} else {
  geneTraitCorA = bicorAndPvalue(dataExpr, design_traits, robustY=robustY)
  geneTraitCor = as.data.frame(geneTraitCorA$bicor)
  geneTraitP   = as.data.frame(geneTraitCorA$p)
}
### lightgreen------
module = "lightgreen"
pheno = "asthma"
modNames = substring(colnames(MEs_col), 3)
# 获取关注的列
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(design_traits))
# 获取模块内的基因
moduleGenes = moduleColors == module

sizeGrWindow(7, 7)
par(mfrow = c(1,1),
    mai = c(1,1,1,0.5))
# 与性状高度相关的基因，也是与性状相关的模型的关键基因
pdf('lightgreen_gene_cor.pdf',w=6,h=6)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()
png('ligthgreen_gene_cor.png',w=600,h=600)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()
### 最相关模块GO KEGG富集
module='lightgreen'
probes=colnames(dataExpr)
inModule=(moduleColors==module)
modProbes=probes[inModule]
GOGenes<-as.data.frame(modProbes)
colnames(GOGenes)<-'lightgreen'
dim(GOGenes)

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(GOplot)

gene_transform <- bitr(GOGenes$lightgreen,
                       fromType = "SYMBOL",
                       toType = c("ENTREZID",'ENSEMBL'),
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
go_bar<-barplot(ego, showCategory=5, split="ONTOLOGY",label_format = 50) +
  facet_grid(ONTOLOGY ~ ., scales = "free")
go_bar
## GO 圈图
go_result<-read.table('GO.xls',header = T,sep = '\t',check.names = F)
go2=data.frame(Category='ALL',ID=go_result$ID,Term=go_result$Description,Genes=gsub("/", ", ", go_result$geneID), adj_pval = go_result$p.adjust)
df<-diff[GOGenes$lightgreen,]
idfc2<-na.omit(df)
genelist<-data.frame(ID=rownames(idfc2),logFC=idfc2$logFC)
rownames(genelist)=genelist[,1]
circ<-circle_dat(go2,genelist)
circ<-na.omit(circ)
circ$logFC<-as.numeric(circ$logFC)
go_cir<-GOCircle(circ,rad1 = 2.5,rad2 = 3.5,label.size = 4,nsub = 10)
ggsave('GO_cir.png',go_cir,width =10,height = 6)
ggsave('GO_cir.pdf',go_cir,width =10,height = 6)
## KEGG富集分析（气泡图）
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "hsa",
                 pAdjustMethod = "none",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.2)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk@result,file = "KEGG.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk, showCategory=15,label_format = 50)
kk_dot
ggsave('KK_dot.png',kk_dot,width =8,height = 5)
ggsave('KK_dot.pdf',kk_dot,width =8,height = 5)
kk2<-data.frame(Category = "ALL",ID = kk$ID,Term = kk$Description, Genes = gsub("/", ", ", kk$geneID), adj_pval = kk$p.adjust)
sigind<-diff[rownames(diff)%in%GOGenes$lightgreen,]
genelist<-data.frame(ID=rownames(sigind),
                     logFC=sigind$logFC)
genelist$logFC<-as.numeric(genelist$logFC)
rownames(genelist)<-genelist$ID
genelist<-genelist[order(genelist$logFC,decreasing = T),]

## 读取logFC文件
circ<-circle_dat(kk2,genelist)
# GOBubble(circ, labels = 3,table.legend =F)
# GOCircle(circ,rad1=2.5,rad2=3.5,label.size=4,nsub=10) 
#termNum<-20
#geneNum<-nrow(genelist)
chord<-chord_dat(circ,genelist)
trace(GOChord,edit = T)
kegg_chord<-GOChord(chord,
                    gene.order = 'logFC',
                    gene.space = 0.3,
                    gene.size = 4,
                    space = 0.01,
                    lfc.col = c('red','white','blue'),
                    process.label = 10)
kegg_chord

### blue---------
module = "blue"
pheno = "asthma"
modNames = substring(colnames(MEs_col), 3)
# 获取关注的列
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(design_traits))
# 获取模块内的基因
moduleGenes = moduleColors == module

sizeGrWindow(7, 7)
par(mfrow = c(1,1),
    mai = c(1,1,1,0.5))
pdf('blue_gene_cor.pdf',w=6,h=6)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()
png('blue_gene_cor.png',w=600,h=600)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()
### 最相关模块GO KEGG富集
module='blue'
probes=colnames(dataExpr)
inModule=(moduleColors==module)
modProbes=probes[inModule]
GOGenes<-as.data.frame(modProbes)
colnames(GOGenes)<-'blue'
dim(GOGenes)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(GOplot)

gene_transform <- bitr(GOGenes$blue,
                       fromType = "SYMBOL",
                       toType = c("ENTREZID",'ENSEMBL'),
                       OrgDb = "org.Hs.eg.db")
ego <- enrichGO(gene = gene_transform$ENTREZID,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.3,
                qvalueCutoff = 1,
                readable = TRUE)
write.table(ego,file = "GO(blue).xls",sep = "\t",quote = F,row.names = F)
# 展示富集最显著的 GO term
go_bar<-barplot(ego, showCategory=5, split="ONTOLOGY",label_format = 50) +
  facet_grid(ONTOLOGY ~ ., scales = "free")
go_bar
## GO 圈图
go_result<-read.table('GO(blue).xls',header = T,sep = '\t',check.names = F)
go2=data.frame(Category='ALL',ID=go_result$ID,Term=go_result$Description,Genes=gsub("/", ", ", go_result$geneID), adj_pval = go_result$p.adjust)
df<-diff[GOGenes$blue,]
idfc2<-na.omit(df)
genelist<-data.frame(ID=rownames(idfc2),logFC=idfc2$logFC)
rownames(genelist)=genelist[,1]
circ<-circle_dat(go2,genelist)
circ<-na.omit(circ)
circ$logFC<-as.numeric(circ$logFC)
go_cir<-GOCircle(circ,rad1 = 2.5,rad2 = 3.5,label.size = 4,nsub = 10)
ggsave('GO_cir(blue).png',go_cir,width =20,height = 12)
ggsave('GO_cir(blue).pdf',go_cir,width =20,height = 12)
## KEGG富集分析（气泡图）
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "hsa",
                 pAdjustMethod = "none",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.2)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk@result,file = "KEGG(blue).xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk, showCategory=15,label_format = 50)
kk_dot
# ggsave('KK_dot.png',kk_dot,width =8,height = 5)
# ggsave('KK_dot.pdf',kk_dot,width =8,height = 5)
kk2<-data.frame(Category = "ALL",ID = kk$ID,Term = kk$Description, Genes = gsub("/", ", ", kk$geneID), adj_pval = kk$p.adjust)
sigind<-diff[rownames(diff)%in%GOGenes$blue,]
genelist<-data.frame(ID=rownames(sigind),
                     logFC=sigind$logFC)
genelist$logFC<-as.numeric(genelist$logFC)
rownames(genelist)<-genelist$ID
genelist<-genelist[order(genelist$logFC,decreasing = T),]

## 读取logFC文件
circ<-circle_dat(kk2,genelist)
circ<-circ[c(1:255),]
# GOBubble(circ, labels = 3,table.legend =F)
# GOCircle(circ,rad1=2.5,rad2=3.5,label.size=4,nsub=10) 
#termNum<-20
#geneNum<-nrow(genelist)
chord<-chord_dat(circ,genelist)
trace(GOChord,edit = T)
kegg_chord<-GOChord(chord,
                    gene.order = 'logFC',
                    gene.space = 0.3,
                    gene.size = 4,
                    space = 0.01,
                    lfc.col = c('red','white','blue'),
                    process.label = 10)
kegg_chord


### magenta---------
module = "magenta"
pheno = "asthma"
modNames = substring(colnames(MEs_col), 3)
# 获取关注的列
module_column = match(module, modNames)
pheno_column = match(pheno,colnames(design_traits))
# 获取模块内的基因
moduleGenes = moduleColors == module

sizeGrWindow(7, 7)
par(mfrow = c(1,1),
    mai = c(1,1,1,0.5))
pdf('magenta_gene_cor.pdf',w=6,h=6)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()
png('magenta_gene_cor.png',w=600,h=600)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, module_column]),
                   abs(geneTraitCor[moduleGenes, pheno_column]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = paste("Gene significance for", pheno),
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
dev.off()
### 最相关模块GO KEGG富集
module='magenta'
probes=colnames(dataExpr)
inModule=(moduleColors==module)
modProbes=probes[inModule]
GOGenes<-as.data.frame(modProbes)
colnames(GOGenes)<-'magenta'
dim(GOGenes)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(GOplot)

gene_transform <- bitr(GOGenes$magenta,
                       fromType = "SYMBOL",
                       toType = c("ENTREZID",'ENSEMBL'),
                       OrgDb = "org.Hs.eg.db")
ego <- enrichGO(gene = gene_transform$ENTREZID,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 1,
                readable = TRUE)
write.table(ego,file = "GO(magenta).xls",sep = "\t",quote = F,row.names = F)
# 展示富集最显著的 GO term
go_bar<-barplot(ego, showCategory=5, split="ONTOLOGY",label_format = 50) +
  facet_grid(ONTOLOGY ~ ., scales = "free")
go_bar
## GO 圈图
go_result<-read.table('GO(magenta).xls',header = T,sep = '\t',check.names = F)
go2=data.frame(Category='ALL',ID=go_result$ID,Term=go_result$Description,Genes=gsub("/", ", ", go_result$geneID), adj_pval = go_result$p.adjust)
df<-diff[GOGenes$magenta,]
idfc2<-na.omit(df)
genelist<-data.frame(ID=rownames(idfc2),logFC=idfc2$logFC)
rownames(genelist)=genelist[,1]
circ<-circle_dat(go2,genelist)
circ<-na.omit(circ)
circ$logFC<-as.numeric(circ$logFC)
go_cir<-GOCircle(circ,rad1 = 2.5,rad2 = 3.5,label.size = 4,nsub = 10)
ggsave('GO_cir(magenta).png',go_cir,width =15,height = 8)
ggsave('GO_cir(magenta).pdf',go_cir,width =15,height = 8)
## KEGG富集分析（气泡图）
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "hsa",
                 pAdjustMethod = "none",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.2)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk@result,file = "KEGG(magenta).xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk, showCategory=15,label_format = 50)
kk_dot
# ggsave('KK_dot.png',kk_dot,width =8,height = 5)
# ggsave('KK_dot.pdf',kk_dot,width =8,height = 5)
kk2<-data.frame(Category = "ALL",ID = kk$ID,Term = kk$Description, Genes = gsub("/", ", ", kk$geneID), adj_pval = kk$p.adjust)
sigind<-diff[rownames(diff)%in%GOGenes$magenta,]
genelist<-data.frame(ID=rownames(sigind),
                     logFC=sigind$logFC)
genelist$logFC<-as.numeric(genelist$logFC)
rownames(genelist)<-genelist$ID
genelist<-genelist[order(genelist$logFC,decreasing = T),]

## 读取logFC文件
circ<-circle_dat(kk2,genelist)
# GOBubble(circ, labels = 3,table.legend =F)
# GOCircle(circ,rad1=2.5,rad2=3.5,label.size=4,nsub=10) 
#termNum<-20
#geneNum<-nrow(genelist)
chord<-chord_dat(circ,genelist)
trace(GOChord,edit = T)
kegg_chord<-GOChord(chord,
                    gene.order = 'logFC',
                    gene.space = 0.3,
                    gene.size = 4,
                    space = 0.01,
                    lfc.col = c('red','white','magenta'),
                    process.label = 10)
kegg_chord

