rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-406-12/")
if (! dir.exists("./16_scRNA")){
  dir.create("./16_scRNA")
}
setwd("./16_scRNA")
library(Seurat)
library(GEOquery)
library(data.table)
library(tidyverse)
library(lance)
## PART A导入数据--------
## 导入数据，创建seurat对象
dat <- fread('GSE131907_Lung_Cancer_raw_UMI_matrix.txt.gz',header = T)%>%
  column_to_rownames(var = 'Index')
colnames(dat)
##只要11个lung T样本 
smps <- data.frame(cell=colnames(dat),sample=colnames(dat) %>% str_sub(18,26))
sample.list = c('LUNG_T06','LUNG_T08','LUNG_T09','LUNG_T18','LUNG_T19','LUNG_T20','LUNG_T25','LUNG_T28','LUNG_T30','LUNG_T31','LUNG_T34')
sample.final <- smps[smps$sample%in%sample.list,]
dat <- dat[,colnames(dat)%in%sample.final$cell]

scRNA <- CreateSeuratObject(counts = dat,min.cells = 3,min.features = 200)
##初步过滤，<=3个细胞中表达的基因，>=200个基因的细胞
scRNA

## ##PART B质控-------
setwd(('/data/nas1/luchunlin/project/BJTC-406-12/16_scRNA/'))
dir.create("01_QC")
setwd('../16_scRNA/01_QC/')
### 质控前 ###

scRNA[['percent.mt']]<-PercentageFeatureSet(scRNA,pattern = "^MT-")
head(scRNA@meta.data)
summary(scRNA@meta.data)
# 计算基因含量，MT为线粒体
theme.set = theme(
  #axis.title.x=element_blank(),
  axis.title = element_text(size = 20, face = "bold", family = "Times"),
  axis.text.x = element_text(size = 12,  family = "Times"),
  axis.text.y = element_text(size = 14,  family = "Times"),
  legend.text = element_text(size = 14, family = "Times"),
  legend.title = element_text(size = 16,face='bold',family = "Times"),
  text = element_text(family = "Times"))
plots <- VlnPlot(scRNA,
                 features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3) + NoLegend() #+theme.set 
plots
# violin <- wrap_plots(plots = plots, nrow=1,ncol = 3)    
# violin

## 可视化。线粒体基因含量占比5%以下的细胞保留。
## 可视化RNA-基因含量（RNA-feature）
plot1<-FeatureScatter(scRNA,feature1 = 'nCount_RNA',feature2 = 'percent.mt')
plot2<-FeatureScatter(scRNA,feature1 = 'nCount_RNA',feature2 = 'nFeature_RNA')
plot1+plot2

#去除线粒体基因表达比例过高的细胞，和一些极值细胞
scRNA <- subset(scRNA,
                subset = nFeature_RNA > 200 & nFeature_RNA < 2500&
                  percent.mt < 5)
length(colnames(scRNA))#29747细胞数量
length(rownames(scRNA))# 24321基因数量
VlnPlot(scRNA,
        features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3) + NoLegend()

pdf(file = '01.QC.pdf',w=8,h=5)
VlnPlot(scRNA,
        features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3) + NoLegend()
dev.off()
png(file = '01.QC.png',w=700,h=400)
VlnPlot(scRNA,
        features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3) + NoLegend()

dev.off()

## PART C数据降维--------
setwd(('/data/nas1/luchunlin/project/BJTC-406-12/16_scRNA/'))
if (! dir.exists("./02_PCA")){dir.create("./02_PCA")}
setwd("./02_PCA")
#标准化
scRNA.norm<-NormalizeData(scRNA,normalization.method = "LogNormalize",scale.factor = 10000)
##scRNA.norm<-scRNA
## PART D PCA UMAP降维聚类-------
## 缩放数据
## 筛选2000个高变基因
all.genes<-rownames(scRNA.norm)
scRNA.norm<-FindVariableFeatures(scRNA.norm,selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(scRNA.norm), 10)
# plot variable features with and without labels
plot3 <- VariableFeaturePlot(scRNA.norm)
plot4 <- LabelPoints(plot = plot3, 
                     points = top10, 
                     repel = T)
plot4
library(Ipaper)
write_fig(plot4,
          file = "01.feature_selection.pdf",
          width = 5,
          height = 6,
          devices = NULL,
          res = 600,
          show = F)
write_fig(plot4,
          file = "01.feature_selection.png",
          width = 5,
          height = 6,
          devices = NULL,
          res = 300,
          show = F)

scRNA.nor.sca<-ScaleData(scRNA.norm)
## PCA降维
scRNA.norm.pca<-RunPCA(scRNA.nor.sca,features = VariableFeatures(object = scRNA.nor.sca))
print(scRNA.norm.pca[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(scRNA.norm.pca, dims = 1:2, reduction = "pca")
DimPlot(scRNA.norm.pca, reduction = "pca")
pdf("02.PCA.pdf",w=7,h=5)
DimPlot(scRNA.norm.pca, reduction = "pca")
dev.off()
png("02.PCA.png",w=600,h=400)
DimPlot(scRNA.norm.pca, reduction = "pca")
dev.off()
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
library(Ipaper)
write_fig(plot5,
          file = "03.pca_cluster.pdf",
          width = 13,
          height = 8,
          devices = NULL,
          res = 600,
          show = F)
write_fig(plot5,
          file = "03.pca_cluster.png",
          width = 13,
          height = 8,
          devices = NULL,
          res = 300,
          show = F)
plot6 <- ElbowPlot(scRNA.norm.pca, ndims = 50)
plot6
write_fig(plot6,
          file = "04.pca_sd.pdf",
          width = 8,
          height = 6,
          devices = NULL,
          res = 600,
          show = F)
write_fig(plot6,
          file = "04.pca_sd.png",
          width = 8,
          height = 6,
          devices = NULL,
          res = 300,
          show = F)
scRNA.norm.pca.c<-FindNeighbors(scRNA.norm.pca,dims = 1:50)
scRNA.norm.pca.c<-FindClusters(scRNA.norm.pca.c,resolution = 0.4)
UMAP<-RunUMAP(scRNA.norm.pca.c,dims = 1:50)
DimPlot(UMAP,reduction = 'umap')
pdf('05.UMAP.pdf',w=8,h=7)
DimPlot(UMAP,reduction = 'umap')
dev.off()
png('05.UMAP.png',w=700,h=600)
DimPlot(UMAP,reduction = 'umap')
dev.off()
#-tsne-----
# set.seed(123)
# TSNE <- RunTSNE(scRNA.norm.pca.c, dims = 1:27)
# head(TSNE@reductions$tsne@cell.embeddings)
# p5 <- DimPlot(TSNE, reduction = "tsne", label = T)
# p5
# saveRDS(TSNE,file = 'TSNE.rds')
## 寻找markergene------
# UMAP$group
# DefaultAssay(UMAP) = "RNA"
# marker.genes = lapply(0:21, function(x){
#   FindConservedMarkers(UMAP, ident.1 = x, grouping.var = "group", 
#                        logfc.threshold = 1, min.pct = 0.1, only.pos = T)
# })
# saveRDS(marker.genes, "marker.genes.rds")
# names(marker.genes) = paste0("cluster", seq(0,21))
# res.marker = lapply(marker.genes, function(x){
#   cbind(gene = rownames(x), x)
# })
# #install.packages('WriteXLS')
# library(WriteXLS)
# WriteXLS(res.marker, "MarkerGenes.xls")
# cluster1.markers<-FindMarkers(UMAP,ident.1 = 1,min.pct = 0.25)
# head(cluster1.markers,n=5)
all.markers<-FindAllMarkers(UMAP,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.25,test.use = 'wilcox',return.thresh = 0.05)
all.markers%>%group_by(cluster)%>%top_n(n=2,wt=avg_log2FC)

write.table(all.markers,file = 'DEGs.xls',sep = '\t',row.names = F,quote = F)
top10<-all.markers%>%group_by(cluster)%>%top_n(n=5,wt=avg_log2FC)
pdf('06.topheat.pdf',w=5,h=4)
DoHeatmap(UMAP,features = top10$gene,label = F)+NoLegend()+
  labs(title="", x="Cells separated by clusters", y = "5 most-sig.DEGs in each cluster",size=40)+
  theme(axis.text.y = element_blank())
dev.off()
png('06.topheat.png',w=600,h=500)
DoHeatmap(UMAP,features = top10$gene,label = F)+NoLegend()+
  labs(title="", x="Cells separated by clusters", y = "5 most-sig.DEGs in each cluster",size=40)+
  theme(axis.text.y = element_blank())
dev.off()
### 每个簇中的标记基因
marker_gene <- all.markers%>%group_by(cluster)%>%top_n(n=1,wt=avg_log2FC)
marker_gene <- marker_gene$gene[!duplicated(marker_gene$gene)]

pdf(file = '07.marker.exp.pdf',w=8,h=6)
DotPlot(UMAP, features = marker_gene,group.by = 'seurat_clusters') + RotatedAxis()
dev.off()
png(file = '07.marker.exp.png',w=600,h=450)
DotPlot(UMAP, features = marker_gene,group.by = 'seurat_clusters') + RotatedAxis()
dev.off()

setwd("/data/nas1/luchunlin/project/BJTC-406-12/16_scRNA/")
if (! dir.exists("./03_SingleR")){
  dir.create("./03_SingleR")
}
setwd("./03_SingleR")
## PART E 细胞类型注释------
###SingleR------
library(SingleR)
scRNA_nor_singleR <- GetAssayData(UMAP, slot = "data")
hpca.se <- HumanPrimaryCellAtlasData()
hpca.se
clusters=UMAP@meta.data$seurat_clusters
#clusters=TSNE@meta.data$seurat_clusters
pred.hesc <- SingleR(test = scRNA_nor_singleR,
                     ref = hpca.se,
                     labels = hpca.se$label.main,
                     method = "cluster", 
                     clusters = clusters, 
                     quantile = 0.2,
                     assay.type.test = "logcounts", 
                     assay.type.ref = "logcounts")

table(pred.hesc$labels)
#9
celltype = data.frame(ClusterID=rownames(pred.hesc), 
                      celltype=pred.hesc$labels, 
                      stringsAsFactors = F) 
UMAP@meta.data$singleR=celltype[match(clusters,celltype$ClusterID),'celltype']
p6 <- DimPlot(UMAP, reduction = "umap", group.by = "singleR", label = T)
p6
phe=UMAP@meta.data
table(phe$singleR)
cell_type_stat <- as.data.frame(sort(table(phe$singleR)))
write.table(cell_type_stat,
            file = "cell_type_stat.txt",
            quote = F)
library(Ipaper)
write_fig(p6,
          file = "01.cell_type_singleR.pdf",
          width = 10,
          height = 8,
          devices = NULL,
          res = 600,
          show = F)
write_fig(p6,
          file = "01.cell_type_singleR.png",
          width = 10,
          height = 8,
          devices = NULL,
          res = 300,
          show = F)

mac.gene <- all.markers[which(all.markers$cluster==8 | all.markers$cluster == 11),]
mac.gene <- mac.gene[which(mac.gene$avg_log2FC>0.7),]
length(unique(mac.gene$gene))
##424
write.table(mac.gene,file = 'mac.gene.xls',sep = '\t',row.names = T,quote = F)
new.cluster.ids = as.character(celltype$celltype)
new.cluster.ids
names(new.cluster.ids) <- celltype$ClusterID

UMAP <- RenameIdents(UMAP,new.cluster.ids)
# 
# marker_gene <- cbind(new.cluster.ids,marker_gene)
# marker_gene <- as.data.frame(marker_gene)
# marker_gene <- marker_gene[!duplicated(marker_gene$new.cluster.ids),]
# marker_gene <- marker_gene$marker_gene
# #DotPlot(UMAP, features = marker_gene,group.by = 'singleR') + RotatedAxis()
# p7 <- VlnPlot(UMAP, features = marker_gene, group.by = "singleR", ncol = 4)
# p7
# ggsave(filename = '02.hub.violin.png',p7,width = 10,height = 9)
# ggsave(filename = '02.hub.violin.pdf',p7,width = 10,height = 9)

setwd("/data/nas1/luchunlin/project/BJTC-406-12/16_scRNA/")
if (! dir.exists("./04_Cellchat")){
  dir.create("./04_Cellchat")
}
setwd("./04_Cellchat")
# 细胞通讯--------
saveRDS(UMAP,file = 'UMAP.rds')
#安装依赖包：
devtools::install_github("jokergoo/ComplexHeatmap")
devtools::install_github("sqjin/NMF")
devtools::install_github("jokergoo/circlize")
install.packages('reshape2')
#pip install umap-learn

#python中安装umap-learn
#pip install umap-learn
#安装cellchat
devtools::install_github("sqjin/CellChat")
library(CellChat) 
#读入seurat处理后的rds文件
scRNA<-readRDS(file = 'scRNA.rds')
#data.input = as.data.frame(GetAssayData(subset(scRNA), slot='counts'))# raw count
data.input<-scRNA@assays$RNA@data # normalized data matrix
meta = scRNA@meta.data # a dataframe with rownames containing cell mata data
if (!dir.exists("cellchat")) {dir.create("cellchat")}
setwd('cellchat')
save(data.input,meta,file='01.inputdata.rda')
# table(meta$Type) ## AMLL 15058    7337 Healthy 
# cell.use = rownames(meta)[meta$Type == "AML"] # 提取正常组织
# data.input = data.input[, cell.use]
# meta = meta[cell.use, ]
unique(meta$singleR)

##CellChat的输入需要的matrix和meta我们已经准备好，下面开始创建
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "singleR")
cellchat <- addMeta(cellchat, meta = meta)##增加其他meta信息
cellchat <- setIdent(cellchat, ident.use = "singleR") # 将 "labels" 设为默认细胞标记类型，这个可以根据自己的数据自定义
levels(cellchat@idents)
groupSize <- as.numeric(table(cellchat@idents)) # 每组细胞的数量

##基于配受体分析的数据库
CellChatDB <- CellChatDB.human # 包括人和老鼠的，这里我们用human
showDatabaseCategory(CellChatDB)###作者提供了可视化的代码，可以看到该数据库中“Secreted Signaling”占比过半
library(tidyverse)
dplyr::glimpse(CellChatDB$interaction)  ###看一下CellChatDB的基本结构
CellChatDB_interaction <- CellChatDB$interaction
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # 我们这里使用“Secreted Signaling”部分做后续的细胞通讯分析
cellchat@DB <- CellChatDB.use

##表达数据做进一步预处理,节省算力
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)###首先识别过表达基因（配体——受体）
cellchat <- identifyOverExpressedInteractions(cellchat)###然后识别过表达配受体之间过表达的相互作用   （绕绕绕绕绕....）
# project gene expression data onto PPI network (optional)
cellchat <- projectData(cellchat, PPI.human)

##计算胞间通讯概率，预测通讯网络
cellchat <- computeCommunProb(cellchat)#
#cellchat <- filterCommunication(cellchat, min.cells = 10)##过滤掉小于10个细胞的胞间通讯网络

##胞间通讯网络的输出代码
df.net <- subsetCommunication(cellchat)  
write.table(df.net,'01.df.net.xls',sep = '\t',row.names = F,quote = F)
##信号通路的水平进一步推测胞间通讯，计算聚合网络
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
##可视化细胞互作的结果
groupSize <- as.numeric(table(cellchat@idents))
## 相互作用数目
pdf('01.number_of_interactions.pdf',w=6,h=5,family='serif')
netVisual_heatmap(cellchat)
dev.off()
png('01.number_of_interactions.png',w=500,h=400,family='serif')
netVisual_heatmap(cellchat)
dev.off()


setwd("/data/nas1/luchunlin/project/BJTC-406-12/16_scRNA/")
if (! dir.exists("./05_subtype")){
  dir.create("./05_subtype")
}
setwd("./05_subtype")
# Step.F 巨噬细胞亚群注释 ####
UMAP.som = UMAP[,UMAP@meta.data$singleR%in%c('Macrophage')]
sce=UMAP.som
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
sce <- FindVariableFeatures(sce, selection.method = 'vst', nfeatures = 2000)
sce <- ScaleData(sce)
sce <- RunPCA(sce, features = VariableFeatures(object = sce))

sce <- JackStraw(sce, num.replicate = 100, dims = 50)
sce <- ScoreJackStraw(sce, dims = 1:50)
JackStrawPlot(sce, dims = 1:50)
ElbowPlot(sce, ndims = 50)

sce <- FindNeighbors(sce, dims = 1:28)
sce <- FindClusters(sce, resolution = 0.4)

DimPlot(sce, reduction = "umap", group.by = "seurat_clusters", label = T)

# pdf(file = '01.subtype.pdf',w=6,h=5)
# DimPlot(sce, reduction = "umap", group.by = "seurat_clusters", label = T)
# dev.off() 
# png(file = '01.subtype.png',w=500,h=400)
# DimPlot(sce, reduction = "umap", group.by = "seurat_clusters", label = T)
# dev.off() 
##注释
library(SingleR)
scRNA.som_nor_singleR <- GetAssayData(sce, slot = "data")
sub.clusters=sce@meta.data$seurat_clusters
#clusters=TSNE@meta.data$seurat_clusters
sub.pred.hesc <- SingleR(test = scRNA.som_nor_singleR,
                     ref = hpca.se,
                     labels = hpca.se$label.fine,
                     method = "cluster", 
                     clusters = sub.clusters, 
                     quantile = 0.4,
                     assay.type.test = "logcounts", 
                     assay.type.ref = "logcounts")

table(sub.pred.hesc$labels)

subcelltype = data.frame(ClusterID=rownames(sub.pred.hesc), 
                      celltype=sub.pred.hesc$labels, 
                      stringsAsFactors = F) 
sce@meta.data$singleR=subcelltype[match(sub.clusters,subcelltype$ClusterID),'celltype']
subp <- DimPlot(sce, reduction = "umap", group.by = "singleR", label = F)
subp

sub.phe=sce@meta.data
table(sub.phe$singleR)
subcell_type_stat <- as.data.frame(sort(table(sub.phe$singleR)))
write.table(subcell_type_stat,
            file = "subcell_type_stat.txt",
            quote = F)

write_fig(subp,
          file = "01.subcell_type_singleR.pdf",
          width = 10,
          height = 8,
          devices = NULL,
          res = 600,
          show = F)
write_fig(subp,
          file = "01.subcell_type_singleR.png",
          width = 10,
          height = 8,
          devices = NULL,
          res = 300,
          show = F)
# PRAR G Trajectory------
setwd("/data/nas1/luchunlin/project/BJTC-406-12/16_scRNA/")
if (! dir.exists("./06_Trajectory")){
  dir.create("./06_Trajectory")
}
setwd("./06_Trajectory")
library(monocle)

Mono_tj = sce
Mono_matrix = GetAssayData(Mono_tj, slot = "count", assay = "RNA")
feature_ann = data.frame(gene_id=rownames(Mono_matrix),gene_short_name=rownames(Mono_matrix))
rownames(feature_ann) = rownames(Mono_matrix)
Mono_fd = new("AnnotatedDataFrame", data = feature_ann)
sample_ann = Mono_tj@meta.data
cell.type = Idents(Mono_tj) %>% data.frame
sample_ann = cbind(sample_ann, cell.type)
colnames(sample_ann)[ncol(sample_ann)] = "Cell_Type"
Mono_pd = new("AnnotatedDataFrame", data =sample_ann)
Mono.cds = newCellDataSet(Mono_matrix, phenoData =Mono_pd, featureData=Mono_fd, expressionFamily=negbinomial.size())
rm(Mono_tj)
Mono.cds = estimateSizeFactors(Mono.cds)
Mono.cds = estimateDispersions(Mono.cds)
WGCNA::collectGarbage()

#plot_ordering_genes(Mono.cds)
#8 9+3 6 1+2
Mono.cds = reduceDimension(Mono.cds, max_components = 2, verbose = T, 
                           norm_method = "log"#, residualModelFormulaStr = "~nFeature_RNA+group"
                           )
Mono.cds = orderCells(Mono.cds)
save(Mono.cds, file = "mycds_order.RData")
table(Mono.cds$singleR)
Mono.cds$singleR = factor(Mono.cds$singleR, levels = c("Macrophage:Alveolar","Macrophage:Alveolar:B._anthacis_spores",
                                                           "Macrophage:monocyte-derived:M-CSF","NK_cell",
                                                           "NK_cell:CD56hiCD62L+"))
Mono.cds$singleR

p = plot_cell_trajectory(Mono.cds,color_by="Pseudotime",cell_size = 1, theta = 180,
                         size=1, show_backbone=TRUE, show_branch_points = F) +
 #  scale_color_manual(values = c("red","orange","green","blue","yellow",'brown','darkgreen','midnightblue','pink','purple')) +
  theme(text = element_text(size = 14))
p
ggsave("01.PSD_Trajectory.png", width = 6, height = 6, units = "in", dpi = 300)
ggsave("01.PSD_Trajectory.pdf", width = 6, height = 6, units = "in", dpi = 300)
# saveRDS(Mono.cds, "mono.RDS")
p2 = plot_cell_trajectory(Mono.cds,color_by="singleR",cell_size = 1, theta = 180,
                          size=1, show_backbone=TRUE, show_branch_points = F) +
  scale_color_manual(values = c("red","orange","blue",'darkgreen','purple')) +
  theme(text = element_text(size = 14))
p2
ggsave("02.cell_type_Trajectory.png", width = 6, height = 6, units = "in", dpi = 300)
ggsave("02.cell_type_Trajectory.pdf", width = 6, height = 6, units = "in", dpi = 300)

###预后基因表达------
lasso.gene<-read.delim2('/data/nas1/luchunlin/project/BJTC-406-12/08_Lasso/lasso_genes.csv',header = F)
hubgene<-lasso.gene$V1
sce$singleR



p3 <- VlnPlot(sce, features = hubgene, group.by = "singleR", ncol = 4)
p3 <- VlnPlot(sce, features = hubgene, group.by = "singleR", ncol = 4)
p3
# ggsave(filename = '03.hub.violin.png',p3,width = 12,height = 10)
# ggsave(filename = '03.hub.violin.pdf',p3,width = 12,height = 10)


##一致性聚类-------
setwd("/data/nas1/luchunlin/project/BJTC-406-12/16_scRNA/")
if (! dir.exists("./07_consensus")){
  dir.create("./07_consensus")
}
setwd("./07_consensus")
##提取关键细胞簇
##3个巨噬细胞相关的
table(sce$singleR)
UMAP.cluster = sce[,sce@meta.data$singleR%in%c('Macrophage:Alveolar','Macrophage:Alveolar:B._anthacis_spores','Macrophage:monocyte-derived:M-CSF')]
saveRDS(UMAP.cluster,file = 'UMAP.cluster.rds')
dat.cluster <- UMAP.cluster@assays$RNA@data%>%as.data.frame()
tpm_dat <- readRDS("/data/nas1/luchunlin/project/BJTC-406-12/16_scRNA/GSE131907_Lung_Cancer_normalized_log2TPM_matrix.rds")
View(tpm_dat[1:3,1:3])
tpm_dat <- tpm_dat[,colnames(tpm_dat)%in%colnames(dat.cluster)]
tpm_dat <- tpm_dat[rownames(tpm_dat)%in%lasso.gene$V1,]


#dat.cluster<-log2(dat.cluster+1)
#cluster_exp<-log2(cluster_exp+0.0001)
cluster_exp<-as.matrix(tpm_dat)
library(ConsensusClusterPlus)
#sweep函数减去中位数进行标准化
df<-sweep(cluster_exp,1, apply(cluster_exp,1,median,na.rm=T))
maxK <-  7 #最多分成几组
results <-  ConsensusClusterPlus(df,
                                 maxK = maxK,
                                 reps = 1000,              # 抽样次数(一般1000或更多)
                                 pItem = 0.8,             # 抽样比例
                                 pFeature = 1,
                                 clusterAlg = "pam",      # 聚类方法
                                 seed = 123,
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
## [2]
cluster<-results[[2]]$consensusClass
cluster
group <- data.frame(cluster)%>%rownames_to_column(var = 'sample')
group$cluster <- ifelse(group$cluster=='1','cluster 1','cluster 2')

write.table(group,file = 'cluster.xls',sep = '\t',row.names = F,quote = F)

###特征基因表达水平------
hub_exp<-cluster_exp
hub_exp2<-hub_exp%>%as.data.frame()
hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))

control.sample <- group$sample[which(group$cluster=='cluster 1')]
## 样本分组
hub_exp2$Group<-ifelse(hub_exp2$sample%in%control.sample,'cluster 1','cluster 2')

##分面图形
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(rstatix)
stat.test<-hub_exp2%>%
  group_by(Symbol)%>%
  wilcox_test(expr ~ Group)%>%
  adjust_pvalue(method = 'fdr')

stat.test$p.adj<-round(stat.test$p.adj,digits = 3)
stat.test$p.adj<-ifelse(stat.test$p.adj<0.001,"***",
                        ifelse(stat.test$p.adj<0.05,"**",
                               ifelse(stat.test$p.adj<0.05,"*",'ns')))
exp_plot <- ggplot(hub_exp2,aes(x = Group, y = expr, color = Group)) +
  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.4,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#CD3700","#4682B4"), name = "Group")+
  labs(title="Expression level", x="", y = "",size=20) +
  stat_pvalue_manual(stat.test,
                     y.position = c(7,7,7,5,4,8,4),
                     size = 3.2,
                     family = "Times",
                     label = "p.adj",
                     #parse = T,
                     face = "bold")+
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=15),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=12), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  guides(fill='none')+
  facet_wrap(~Symbol,scales = "free",nrow = 2) 
exp_plot
ggsave('02.expression.pdf',exp_plot,width = 7,height = 6)
ggsave('02.expression.png',exp_plot,width = 7,height = 6)

###GSEA------
library(GSVA)
library(GSEABase)
library(DESeq2)
group.gsea <- group[order(group$cluster),]
gsea_exp<-dat.cluster[,group$sample]%>%round(digits = 0)
library(DESeq2)
colData<-group.gsea
colnames(colData)<-c('sample','group')
table(colData$group)
colData$group<-factor(colData$group,levels = c('cluster 1','cluster 2'))
dds<-DESeqDataSetFromMatrix(countData = gsea_exp,colData=colData,design = ~group)
dds = dds[rownames(counts(dds)) > 10,]
dds<-DESeq(dds)
dds
res =results(dds, contrast = c("group","cluster 1","cluster 2"))
res =res[order(res$padj),]
head(res)
summary(res)
table(res$padj<0.05)
allGeneSets<-as.data.frame(res)
allGeneSets<-na.omit(allGeneSets)
logFCcutoff <- 0
allGeneSets$change = as.factor(
  ifelse(allGeneSets$padj < 0.05 & abs(allGeneSets$log2FoldChange) > logFCcutoff,
         ifelse(allGeneSets$log2FoldChange > logFCcutoff,'UP','DOWN'),'NOT')
)

DEGeneSets <- subset(allGeneSets,
                     allGeneSets$padj < 0.05 & abs(allGeneSets$log2FoldChange) > 0)
genelist <- allGeneSets$log2FoldChange
names(genelist) <- rownames(allGeneSets)
geneList <- sort(genelist, decreasing = T)
# DEGeneSets <- DEGeneSets[order(DEGeneSets$padj),]
# dim(DEGeneSets)
# write.table(DEGeneSets,'DEGeneSet.xls',sep = '\t',row.names = T,quote = F)
# write.table(allGeneSets,'allGeneSet.xls',sep = '\t',row.names = T,quote = F)
## GSEA KEGG-----
library(clusterProfiler)
library(enrichplot)
kegg_set<- read.gmt("/data/nas1/luchunlin/pipeline/GSVA/c2.cp.kegg.v7.4.symbols.gmt")
set.seed(1)
kegg_gsea <- GSEA(geneList, TERM2GENE = kegg_set, pvalueCutoff = 0.25)
kegg_result <- kegg_gsea@result
write.table(kegg_result,file = 'GSEA(KEGG).result.xls',sep = '\t',row.names = F,quote = F)
dim(kegg_result)
pdf(file = '03.GSEA(KEGG).pdf',w=9,h=6)
gseaplot2(kegg_gsea,c(1:5),color = c('#7B68EE','#CD3333','#20B2AA','#FF8C00','#FF6666'),
          title = 'KEGG GSEA',
          base_size = 11,
          rel_heights = c(1.5, 0.3, 0.5))
dev.off()
png(file = '03.GSEA(KEGG).png',w=600,h=400)
gseaplot2(kegg_gsea,c(1:5),color = c('#7B68EE','#CD3333','#20B2AA','#FF8C00','#FF6666'),
          title = 'KEGG GSEA',
          base_size = 11,
          rel_heights = c(1.5, 0.3, 0.5))
dev.off()



go_set<- read.gmt("/data/nas1/luchunlin/pipeline/GSVA/c5.go.bp.v7.4.symbols.gmt")
set.seed(1)
go_gsea <- GSEA(geneList, TERM2GENE = go_set, pvalueCutoff = 0.3)
go_result <- go_gsea@result
write.table(go_result,file = 'GSEA(GO).result.xls',sep = '\t',row.names = F,quote = F)
dim(go_result)
pdf(file = '03.GSEA(GO).pdf',w=9,h=6)
gseaplot2(go_gsea,c(1:5),color = c('#7B68EE','#CD3333','#20B2AA','#FF8C00','#FF6666'),
          title = 'GO GSEA',
          base_size = 11,
          rel_heights = c(1.5, 0.3, 0.5))
dev.off()
png(file = '03.GSEA(GO).png',w=600,h=400)
gseaplot2(go_gsea,c(1:5),color = c('#7B68EE','#CD3333','#20B2AA','#FF8C00','#FF6666'),
          title = 'GO GSEA',
          base_size = 11,
          rel_heights = c(1.5, 0.3, 0.5))
dev.off()


##pyscenic
##cluster
pbmc.all<-t(as.matrix(UMAP.cluster@assays$RNA@counts))
write.csv(pbmc.all,file = 'pbmc.all.csv')
# pbmc.all2<-as.matrix(UMAP.cluster@assays$RNA@counts)
# write.csv(pbmc.all2,file = 'pbmc2.all.csv')
