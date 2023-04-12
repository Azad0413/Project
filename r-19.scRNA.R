rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/HF-0103-1/")
if (! dir.exists("./19_scRNA")){
  dir.create("./19_scRNA")
}
setwd("./19_scRNA")
setwd("/data/nas1/luchunlin/project/HF-0103-1/19_scRNA/")
if (! dir.exists("./00_rawdata")){
  dir.create("./00_rawdata")
}
setwd("./00_rawdata")
library(Seurat)
library(GEOquery)
library(data.table)
library(tidyverse)
library(lance)
## PART A导入数据--------
## 导入数据，创建seurat对象
dat <- fread('GSE182434_raw_count_matrix.txt.gz',header = T)%>%
  column_to_rownames(var = 'Gene')
colnames(dat)
##
smps <- data.frame(cell=colnames(dat),sample=colnames(dat) %>% str_sub(18,26))
table(smps$sample)
##保留DLBCL
sample.list = c('DLBCL002B','DLBCL002N','DLBCL007B','DLBCL007N','DLBCL008B','DLBCL008N','DLBCL111B','DLBCL111N')
sample.final <- smps[smps$sample%in%sample.list,]
dat <- dat[,colnames(dat)%in%sample.final$cell]

scRNA <- CreateSeuratObject(counts = dat,min.cells = 3,min.features = 200)
##初步过滤，<=3个细胞中表达的基因，>=200个基因的细胞
scRNA
system.time(save(scRNA, file = "01.scRNA_orig.Rdata"))  

##QC------
rm(list = ls())
setwd(('/data/nas1/luchunlin/project/HF-0103-1/19_scRNA/'))
if (! dir.exists("./01_QC/")){
  dir.create("01_QC")
}
setwd('../01_QC/')
### 质控前 ###
load('../19_scRNA/00_rawdata/01.scRNA_orig.Rdata')
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
                subset = nFeature_RNA > 200 & nFeature_RNA < 6000&
                  percent.mt < 5)
length(colnames(scRNA))#11934细胞数量
length(rownames(scRNA))# 21619基因数量
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

system.time(save(scRNA, file = "scRNA_qc.Rdata"))  


rm(list = ls())
## PART C数据降维--------
setwd(('/data/nas1/luchunlin/project/HF-0103-1//19_scRNA/'))
if (! dir.exists("./02_PCA")){dir.create("./02_PCA")}
setwd("./02_PCA")

load('../01_QC/scRNA_qc.Rdata')
length(colnames(scRNA))
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
all.markers<-FindAllMarkers(UMAP,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 1,test.use = 'wilcox',return.thresh = 0.05)
all.markers%>%group_by(cluster)%>%top_n(n=2,wt=avg_log2FC)

write.table(all.markers,file = 'All.markers.xls',sep = '\t',row.names = F,quote = F)
top10<-all.markers%>%group_by(cluster)%>%top_n(n=5,wt=avg_log2FC)
# pdf('06.topheat.pdf',w=5,h=4)
# DoHeatmap(UMAP,features = top10$gene,label = F)+NoLegend()+
#   labs(title="", x="Cells separated by clusters", y = "5 most-sig.DEGs in each cluster",size=40)+
#   theme(axis.text.y = element_blank())
# dev.off()
# png('06.topheat.png',w=600,h=500)
# DoHeatmap(UMAP,features = top10$gene,label = F)+NoLegend()+
#   labs(title="", x="Cells separated by clusters", y = "5 most-sig.DEGs in each cluster",size=40)+
#   theme(axis.text.y = element_blank())
# dev.off()
### 每个簇中的标记基因
marker_gene <- all.markers%>%group_by(cluster)%>%top_n(n=1,wt=avg_log2FC)
marker_gene <- marker_gene$gene[!duplicated(marker_gene$gene)]

# pdf(file = '07.marker.exp.pdf',w=8,h=6)
# DotPlot(UMAP, features = marker_gene,group.by = 'seurat_clusters') + RotatedAxis()
# dev.off()
# png(file = '07.marker.exp.png',w=600,h=450)
# DotPlot(UMAP, features = marker_gene,group.by = 'seurat_clusters') + RotatedAxis()
# dev.off()

setwd("/data/nas1/luchunlin/project/HF-0103-1/19_scRNA/")
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
                     quantile = 0.8,
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

pdf(file = '02.marker.exp.pdf',w=8,h=6)
DotPlot(UMAP, features = marker_gene,group.by = 'singleR') + RotatedAxis()
dev.off()
png(file = '02.marker.exp.png',w=600,h=450)
DotPlot(UMAP, features = marker_gene,group.by = 'singleR') + RotatedAxis()
dev.off()


hubgene <- read.delim2('/data/nas1/luchunlin/project/HF-0103-1/06_Lasso/lasso_genes.csv',header = F)
pdf(file = '03.expression.pdf',w=8,h=5)
FeaturePlot(UMAP,features = hubgene$V1,ncol = 3)
dev.off()

png(file = '03.expression.png',w=650,h=400)
FeaturePlot(UMAP,features = hubgene$V1,ncol = 3)
dev.off()

## PART F 细胞轨迹分析--------
setwd("/data/nas1/luchunlin/project/HF-0103-1/19_scRNA/")
if (! dir.exists("./04_trajectory")){
  dir.create("./04_trajectory")
}
setwd("./04_trajectory")

library(future)
library(monocle3)
library(monocle)
data <- as(as.matrix(UMAP@assays$RNA@counts), 'sparseMatrix')
# count矩阵
pd <- new('AnnotatedDataFrame', data = UMAP@meta.data)
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
sce.markers <- FindAllMarkers(UMAP)
all.markers = sce.markers %>% dplyr::select(gene, everything()) %>% 
  subset(p_val<0.05 & abs(sce.markers$avg_log2FC) > 0.5)
# 
# write.table(all.markers, file = "All_Markers.xls",
#             quote = F, row.names = F)
markers.gene <- all.markers$gene
mycds <- setOrderingFilter(mycds, markers.gene)
#降维
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
save(mycds, file = "mycds_reduce.RData")
#耗时，耗内存
#排序
#load("mycds_reduce.RData")
mycds <- orderCells(mycds)
save(mycds, file = "mycds_order.RData")
load("mycds_order.RData")
table(mycds$singleR)
mycds$singleR <- factor(mycds$singleR, levels = c("B_cell",
                                                  "DC" ,
                                                  "Monocyte",
                                                  "NK_cell",
                                                  "T_cells"))
# pdf("sub_cell_trajectory_raw.pdf",w=10,h=8)
# plot_cell_trajectory(mycds, color_by = "seurat_clusters") +
#   theme(legend.position = "right") +
#   theme(axis.title.x =element_text(size=22,color='black', face = "bold",family='Times'),
#         axis.text.x =element_text(size=18, family='Times'),
#         axis.title.y =element_text(size=22,color='black', face = "bold",family='Times'),
#         axis.text.y=element_text(size=18,  family='Times'),
#         legend.title=element_text(size=20, color='black', face = "bold",family='Times'),
#         legend.text=element_text(size=18,  face = "bold",family='Times'),
#         title=element_text(size=20, color='black', face = "bold",family='Times'))+
#   theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
# dev.off()
# 
# png("sub_cell_trajectory_raw.png",w=10,h=8,units='in',res=600)
# plot_cell_trajectory(mycds, color_by = "seurat_clusters") +
#   theme(legend.position = "right") +
#   theme(axis.title.x =element_text(size=22,color='black', face = "bold",family='Times'),
#         axis.text.x =element_text(size=18, family='Times'),
#         axis.title.y =element_text(size=22,color='black', face = "bold",family='Times'),
#         axis.text.y=element_text(size=18,  family='Times'),
#         legend.title=element_text(size=20, color='black', face = "bold",family='Times'),
#         legend.text=element_text(size=18,  face = "bold",family='Times'),
#         title=element_text(size=20, color='black', face = "bold",family='Times'))+
#   theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
# dev.off()

pdf('01.sub_cell_trajectory_singleR.pdf',w=10,h=8)
plot_cell_trajectory(mycds, color_by = "singleR") +
  scale_color_discrete(breaks = c("B_cell",
                                  "DC" ,
                                  "Monocyte",
                                  "NK_cell",
                                  "T_cells" 
  ))+
  theme(legend.position = "right") +
  theme(axis.title.x =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.x =element_text(size=18, family='Times'),
        axis.title.y =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.y=element_text(size=18,  family='Times'),
        legend.title=element_text(size=20, color='black', face = "bold",family='Times'),
        legend.text=element_text(size=18,  face = "bold",family='Times'),
        title=element_text(size=20, color='black', face = "bold",family='Times'))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
dev.off()

png('01.sub_cell_trajectory_singleR.png',w=10,h=8,units='in',res=600)
plot_cell_trajectory(mycds, color_by = "singleR") +
  scale_color_discrete(breaks = c("B_cell",
                                  "DC" ,
                                  "Monocyte",
                                  "NK_cell",
                                  "T_cells" 
  ))+
  theme(legend.position = "right") +
  theme(axis.title.x =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.x =element_text(size=18, family='Times'),
        axis.title.y =element_text(size=22,color='black', face = "bold",family='Times'),
        axis.text.y=element_text(size=18,  family='Times'),
        legend.title=element_text(size=20, color='black', face = "bold",family='Times'),
        legend.text=element_text(size=18,  face = "bold",family='Times'),
        title=element_text(size=20, color='black', face = "bold",family='Times'))+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
dev.off()
p = plot_cell_trajectory(mycds,color_by="Pseudotime",cell_size = 1, theta = 180,
                         size=1, show_backbone=TRUE, show_branch_points = F) +
  #  scale_color_manual(values = c("red","orange","green","blue","yellow",'brown','darkgreen','midnightblue','pink','purple')) +
  theme(text = element_text(size = 14))
p
ggsave("01.PSD_Trajectory.png", width = 6, height = 6, units = "in", dpi = 300)
ggsave("01.PSD_Trajectory.pdf", width = 6, height = 6, units = "in", dpi = 300)

p
ggsave("01.PSD_Trajectory.png", width = 6, height = 6, units = "in", dpi = 300)
ggsave("01.PSD_Trajectory.pdf", width = 6, height = 6, units = "in", dpi = 300)
# saveRDS(Mono.cds, "mono.RDS")
p2 = plot_cell_trajectory(mycds,color_by="singleR",cell_size = 1, theta = 180,
                          size=1, show_backbone=TRUE, show_branch_points = F) +
  scale_color_manual(values = c("red","orange","blue",'darkgreen','purple')) +
  theme(text = element_text(size = 14))
p2
ggsave("02.cell_type_Trajectory.png", width = 6, height = 6, units = "in", dpi = 300)
ggsave("02.cell_type_Trajectory.pdf", width = 6, height = 6, units = "in", dpi = 300)



plot16 <- plot_genes_in_pseudotime(mycds[hubgene$V1,],
                                   color_by = "singleR",
                                   ncol = 3)+ 
  theme(axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))+
  scale_color_manual(values = c("red","orange","blue",'darkgreen','purple')) 
plot16


write_fig(plot16,
          file = "03.pseudotim.gene.pdf",
          width = 8,
          height = 6,
          devices = NULL,
          res = 600,
          show = F)
write_fig(plot16,
          file = "03.pseudotime.gene.png",
          width = 8,
          height = 6,
          devices = NULL,
          res = 300,
          show = F)




