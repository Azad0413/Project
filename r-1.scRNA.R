rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/JNZK-205/")
if (! dir.exists("./01_scRNA")){
  dir.create("./01_scRNA")
}
setwd("./01_scRNA")
library(Seurat)
library(GEOquery)
## 导入数据，创建seurat对象
## GSE158937

dir='GSE158937_RAW/'
list.files(dir)
GSM4816045<-Read10X(data.dir ='GSE158937_RAW/GSM4816045_RAW/')
GSM4816046<-Read10X(data.dir ='GSE158937_RAW/GSM4816046_RAW/')
GSM4816047<-Read10X(data.dir ='GSE158937_RAW/GSM4816047_RAW/')
count1<-CreateSeuratObject(counts = GSM4816045,min.cells = 3,min.features = 200)
count2<-CreateSeuratObject(counts = GSM4816046,min.cells = 3,min.features = 200)
count3<-CreateSeuratObject(counts = GSM4816047,min.cells = 3,min.features = 200)

scRNA<-merge(x=count1,y=c(count2,count3),
             add.cell.ids=c('GSM4816045','GSM4816046','GSM4816047'),
             merge.data=T)
##初步过滤，<=3个细胞中表达的基因，>=200个基因的细胞
scRNA
## 质控------
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
                subset = nFeature_RNA > 200 & nFeature_RNA < 2500 &
                  percent.mt < 5)

VlnPlot(scRNA,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
pdf(file = '01.QC.pdf',w=8,h=6)
VlnPlot(scRNA,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
dev.off()
png(file = '01.QC.png',w=800,h=600)
VlnPlot(scRNA,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
dev.off()

#标准化
scRNA.norm<-NormalizeData(scRNA,normalization.method = "LogNormalize",scale.factor = 10000)
##scRNA.norm<-scRNA
## PCA UMAP降维聚类-------
## 缩放数据
## 筛选1500个高变基因
all.genes<-rownames(scRNA.norm)
scRNA.norm<-FindVariableFeatures(scRNA.norm,selection.method = "vst", nfeatures = 2000)
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
          file = "02.feature_selection.pdf",
          width = 5,
          height = 6,
          devices = NULL,
          res = 600,
          show = F)
write_fig(plot4,
          file = "02.feature_selection.png",
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

# write_fig(plot5,
#           file = "pca_cluster.pdf",
#           width = 15,
#           height = 8,
#           devices = NULL,
#           res = 600,
#           show = F)
# write_fig(plot5,
#           file = "pca_cluster.png",
#           width = 15,
#           height = 8,
#           devices = NULL,
#           res = 300,
#           show = F)
plot6 <- ElbowPlot(scRNA.norm.pca, ndims = 50)
plot6
# write_fig(plot6,
#           file = "pca_sd.pdf",
#           width = 15,
#           height = 8,
#           devices = NULL,
#           res = 600,
#           show = F)
# write_fig(plot6,
#           file = "pca_sd.png",
#           width = 15,
#           height = 8,
#           devices = NULL,
#           res = 300,
#           show = F)
scRNA.norm.pca.c<-FindNeighbors(scRNA.norm.pca,dims = 1:27)
scRNA.norm.pca.c<-FindClusters(scRNA.norm.pca.c,resolution = 0.4)
UMAP<-RunUMAP(scRNA.norm.pca.c,dims = 1:27)
DimPlot(UMAP,reduction = 'umap')

#-tsne-----
# set.seed(123)
# TSNE <- RunTSNE(scRNA.norm.pca.c, dims = 1:27)
# head(TSNE@reductions$tsne@cell.embeddings)
# p5 <- DimPlot(TSNE, reduction = "tsne", label = T)
# p5
# saveRDS(TSNE,file = 'TSNE.rds')
## 细胞类型注释------
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
                    # quantile = 0.7,
                     assay.type.test = "logcounts", 
                     assay.type.ref = "logcounts")
table(pred.hesc$labels)
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

write_fig(p6,
          file = "03.cell_type_singleR.pdf",
          width = 10,
          height = 8,
          devices = NULL,
          res = 600,
          show = F)
write_fig(p6,
          file = "03.cell_type_singleR.png",
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

## 细胞轨迹分析------
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








