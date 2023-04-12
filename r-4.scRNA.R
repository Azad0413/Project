rm(list = ls())
setwd("/data/nas1/luchunlin/project/JNZK-214-8/")
if (! dir.exists("./04_scRNA")){
  dir.create("./04_scRNA")
}
setwd("./04_scRNA")
library(Seurat)
library(GEOquery)
library(data.table)
## PART A导入数据--------
## 导入数据，创建seurat对象
dat <- fread('GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt',header = T)%>%
  column_to_rownames(var = 'Index')
scRNA <- CreateSeuratObject(counts = dat,min.cells = 3,min.features = 200)
##初步过滤，<=3个细胞中表达的基因，>=200个基因的细胞
scRNA
## ##PART B质控-------
setwd(('/data/nas1/luchunlin/project/JNZK-214-8/04_scRNA/'))
dir.create("01_QC")
setwd('../04_scRNA/01_QC/')
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
        features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 1) + NoLegend() #+theme.set 
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
length(colnames(scRNA))#31383细胞数量
length(rownames(scRNA))# 25655基因数量
VlnPlot(scRNA,
        features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 1) + NoLegend()

pdf(file = '01.QC.pdf',w=10,h=9)
VlnPlot(scRNA,
        features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 1) + NoLegend()
dev.off()
png(file = '01.QC.png',w=900,h=800)
VlnPlot(scRNA,
        features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 1) + NoLegend()

dev.off()

##添加分组
smps = Cells(scRNA) %>% str_sub(1,7)
table(smps)
grps = read.delim2('meta.txt',header = T)
group.label = grps$group[match(smps,grps$sample)]
names(group.label) = Cells(scRNA)
group.label
scRNA = AddMetaData(scRNA,group.label,col.name = 'group')
## PART C数据降维--------
setwd(('/data/nas1/luchunlin/project/JNZK-214-8/04_scRNA/'))
if (! dir.exists("./02_PCA")){dir.create("./02_PCA")}
setwd("./02_PCA")
#标准化
scRNA.norm<-NormalizeData(scRNA,normalization.method = "LogNormalize",scale.factor = 10000)
##scRNA.norm<-scRNA
## PART D PCA UMAP降维聚类-------
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
UMAP$group
DefaultAssay(UMAP) = "RNA"
marker.genes = lapply(0:21, function(x){
  FindConservedMarkers(UMAP, ident.1 = x, grouping.var = "group", 
                       logfc.threshold = 1, min.pct = 0.1, only.pos = T)
})
saveRDS(marker.genes, "marker.genes.rds")
names(marker.genes) = paste0("cluster", seq(0,21))
res.marker = lapply(marker.genes, function(x){
  cbind(gene = rownames(x), x)
})
#install.packages('WriteXLS')
library(WriteXLS)
WriteXLS(res.marker, "MarkerGenes.xls")

# cluster1.markers<-FindMarkers(UMAP,ident.1 = 1,min.pct = 0.25)
# head(cluster1.markers,n=5)
all.markers<-FindAllMarkers(UMAP,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 0.5,test.use = 'wilcox',return.thresh = 0.01)
all.markers%>%group_by(cluster)%>%top_n(n=2,wt=avg_log2FC)
# #write.table(all.markers,file = 'DEGs.xls',sep = '\t',row.names = F,quote = F)
# top10<-all.markers%>%group_by(cluster)%>%top_n(n=5,wt=avg_log2FC)
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
marker_gene <- c('IL7R','GZMK','CXCL13','IL1B','CD79A','TNFRSF4','IGKC','TFF3',
                 'HSPA1B','IGLC3','XCL2','IGHG1','ADAMDEC1','TAGLN','PLVAP',
                 'STMN1','PLAC8','SPINK4','CRYAB','HCAR3','TPSB2','ITLN1')
pdf(file = '07.marker.exp.pdf',w=8,h=6)

DotPlot(UMAP, features = marker_gene,group.by = 'seurat_clusters') + RotatedAxis()
dev.off()

png(file = '07.marker.exp.png',w=600,h=450)
DotPlot(UMAP, features = marker_gene,group.by = 'seurat_clusters') + RotatedAxis()
dev.off()

setwd("/data/nas1/luchunlin/project/JNZK-214-8/04_scRNA/")
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
#10
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

new.cluster.ids = as.character(celltype$celltype)
new.cluster.ids
names(new.cluster.ids) <- celltype$ClusterID

UMAP <- RenameIdents(UMAP,new.cluster.ids)

marker_gene <- cbind(new.cluster.ids,marker_gene)
marker_gene <- as.data.frame(marker_gene)
marker_gene <- marker_gene[!duplicated(marker_gene$new.cluster.ids),]
marker_gene <- marker_gene$marker_gene
#DotPlot(UMAP, features = marker_gene,group.by = 'singleR') + RotatedAxis()
p7 <- VlnPlot(UMAP, features = marker_gene, group.by = "singleR", ncol = 4)
p7
ggsave(filename = '02.hub.violin.png',p7,width = 10,height = 9)
ggsave(filename = '02.hub.violin.pdf',p7,width = 10,height = 9)
setwd("/data/nas1/luchunlin/project/JNZK-214-8/04_scRNA/")
if (! dir.exists("./06_DEGs")){
  dir.create("./06_DEGs")
}
setwd("./06_DEGs")
# Step.6 DEG ####
immune.cells = c("B_cell","Monocyte","NK_cell","T_cells")
UMAP.som = subset(UMAP, idents = c("B_cell","Monocyte","NK_cell","T_cells"))
saveRDS(UMAP.som, "immune.RDS")
UMAP.som$cell.type = Idents(UMAP.som)
UMAP.som$cell.type.group = paste(UMAP.som$cell.type, UMAP.som$group, sep = "_")
Idents(UMAP.som) <- "cell.type.group"
all.genes = rownames(UMAP.som)
all.genes = str_replace_all(all.genes, fixed("."), "-")
library(org.Hs.eg.db)
df.annot = select(org.Hs.eg.db, all.genes, c("ENTREZID"), "SYMBOL")
df.annot = na.omit(df.annot)
all.genes = all.genes[all.genes %in% df.annot$SYMBOL]
UMAP.som = UMAP.som[all.genes,]
DefaultAssay(UMAP.som) = "RNA"
res.deg = lapply(immune.cells,function(x){
  FindMarkers(UMAP.som, ident.1 = paste0(x,"_control"), ident.2 = paste0(x,"_tumor"),
              logfc.threshold = 0.1, min.pct = 0.1)
})

names(res.deg) = immune.cells
res.deg = lapply(res.deg, function(x){
  cbind(gene = rownames(x), x)
})
library(WriteXLS)
WriteXLS(res.deg, "DEG.xls")
# 
# top.genes = lapply(res.deg, function(x){
#   c(x[order(x$avg_log2FC, decreasing = T),]$gene[1:5],
#     x[order(x$avg_log2FC, decreasing = F),]$gene[1:5])
# }) %>% unlist %>% unique

Idents(UMAP.som) = "cell.type"
res.deg$B_cell$gene
gene1 <- data.frame(symbol=res.deg$B_cell$gene)
gene2 <- data.frame(symbol=res.deg$Monocyte$gene)
gene3 <- data.frame(symbol=res.deg$NK_cell$gene)
gene4 <- data.frame(symbol=res.deg$T_cells$gene)

allDEG <- rbind(gene1,gene2,gene3,gene4)
allDEG <- allDEG[!duplicated(allDEG$symbol),]%>%as.data.frame()
colnames(allDEG) <- 'symbol'
write.table(allDEG,file = 'scDEGs.xls',sep = '\t',row.names = F,quote = F)
# PRAR F Trajectory------
setwd("/data/nas1/luchunlin/project/JNZK-214-8/04_scRNA/")
if (! dir.exists("./04_Trajectory")){
  dir.create("./04_Trajectory")
}
setwd("./04_Trajectory")
library(monocle)
table(UMAP$singleR)
Mono_tj = subset(UMAP, idents = c("B_cell","Endothelial_cells","Epithelial_cells","Monocyte",
                                  "Neurons","Neutrophils","NK_cell","Smooth_muscle_cells","T_cells","Tissue_stem_cells"))
Mono_matrix = GetAssayData(Mono_tj, slot = "count", assay = "RNA")
#genes = UMAP@assays$integrated@var.features
genes = lapply(marker.genes, rownames) %>% unlist %>% unique
#genes = intersect(genes, UMAP@assays$integrated@var.features)
Mono_matrix = Mono_matrix[genes,]
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
                           norm_method = "log", residualModelFormulaStr = "~nFeature_RNA+group")
Mono.cds = orderCells(Mono.cds)

Mono.cds$Cell_Type = factor(Mono.cds$Cell_Type, levels = c("B_cell","Endothelial_cells","Epithelial_cells","Monocyte",
                                                           "Neurons","Neutrophils","NK_cell","Smooth_muscle_cells","T_cells","Tissue_stem_cells"))
p = plot_cell_trajectory(Mono.cds,color_by="Pseudotime",cell_size = 1, theta = 180,
                    size=1, show_backbone=TRUE, show_branch_points = F) +
# scale_color_manual(values = c("red","orange","green","blue","yellow",'brown','darkgreen','midnightblue','pink','purple')) +
 theme(text = element_text(size = 14))
p
ggsave("01.PSD_Trajectory.png", width = 6, height = 6, units = "in", dpi = 300)
ggsave("01.PSD_Trajectory.pdf", width = 6, height = 6, units = "in", dpi = 300)
saveRDS(Mono.cds, "mono.RDS")
p2 = plot_cell_trajectory(Mono.cds,color_by="Cell_Type",cell_size = 1, theta = 180,
                         size=1, show_backbone=TRUE, show_branch_points = F) +
  scale_color_manual(values = c("red","orange","green","blue","yellow",'brown','darkgreen','midnightblue','pink','purple')) +
  theme(text = element_text(size = 14))
p2
ggsave("02.cell_type_Trajectory.png", width = 6, height = 6, units = "in", dpi = 300)
ggsave("02.cell_type_Trajectory.pdf", width = 6, height = 6, units = "in", dpi = 300)

hubgene <- read.delim2('/data/nas1/luchunlin/project/JNZK-214-8/08_Lasso/lasso_genes.csv',header = F)

plot16 <- plot_genes_in_pseudotime(Mono.cds[hubgene$V1,],
                                   color_by = "Cell Type",
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





my_cds_subset=Mono.cds
# 拟时序数据和细胞位置在pData 中
head(pData(my_cds_subset))

# 这个differentialGeneTest会比较耗费时间，测试每个基因的拟时序表达
my_pseudotime_de <- differentialGeneTest(my_cds_subset,fullModelFormulaStr = "~sm.ns(Pseudotime)",cores = 4 )#cores调用的核心数

head(my_pseudotime_de)

# cluster the top 50 genes that vary as a function of pseudotime
my_pseudotime_de %>% arrange(qval) %>% head(50) %>% select(gene_short_name) -> gene_to_cluster
gene_to_cluster <- gene_to_cluster[,1]
gene_to_cluster
colnames(pData(my_cds_subset))
table(pData(my_cds_subset)$seurat_clusters,pData(my_cds_subset)$State) 
ac=pData(my_cds_subset)[c('Cell_Type','State','Pseudotime')]
head(ac)

# 这个热图绘制的并不是纯粹的细胞基因表达量矩阵，而是被 Pseudotime 好了的100列，50行的矩阵

my_pseudotime_cluster <- plot_pseudotime_heatmap(my_cds_subset[gene_to_cluster,],# num_clusters = 2, # add_annotation_col = ac,
                                                 show_rownames = TRUE,
                                                 return_heatmap = TRUE)

my_pseudotime_cluster
# 
# pdf('monocle_top50_heatmap.pdf')
# print(my_pseudotime_cluster)
# dev.off()
# my_branched_heatmap <- plot_genes_branched_heatmap(my_cds_subset[row.names(subset(BEAM_res, qval < 1e-4)),],branch_point = 1,num_clusters = 4, use_gene_short_name = TRUE,show_rownames = F,return_heatmap = TRUE)
# #将所做热图的基因和cluster提取出来
# head(my_branched_cluster$annotation_row)
# table(my_branched_heatmap$annotation_row$Cluster) 
# my_row <- my_branched_heatmap$annotation_row
# my_row <- data.frame(cluster = my_row$Cluster,
#                      gene = row.names(my_row),
#                      stringsAsFactors = FALSE)
# 
# head(my_row[my_row$cluster == 3,'gene']) 
# 
# my_gene <- row.names(subset(fData(my_cds_subset),
#                             gene_short_name %in% head(my_row[my_row$cluster == 2,'gene'])))
# my_gene

# cds$id <- rownames(cds)
# library(dplyr)
# p3 <- plot_pseudotime_heatmap(Mono.cds,
# #                        num_clusters = 3,
#                         cores = 1,
# # 
# show_rownames = T,return_heatmap = T)
# p3$tree_row
# cluster <- cutree(p3$tree_row,k=3)
# clustering <- data.frame(clusters)
# clustering[,1] <- as.character(clustering[,1])
# colnames(clustering) <- "Gene_Clusters"
# table(clustering)
## 基因名

# cluster <- p3$tree_row$labels

#UMAP$cell.type = Idents(UMAP)

setwd("/data/nas1/luchunlin/project/JNZK-214-8/04_scRNA/")
if (! dir.exists("./05_CellChat")){
  dir.create("./05_CellChat")
}
setwd("./05_CellChat")
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
#> Do heatmap based on a merged object
# netVisual_heatmap(cellchat, measure = "weight")
#> Do heatmap based on a merged object


#### circle ####
pdf('02.net_circle_number.pdf',w=4,h=5,family='serif')
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

png('02.net_circle_number.png',w=4,h=5,units='in',res=600,family='serif')
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

mat <- cellchat@net$count
pdf('03.single_circle.pdf',w=8,h=8)
par(mfrow=c(3,3),xpd=T)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2,vertex.weight = groupSize,weight.scale = T,arrow.width = 0.2,
                   arrow.size = 0.1,edge.weight.max = max(mat),title.name = rownames(mat)[i])
}
dev.off()
png('03.single_circle.png',w=700,h=700)
par(mfrow=c(3,3),xpd=T)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0,nrow = nrow(mat),ncol = ncol(mat),dimnames = dimnames(mat))
  mat2[i,] <- mat[i,]
  netVisual_circle(mat2,vertex.weight = groupSize,weight.scale = T,arrow.width = 0.2,
                   arrow.size = 0.1,edge.weight.max = max(mat),title.name = rownames(mat)[i])
}
dev.off()


setwd("/data/nas1/luchunlin/project/JNZK-214-8/04_scRNA/")
if (! dir.exists("./07_ssGSEA")){
  dir.create("./07_ssGSEA")
}
setwd("./07_ssGSEA")