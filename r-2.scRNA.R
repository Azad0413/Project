rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-370-8/")
if (! dir.exists("./02_scRNA")){
  dir.create("./02_scRNA")
}
setwd("./02_scRNA")
library(Seurat)
library(GEOquery)
library(data.table)

## PART A导入数据--------
## 导入数据，创建seurat对象
## GSE166173
dir='GSE166173_RAW/GSE166173_RAW/'
list.files(dir)
GSM5065164<-Read10X(data.dir ='GSE166173_RAW/GSE166173_RAW/GSM5065164_8_RAW/')
GSM5065165<-Read10X(data.dir ='GSE166173_RAW/GSE166173_RAW/GSM5065165_9_RAW/')
GSM5065166<-Read10X(data.dir ='GSE166173_RAW/GSE166173_RAW/GSM5065166_10_RAW/')
GSM5065167<-Read10X(data.dir ='GSE166173_RAW/GSE166173_RAW/GSM5065167_11_RAW/')
count1<-CreateSeuratObject(counts = GSM5065164,min.cells = 3,min.features = 200)
count2<-CreateSeuratObject(counts = GSM5065165,min.cells = 3,min.features = 200)
count3<-CreateSeuratObject(counts = GSM5065166,min.cells = 3,min.features = 200)
count4<-CreateSeuratObject(counts = GSM5065167,min.cells = 3,min.features = 200)
scRNA<-merge(x=count1,y=c(count2,count3,count4),
             add.cell.ids=c('GSM5065164','GSM5065165','GSM5065166','GSM5065167'),
             merge.data=T)
##初步过滤，<=3个细胞中表达的基因，>=200个基因的细胞
scRNA
length(colnames(scRNA))#16988细胞数量
length(rownames(scRNA))# 24254基因数量
## ##PART B质控-------
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
                subset = nFeature_RNA > 200 & nFeature_RNA < 2500&
                  percent.mt < 5)

VlnPlot(scRNA,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
if (! dir.exists("./01_QC")){dir.create("./01_QC")}
setwd("./01_QC")
pdf(file = '01.QC.pdf',w=8,h=6)
VlnPlot(scRNA,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
dev.off()
png(file = '01.QC.png',w=800,h=600)
VlnPlot(scRNA,features = c('nFeature_RNA','nCount_RNA','percent.mt'),ncol = 3)
dev.off()
length(colnames(scRNA))#3194细胞数量
length(rownames(scRNA))# 24254基因数量
## PART C数据降维--------
setwd(('/data/nas1/luchunlin/project/BJTC-370-8/02_scRNA/'))
if (! dir.exists("./02_PCA")){dir.create("./02_PCA")}
setwd("./02_PCA")
#标准化
scRNA.norm<-NormalizeData(scRNA,normalization.method = "LogNormalize",scale.factor = 10000)
##scRNA.norm<-scRNA
## PART D PCA UMAP降维聚类-------
## 缩放数据
## 筛选1500个高变基因
all.genes<-rownames(scRNA.norm)
scRNA.norm<-FindVariableFeatures(scRNA.norm,selection.method = "vst", nfeatures = 3000)
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
          width = 6,
          height = 6,
          devices = NULL,
          res = 600,
          show = F)
write_fig(plot4,
          file = "01.feature_selection.png",
          width = 6,
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
library(Ipaper)
write_fig(plot5,
          file = "02.pca_cluster.pdf",
          width = 13,
          height = 8,
          devices = NULL,
          res = 600,
          show = F)
write_fig(plot5,
          file = "02.pca_cluster.png",
          width = 13,
          height = 8,
          devices = NULL,
          res = 300,
          show = F)
plot6 <- ElbowPlot(scRNA.norm.pca, ndims = 50)
plot6
write_fig(plot6,
          file = "03.pca_sd.pdf",
          width = 8,
          height = 6,
          devices = NULL,
          res = 600,
          show = F)
write_fig(plot6,
          file = "03.pca_sd.png",
          width = 8,
          height = 6,
          devices = NULL,
          res = 300,
          show = F)
scRNA.norm.pca.c<-FindNeighbors(scRNA.norm.pca,dims = 1:48)
scRNA.norm.pca.c<-FindClusters(scRNA.norm.pca.c,resolution = 0.4)
UMAP<-RunUMAP(scRNA.norm.pca.c,dims = 1:48)
DimPlot(UMAP,reduction = 'umap')

pdf('04.UMAP.pdf',w=8,h=7)
DimPlot(UMAP,reduction = 'umap')
dev.off()
png('04.UMAP.png',w=700,h=600)
DimPlot(UMAP,reduction = 'umap')
dev.off()
#-tsne-----
# set.seed(123)
# TSNE <- RunTSNE(scRNA.norm.pca.c, dims = 1:27)
# head(TSNE@reductions$tsne@cell.embeddings)
# p5 <- DimPlot(TSNE, reduction = "tsne", label = T)
# p5
# saveRDS(TSNE,file = 'TSNE.rds')
## 寻找差异基因------
# cluster1.markers<-FindMarkers(UMAP,ident.1 = 1,min.pct = 0.25)
# head(cluster1.markers,n=5)
all.markers<-FindAllMarkers(UMAP,only.pos = TRUE,min.pct = 0.25,logfc.threshold = 1,test.use = 'wilcox',return.thresh = 0.05)
all.markers%>%group_by(cluster)%>%top_n(n=2,wt=avg_log2FC)
write.table(all.markers,file = 'DEGs.xls',sep = '\t',row.names = F,quote = F)
top10<-all.markers%>%group_by(cluster)%>%top_n(n=10,wt=avg_log2FC)

# exp<-as.matrix(UMAP@assays$RNA@counts)
# smps = Cells(UMAP) %>% str_sub(1,10)
# table(smps)
# grps = data.frame(cell=c('GSM5065164','GSM5065165','GSM5065166','GSM5065167'),group=c('control','RB','control','RB'))
# group.label = grps$group[match(smps,grps$cell)]
# names(group.label) = Cells(UMAP)
# group.label
# UMAP = AddMetaData(UMAP,group.label,col.name = 'group ')
# group<-as.data.frame(group.label)
# dat<-exp
# dat <- na.omit(dat)
# library(DESeq2)
# colData<-group
# colData$group<-factor(colData$group,levels = c('control','RB'))
# dds<-DESeqDataSetFromMatrix(countData = dat,colData=colData,design = ~group)
# dds = dds[rownames(counts(dds)) > 1,]
# dds<-estimateSizeFactors(dds)
# ##提取标准化后的数据 
# #normalized_counts <- counts(dds,normalized=T) 
# # write.table(normalized_counts,file = 'normalized.counts.xls',sep = '\t',row.names = T,quote = F)
# dds<-DESeq(dds)
# ## 提取差异结果
# res =results(dds, contrast = c("group","Tumor","Normal"))
# res =res[order(res$padj),]
# head(res)
# summary(res)
# table(res$padj<0.05)
# DEG <- subset(res, padj < 0.05 & abs(log2FoldChange) >1 )
# DEG<-as.data.frame(res)
# DEG<-na.omit(DEG)
# dim(DEG)
# head(DEG)
# ## 添加change列
# logFC_cutoff<-1
# DEG$change=as.factor(
#   ifelse(DEG$padj<0.05&abs(DEG$log2FoldChange)>logFC_cutoff,
#          ifelse(DEG$log2FoldChange>logFC_cutoff,'UP','DOWN'),'NOT'))
# table(DEG$change)
# ## DOWN   NOT    UP 
# ## 2003 14547  3041 
# sig_diff <- subset(DEG,
#                    DEG$padj < 0.05 & abs(DEG$log2FoldChange) >= logFC_cutoff)
# ## 5044
# DEG_write <- cbind(GeneSymbol=rownames(DEG), DEG)
# write.table(DEG_write, file = "DEG_all(mRNA).xls",
#             quote = F,
#             sep = "\t",
#             row.names = F)
# sig_diff_write <- cbind(GeneSymbol=rownames(sig_diff), sig_diff)
# write.table(sig_diff_write, file = "DEG_sig(mRNA).xls",
#             quote = F,
#             sep = "\t",
#             row.names = F)
# ### 火山图---------
# #devtools::install_github("kongdd/Ipaper")
# library(ggplot2)
# library(ggthemes)
# library(RColorBrewer)
# library(Ipaper)
# library(scales)
# library(ggrepel)
# dat_rep<-DEG[rownames(DEG)%in%
#                rownames(rbind(head(sig_diff[order(sig_diff$log2FoldChange,decreasing = T),],10),
#                               head(sig_diff[order(sig_diff$log2FoldChange,decreasing = F),],10))),]
# volcano_plot<- ggplot(data = DEG, 
#                       aes(x = log2FoldChange,
#                           y = -log10(padj), 
#                           color =change)) +
#   scale_color_manual(values = c("blue", "darkgray","red")) +
#   scale_x_continuous(breaks = c(-10,-6,-3,-1,0,1,3,6,10)) +
#   scale_y_continuous(trans = "log1p",
#                      breaks = c(0,1,5,10,20,50, 100,200)) +
#   geom_point(size = 1.2, alpha = 0.4, na.rm=T) +
#   theme_bw(base_size = 12, base_family = "Times") +
#   geom_vline(xintercept = c(-1,1),
#              lty = 4,
#              col = "darkgray",
#              lwd = 0.6)+
#   geom_hline(yintercept = -log10(0.05),
#              lty = 4,
#              col = "darkgray",
#              lwd = 0.6)+
#   theme(legend.position = "right",
#         panel.grid = element_blank(),
#         legend.title = element_blank(),
#         legend.text = element_text(face="bold",
#                                    color="black",
#                                    family = "Times",
#                                    size=13),
#         plot.title = element_text(hjust = 0.5),
#         axis.text.x = element_text(face = "bold",
#                                    color = "black",
#                                    size = 15),
#         axis.text.y = element_text(face = "bold",
#                                    color = "black",
#                                    size = 15),
#         axis.title.x = element_text(face = "bold",
#                                     color = "black",
#                                     size = 15),
#         axis.title.y = element_text(face = "bold",
#                                     color = "black",
#                                     size = 15)) +
#   geom_label_repel(
#     data = dat_rep,
#     aes(label = rownames(dat_rep)),
#     max.overlaps = 20,
#     size = 4,
#     box.padding = unit(0.5, "lines"),
#     min.segment.length = 0,
#     point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )+
#   labs(x = "log (Fold Change)",
#        y = "-log10 (adj.P.Val)")
# volcano_plot
# ggsave('volcano(mRNA).png', volcano_plot,width = 8, height = 6)
# ggsave('volcano(mRNA).pdf', volcano_plot,width = 8, height = 6)
setwd("/data/nas1/luchunlin/project/BJTC-370-8/02_scRNA/")
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
                     quantile = 0.5,
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
# ferr <- read.table('/data/nas2/database/gene_set/Ferr_genes.txt',header = T)
# 
# inter <- all.markers[all.markers$gene%in%ferr$Ferr_gene,]
## 112