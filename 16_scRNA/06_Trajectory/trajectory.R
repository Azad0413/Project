setwd("/data/nas1/luchunlin/project/BJTC-406-12/16_scRNA/")
if (! dir.exists("./06_Trajectory")){
  dir.create("./06_Trajectory")
}
setwd("./06_Trajectory")

library(monocle)
library(Seurat)
library(GEOquery)
library(data.table)
library(tidyverse)
library(lance)
load('/data/nas1/luchunlin/project/BJTC-406-12/scRNA.RData')

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