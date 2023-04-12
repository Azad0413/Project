rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-406-12/")
if (! dir.exists("./17_pySCENIC")){
  dir.create("./17_pySCENIC")
}
setwd("./17_pySCENIC")
#devtools::install_github("aertslab/SCENIC")
library(SCENIC)

## SCENIC需要一些依赖包，先安装好
# BiocManager::install("AUCell")
#BiocManager::install(c("GENIE3"),ask = F,update = F)
# BiocManager::install("RcisTarget",force = TRUE)

library(AUCell)
library(GENIE3)
library(RcisTarget)


##可视化
library(Seurat)
library(SCopeLoomR)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
#BiocManager::install('scRNAseq')
library(scRNAseq)

## 提取 out_SCENIC.loom 信息
#inputDir='./outputs/'
#scenicLoomPath=file.path(inputDir,'out_SCENIC.loom')
library(SCENIC)
loom <- open_loom('out_SCENIC.loom') 

regulons_incidMat <- get_regulons(loom, column.attr.name = 'Regulons')
regulons_incidMat[1:4,1:4] 
regulons <- regulonsToGeneLists(regulons_incidMat)

regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds <- get_regulon_thresholds(loom)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])

embeddings <- get_embeddings(loom)  
close_loom(loom)

rownames(regulonAUC)
names(regulons)


sce<-readRDS(file = '../16_scRNA/07_consensus/UMAP.cluster.rds')
table(sce$seurat_clusters)
sce$singleR

library(ggplot2) 
genes_to_check = read.delim2(file = '../08_Lasso/lasso_genes.csv',header = F)
library(stringr)  
genes_to_check=str_to_upper(genes_to_check$V1)
genes_to_check
colnames(sce)
###2个亚型分开
group <- read.delim2('../16_scRNA/07_consensus/cluster.xls')
table(group$cluster)
cluster1.sample <- group$sample[which(group$cluster=='cluster 1')]
cluster2.sample <- group$sample[which(group$cluster=='cluster 2')]
sce1 <- sce[,colnames(sce)%in%cluster1.sample]
sce2 <- sce[,colnames(sce)%in%cluster2.sample]
sub_regulonAUC1 <- regulonAUC[,match(cluster1.sample,colnames(regulonAUC))]
dim(sub_regulonAUC1)
sub_regulonAUC2 <- regulonAUC[,match(cluster2.sample,colnames(regulonAUC))]
dim(sub_regulonAUC2)
#确认是否一致
identical(colnames(sub_regulonAUC1), colnames(sce1))
identical(colnames(sub_regulonAUC2), colnames(sce2))
#[1] TRUE 

## Cluster1------------
cellClusters1 <- data.frame(row.names = colnames(sce1), 
                           seurat_clusters = as.character(sce1$seurat_clusters))
cellTypes1 <- data.frame(row.names = colnames(sce1), 
                        celltype = sce1$singleR)
cellTypes1$celltype <- 'Macrophage'
head(cellTypes1)
head(cellClusters1)
sub_regulonAUC1[1:4,1:4] 
save(sub_regulonAUC1,cellTypes1,cellClusters1,sce1,
     file = 'for_rss_and_visual1.Rdata')

dim(sub_regulonAUC1)

# Split the cells by cluster:
selectedResolution <- "celltype" # select resolution
cellsPerGroup1 <- split(rownames(cellTypes1), 
                       cellTypes1[,selectedResolution]) 

# 去除extened regulons
sub_regulonAUC1 <- sub_regulonAUC1[onlyNonDuplicatedExtended(rownames(sub_regulonAUC1)),] 
dim(sub_regulonAUC1)
# Calculate average expression:
regulonActivity_byGroup1 <- sapply(cellsPerGroup1,
                                  function(cells) 
                                    rowMeans(getAUC(sub_regulonAUC1)[,cells]))

# Scale expression. 
# Scale函数是对列进行归一化，所以要把regulonActivity_byGroup转置成细胞为行，基因为列
# 参考：https://www.jianshu.com/p/115d07af3029
# regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
#                                           center = T, scale=T)) 
regulonActivity_byGroup_Scaled1 <- scale(regulonActivity_byGroup1)
# 同一个regulon在不同cluster的scale处理
dim(regulonActivity_byGroup_Scaled1)
#[1] 194   1
regulonActivity_byGroup_Scaled1=regulonActivity_byGroup_Scaled1[]
regulonActivity_byGroup_Scaled1=na.omit(regulonActivity_byGroup_Scaled1)
hubAUC<-regulons_incidMat[,genes_to_check]
hubAUC<-hubAUC[rowSums(hubAUC)>0,]
regulonActivity_byGroup_Scaled1<-regulonActivity_byGroup_Scaled1[rownames(hubAUC),]
# pdf(file = '01.heatmap.pdf',w=7,h=9)
regulonActivity_byGroup_Scaled1 <- as.data.frame(regulonActivity_byGroup_Scaled1)
colnames(regulonActivity_byGroup_Scaled1) <- 'rank'
pheatmap(regulonActivity_byGroup_Scaled1)
# dev.off()
# png(file = '01.heatmap.png',w=600,h=800)
# pheatmap(regulonActivity_byGroup_Scaled)
# dev.off()

library(ggrepel)
### RSS  regulon specificity score
rss1 <- calcRSS(AUC = getAUC(sub_regulonAUC1),
                cellAnnotation = cellTypes1[colnames(sub_regulonAUC1),selectedResolution])
rss1=na.omit(rss1)%>%as.data.frame() %>%rownames_to_column(var = 'TF')
colnames(rss1) <- c('TF','RSS')

rss1 = rss1 %>%
  arrange(desc(RSS)) %>%
  mutate("rank" = row_number())
write.table(rss1,file = 'cluster1.rss.xls',sep = '\t',row.names = F,quote = F)
rep.dat <- head(rss1[order(rss1$RSS,decreasing = T),],5)

p1 <- ggplot(rss1,
       aes(
         x = rank,
         y = RSS
       )
)+
  geom_point(color="skyblue")+
  theme_bw()+
  labs(
    x = "Rank",
    y = "Regulon Specificity Score",
    title = 'Cluster 1'
    )+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,1), # 调整y轴的坐标轴样式为10^1,10^2...
                     labels = c(0,0.2,0.4,0.6,0.8,1))+
  geom_label_repel(data = rep.dat,
                   aes(label=TF),
                   max.overlaps = 20)+
  theme(axis.text.x = element_blank())
p1
ggsave(filename = '01.RSS(cluster1).pdf',w=7,h=5)
ggsave(filename = '01.RSS(cluster1).png',w=6,h=4)
top10_1 <- head(rss1[order(rss1$RSS,decreasing = T),],10)
#write.table(top10,file = 'cluster1(top10).xls',sep = '\t',row.names = F,quote = F)

# rssPlot1 <- plotRSS(rss1)
# plotly::ggplotly(rssPlot1$plot)

### UMAP表达图
TF.gene1 <- gsub('(+)','',top10_1$TF,fixed = T)
pdf(file = '02.umap(cluster1).pdf',w=10,h=7)
FeaturePlot(sce,features = TF.gene1)
dev.off()

png(file = '02.umap(cluster1).png',w=750,h=500)
FeaturePlot(sce,features = TF.gene1)
dev.off()

## Cluster2------------
cellClusters2 <- data.frame(row.names = colnames(sce2), 
                            seurat_clusters = as.character(sce2$seurat_clusters))
cellTypes2 <- data.frame(row.names = colnames(sce2), 
                         celltype = sce2$singleR)
cellTypes2$celltype <- 'Macrophage'
head(cellTypes2)
head(cellClusters2)
sub_regulonAUC2[2:4,2:4] 
save(sub_regulonAUC2,cellTypes2,cellClusters2,sce2,
     file = 'for_rss_and_visual2.Rdata')

dim(sub_regulonAUC2)

# Split the cells by cluster:
selectedResolution <- "celltype" # select resolution
cellsPerGroup2 <- split(rownames(cellTypes2), 
                        cellTypes2[,selectedResolution]) 

# 去除extened regulons
sub_regulonAUC2 <- sub_regulonAUC2[onlyNonDuplicatedExtended(rownames(sub_regulonAUC2)),] 
dim(sub_regulonAUC2)
# Calculate average expression:
regulonActivity_byGroup2 <- sapply(cellsPerGroup2,
                                   function(cells) 
                                     rowMeans(getAUC(sub_regulonAUC2)[,cells]))

# Scale expression. 
# Scale函数是对列进行归一化，所以要把regulonActivity_byGroup转置成细胞为行，基因为列
# 参考：https://www.jianshu.com/p/225d07af3029
# regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),
#                                           center = T, scale=T)) 
regulonActivity_byGroup_Scaled2 <- scale(regulonActivity_byGroup2)
# 同一个regulon在不同cluster的scale处理
dim(regulonActivity_byGroup_Scaled2)
#[2] 194   1
regulonActivity_byGroup_Scaled2=regulonActivity_byGroup_Scaled2[]
regulonActivity_byGroup_Scaled2=na.omit(regulonActivity_byGroup_Scaled2)
hubAUC<-regulons_incidMat[,genes_to_check]
hubAUC<-hubAUC[rowSums(hubAUC)>0,]
regulonActivity_byGroup_Scaled2<-regulonActivity_byGroup_Scaled2[rownames(hubAUC),]
# pdf(file = '02.heatmap.pdf',w=7,h=9)
regulonActivity_byGroup_Scaled2 <- as.data.frame(regulonActivity_byGroup_Scaled2)
colnames(regulonActivity_byGroup_Scaled2) <- 'rank'
pheatmap(regulonActivity_byGroup_Scaled2)
# dev.off()
# png(file = '02.heatmap.png',w=600,h=800)
# pheatmap(regulonActivity_byGroup_Scaled)
# dev.off()

library(ggrepel)
### RSS  regulon specificity score
rss2 <- calcRSS(AUC = getAUC(sub_regulonAUC2),
                cellAnnotation = cellTypes2[colnames(sub_regulonAUC2),selectedResolution])
rss2=na.omit(rss2)%>%as.data.frame() %>%rownames_to_column(var = 'TF')
colnames(rss2) <- c('TF','RSS')

rss2 = rss2 %>%
  arrange(desc(RSS)) %>%
  mutate("rank" = row_number())
write.table(rss2,file = 'cluster2.rss.xls',sep = '\t',row.names = F,quote = F)
rep.dat <- head(rss2[order(rss2$RSS,decreasing = T),],5)

p2 <- ggplot(rss2,
             aes(
               x = rank,
               y = RSS
             )
)+
  geom_point(color="skyblue")+
  theme_bw()+
  labs(
    x = "Rank",
    y = "Regulon Specificity Score",
    title = 'Cluster 2'
  )+
  scale_y_continuous(breaks = c(0,0.2,0.4,0.6,0.8,2), # 调整y轴的坐标轴样式为20^2,20^2...
                     labels = c(0,0.2,0.4,0.6,0.8,2))+
  geom_label_repel(data = rep.dat,
                   aes(label=TF),
                   max.overlaps = 20)+
  theme(axis.text.x = element_blank())
p2
ggsave(filename = '03.RSS(cluster2).pdf',w=7,h=5)
ggsave(filename = '03.RSS(cluster2).png',w=6,h=4)
top10_2 <- head(rss2[order(rss2$RSS,decreasing = T),],10)
#write.table(top20,file = 'cluster2(top20).xls',sep = '\t',row.names = F,quote = F)

# rssPlot2 <- plotRSS(rss2)
# plotly::ggplotly(rssPlot2$plot)

### UMAP表达图
TF.gene2 <- gsub('(+)','',top10_2$TF,fixed = T)
pdf(file = '04.umap(cluster2).pdf',w=10,h=7)
FeaturePlot(sce,features = TF.gene2)
dev.off()

png(file = '04.umap(cluster2).png',w=750,h=500)
FeaturePlot(sce,features = TF.gene2)
dev.off()


### 没有在2组中同时出现的TF
diff.tf <- rbind(TF.gene1[!TF.gene1%in%TF.gene2],TF.gene2[!TF.gene2%in%TF.gene1])

library(lance)
dat<-sce@assays$RNA@counts%>%as.data.frame()
dat <- edgeR::cpm(dat)
hub_exp<-dat[diff.tf,]
hub_exp2<-log2(hub_exp+1)%>%as.data.frame()
hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))

## 样本分组
hub_exp2$Group<-ifelse(hub_exp2$sample%in%cluster1.sample,'cluster 1','cluster 2')
hub_exp2 <- hub_exp2[which(hub_exp2$expr>0),]
library(rstatix)
stat.test<-hub_exp2%>%
  group_by(Symbol)%>%
  wilcox_test(expr ~ Group)%>%
  adjust_pvalue(method = 'fdr')

write.table(stat.test,file = 'wilcox.validation.xls',sep = '\t',row.names = F,quote = F)

##分面图形
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
exp_plot <- ggplot(hub_exp2,aes(x = Group, y = expr, fill = Group)) +
  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.4,
               position=position_dodge(0.9),
               outlier.shape = NA,fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#4682B4","#CD3700"), name = "Group")+
  ylim(c(4,12))+
  labs(title="", x="", y = "",size=20) +
  stat_compare_means(data = hub_exp2,
                     mapping = aes(group = Group),
                     label ="p.signif",
                     method = 'wilcox.test',label.x = 1.45) +
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=15),
        axis.text.x=element_text(angle=0,hjust=0.5,colour="black",face="bold",size=12), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+facet_wrap(~Symbol,scales = "free",nrow = 1) 
exp_plot
ggsave(filename = '05.tfexp.pdf',exp_plot,w=6,h=4)
ggsave(filename = '05.tfexp.png',exp_plot,w=6,h=4)
