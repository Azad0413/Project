## rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/BJTC-292")
if (! dir.exists("./00_rawdata")){
  dir.create("./00_rawdata")
}
setwd("./00_rawdata")
##01-1 TCGA LAML---------
#BiocManager::install("EDASeq")
library(TCGAbiolinks)
library(readr)
library(dplyr)
library(EDASeq)

## 读取从xena下载的数据
expr<-read_tsv(file = 'TCGA-LAML.htseq_counts.tsv')
expr<-as.data.frame(expr)
rownames(expr)<-expr[,1]
expr<-expr[,-1]
## xena下载的数据经过了log2+1转化，需要将其还原
expr<-2^expr-1
## 对数据进行id转化
genecode<-read.table(file = 'gencode.v22.annotation.gene.probeMap')
probe2symbol<-genecode[,(1:2)]
colnames(probe2symbol)<-c('ID','symbol')
probe2symbol<-probe2symbol[-1,]
dat<-expr
dat$ID <- rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symbol$ID<-as.character(probe2symbol$ID)
dat<-dat %>%
  inner_join(probe2symbol,by='ID')%>% 
  select(-ID)%>%     ## 去除多余信息
  select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('TCGA',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
dim(dat)

write.table(dat,
            file = "dataPrep.xls",
            quote = F)

#将数据进行标准化
dataNorm.AML <- TCGAanalyze_Normalization(tabDF = dat,
                                           geneInfo = geneInfo,
                                           method = "gcContent")
# 将标准化后的数据再过滤，得到最终的数据
dataFilt.AML.final <- TCGAanalyze_Filtering(tabDF = dataNorm.AML,
                                             method = "quantile", 
                                             qnt.cut =  0.25)
dim(dataFilt.AML.final)
##  13627   151

write.table(dataFilt.AML.final,
            file = "dataFilt.AML.final.xls",
            quote = F)

## 01-2 TARGET-ALM-------
expr_va1<-read_tsv(file = 'TARGET-AML.htseq_counts.tsv')
expr_va1<-as.data.frame(expr_va1)
rownames(expr_va1)<-expr_va1[,1]
expr_va1<-expr_va1[,-1]
dat_va1<-expr_va1
dat_va1$ID <- rownames(dat_va1)
dat_va1$ID<-as.character(dat_va1$ID)
probe2symbol$ID<-as.character(probe2symbol$ID)
dat_va1<-dat_va1 %>%
  inner_join(probe2symbol,by='ID')%>% 
  select(-ID)%>%     ## 去除多余信息
  select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('TCGA',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除

dim(dat_va1)

## 01-3 GSE71014-------
library(GEOquery)
library(Biobase)
## GSE71014
gset_va<-getGEO("GSE71014",
                destdir = '.',
                GSEMatrix = T,
                getGPL = F)
expr_va2<-as.data.frame(exprs(gset_va[[1]]))
a_va=gset_va[[1]]
pd_va<-pData(a_va)
library(AnnoProbe)
gpl<-"GPL10558"
probe2symbol2=idmap(gpl,type = 'pipe')
colnames(probe2symbol2)<-c('ID','symbol')
dat_va2<-expr_va2
dat_va2$ID<-rownames(dat_va2)
dat_va2$ID<-as.character(dat_va2$ID)
probe2symbol2$ID<-as.character(probe2symbol2$ID)
dat_va2<-dat_va2 %>%
  inner_join(probe2symbol2,by='ID')%>% 
  select(-ID)%>%     ## 去除多余信息
  select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除


# 02 Immune cluster------
setwd("/data/nas1/luchunlin/project/BJTC-292")
if (! dir.exists("./01_Immune_cluster")){
  dir.create("./01_Immune_cluster")
}
setwd("./01_Immune_cluster")

## 02-1 ssGSEA估算AML中免疫浸润水平------------

library(GSVA)
gene_set <- read.table("/data/nas1/luchunlin/project/BJTC-292/01_Immune_cluster/mmc3.txt",
                       header = T,
                       sep ="\t")
dat_gsva <- dataFilt.AML.final
gene_list <- split(as.matrix(gene_set)[,1],
                   gene_set[,2])

ssgsea_score = gsva(dat_gsva, gene_list, 
                    method = "ssgsea", 
                    ssgsea.norm = TRUE, 
                    verbose = TRUE)
write.table(ssgsea_score,
            file = "ssgsea_result.xls",
            sep = "\t",
            quote = F)

## 02-2 无监督聚类-------
dat_cluster<-t(ssgsea_score)%>%as.data.frame()
dat_cluster$cluster<-NULL
kc<-kmeans(dat_cluster,2)
kc$size
##  85 66  2个cluster
cluster<-data.frame(cluster=kc$cluster)
#install.packages('factoextra')
#library(factoextra)
fviz_cluster(object=kc,data=dat_cluster,
             ellipse.type = "convex",star.plot=T,repel=T,
             geom = c("point"),palette=c( "#FF6347","#6495ED"),main="",
             ggtheme=theme_bw())

 ## 02-3 ESTIMATE算法验证--------
group<-cluster
group$label<-ifelse(group$cluster==1,'Cluster1','Cluster2')
group_estimate<-group$label%>%as.factor()
design<-model.matrix(~0 + group_estimate)
rownames(design)<-rownames(group)
colnames(design)<-levels(group_estimate)
design<-as.data.frame(design)
Cluster1<-rownames(design)[which(design$Cluster1==1)]
Cluster2<-rownames(design)[which(design$Cluster2==1)]
length(Cluster1)
## 85
length(Cluster2)
## 66
#install.packages("estimate", repos="http://R-Forge.R-project.org")
library(estimate)
dat_estimate<-log2(dat_gsva+1)
write.table(dat_estimate, 
            'dat_estimate_log2.txt', 
            col.names = T, 
            row.names = T, 
            quote = F, sep="\t")
filterCommonGenes(input.f = './dat_estimate_log2.txt', 
                  output.f = 'dat_estimate.gct', 
                  id = 'GeneSymbol')
estimateScore('dat_estimate.gct', 'dat_purity.gct', platform="affymetrix")
# [1] "1 gene set: StromalSignature  overlap= 112"
# [1] "2 gene set: ImmuneSignature  overlap= 140"
es_score <- read.table('dat_purity.gct', skip = 2, header = T)
immu_score <- es_score[,3:length(es_score)]
rownames(immu_score) <- es_score$NAME
write.table(es_score,
            file = "es_score.xls",
            sep = "\t",
            quote = F,
            row.names = F)

## 热图---------
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(gplots)
rt_group <- data.frame(t(immu_score))
rownames(rt_group)<-rownames(group)
rt_group$group <- ifelse(rownames(rt_group) %in% Cluster1,
                           "Cluster1", "Cluster2")
rt_group<-rt_group[order(rt_group$group),]
rt_group$ImmuneScore<-cut(rt_group$ImmuneScore,breaks = c(600,1200,1800,2400,3000,3600),labels = c('900','1500','2100','2700','3300'))
rt_group$StromalScore<-cut(rt_group$StromalScore,breaks = c(-2400,-1800,-1200,-600,0,600),labels = c('-2100','-1500','-900','-300','300'))
rt_group$ESTIMATEScore<-cut(rt_group$ESTIMATEScore,breaks = c(-1400,-400,600,1600,2600,3600),labels = c('-900','100','1100','2100','3100'))
rt_group$TumorPurity<-cut(rt_group$TumorPurity,breaks = c(0.35,0.45,0.55,0.65,0.75,0.85,0.95),labels = c('0.4','0.5','0.6','0.7','0.8','0.9'))
rt_group$StromalScore<-as.factor(rt_group$StromalScore)
rt_group$ImmuneScore<-as.factor(rt_group$ImmuneScore)
rt_group$ESTIMATEScore<-as.factor(rt_group$ESTIMATEScore)
rt_group$TumorPurity<-as.factor(rt_group$TumorPurity)
rt_group$group<-as.factor(rt_group$group)
rt_group<-rt_group[,c(5,1:4)]
rt_dat<-ssgsea_score[,rownames(rt_group)]
annotation_col<-rt_group
ann_colors<-list(
  group = c(Cluster1="#87CEFA",Cluster2="#FFB6C1"),
  ImmuneScore=c('900'='#B0E0E6','1500'='#ADD8E6','2100'='#87CEFA','2700'='#00BFFF','3300'='#6495ED'),
  StromalScore=c('-2100'='#98F898','-1500'='#90EE90','-900'='#00DDAA','-300'='#66CDAA','300'='#20B2AA'),
  ESTIMATEScore=c('-900'='#E6E6FA','100'='#D8BFD8','1100'='#DDA0DD','2100'='#EE82EE','3100'='#BA55D3'),
  TumorPurity=c('0.4'='#EEFFBB','0.5'='#CCFF99','0.6'='#CCFF33','0.7'='#BBFF66','0.8'='#99FF99','0.9'='#66FF66'))

pheatmap(rt_dat,
  color = bluered(100),
  border_color = NA,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  labels_row = NULL,
  clustering_method = 'ward.D2',
  show_rownames = T,
  show_colnames = F,
  fontsize_col = 5,
  cluster_cols = F,
  cluster_rows = T)

## 小提琴图--------
violin_dat <- data.frame(t(immu_score))
rownames(violin_dat)<-rownames(group)
violin_dat$sample <- rownames(violin_dat)
violin_dat$group <- ifelse(violin_dat$sample %in% Cluster1,
                           "Cluster1", "Cluster2")
library(tidyr)
library(ggplot2)
library(ggpubr)
library(Ipaper)
{
  p1 <- ggplot(violin_dat, aes(x=group,y=StromalScore, fill=group))+
    geom_violin() + #绘制小提琴图
    stat_boxplot(geom="errorbar",
                 width=0.1,
                 position = position_dodge(0.9)) +
    geom_boxplot(width=0.4,
                 position=position_dodge(0.9),
                 outlier.shape = NA,
                 fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
    geom_point(aes(fill = group),
               size = 0.5,
               position = position_dodge(0.9))+
    scale_fill_manual(values = c("#FF6347","#6495ED"), name = "Group") + 
    stat_compare_means(aes(group = group),method = 'wilcox')+ 
    theme_bw()+ #背景变为白色
    ylim(-2500,1500) +
    theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
          panel.grid = element_blank()) + 
    theme(axis.text.x=element_text(family="Times",size=12,face="bold"), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
          axis.text.y=element_text(family="Times",size=12,face="bold"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
          axis.title.y=element_text(family="Times",size = 15,face="bold"), #设置y轴标题的字体属性
          axis.title.x=element_text(family="Times",size = 10,face="bold"),
          plot.title = element_text(hjust = 0.5,family="Times",size = 15,face="bold"),
          legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold"))+
    labs(title="Stromal Score", x="", y="Stromal Score")
  p1
  p2 <- ggplot(violin_dat, aes(x=group,y=ImmuneScore, fill=group))+
    geom_violin() + #绘制小提琴图
    stat_boxplot(geom="errorbar",
                 width=0.1,
                 position = position_dodge(0.9)) +
    geom_boxplot(width=0.4,
                 position=position_dodge(0.9),
                 outlier.shape = NA,
                 fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
    geom_point(aes(fill = group),
               size = 0.5,
               position = position_dodge(0.9))+
    scale_fill_manual(values = c("#FF6347","#6495ED"), name = "Group") +
    stat_compare_means(aes(group = group),method = 'wilcox')+ 
    theme_bw()+ #背景变为白色
    ylim(0,4000) +
    theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
          panel.grid = element_blank()) + 
    theme(axis.text.x=element_text(family="Times",size=12,face="bold"), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
          axis.text.y=element_text(family="Times",size=12,face="bold"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
          axis.title.y=element_text(family="Times",size = 15,face="bold"), #设置y轴标题的字体属性
          axis.title.x=element_text(family="Times",size = 10,face="bold"),
          plot.title = element_text(hjust = 0.5,family="Times",size = 15,face="bold"),     
          legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold"))+
    labs(title="Immune Score", x="", y="Immune Score")
  p2
  p3 <- ggplot(violin_dat, aes(x=group, y=ESTIMATEScore, fill=group))+
    geom_violin() + #绘制小提琴图
    stat_boxplot(geom="errorbar",
                 width=0.1,
                 position = position_dodge(0.9)) +
    geom_boxplot(width=0.4,
                 position=position_dodge(0.9),
                 outlier.shape = NA,
                 fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
    geom_point(aes(fill = group),
               size = 0.5,
               position = position_dodge(0.9))+
    scale_fill_manual(values = c("#FF6347","#6495ED"), name = "Group") +
    stat_compare_means(aes(group = group),method = 'wilcox')+ 
    theme_bw()+ #背景变为白色
    theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
          panel.grid = element_blank()) + 
    ylim(-2000, 4500) +
    theme(axis.text.x=element_text(family="Times",size=12,face="bold"), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
          axis.text.y=element_text(family="Times",size=12,face="bold"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
          axis.title.y=element_text(family="Times",size = 15,face="bold"), #设置y轴标题的字体属性
          axis.title.x=element_text(family="Times",size = 10,face="bold"),
          plot.title = element_text(hjust = 0.5,family="Times",size = 15,face="bold"),
          legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold"))+
    labs(title="ESTIMATE Score", x="", y="ESTIMATE Score")
  p3
  p4 <- ggplot(violin_dat, aes(x=group, y=TumorPurity, fill=group))+
    geom_violin() + #绘制小提琴图
    stat_boxplot(geom="errorbar",
                 width=0.1,
                 position = position_dodge(0.9)) +
    geom_boxplot(width=0.4,
                 position=position_dodge(0.9),
                 outlier.shape = NA,
                 fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
    geom_point(aes(fill = group),
               size = 0.5,
               position = position_dodge(0.9))+
    scale_fill_manual(values = c("#FF6347","#6495ED"), name = "Group") +
    stat_compare_means(aes(group = group),method = 'wilcox')+ 
    theme_bw()+ #背景变为白色
    theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
          panel.grid = element_blank()) + 
    ylim(0.3, 1.1) +
    theme(axis.text.x=element_text(family="Times",size=12,face="bold"), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
          axis.text.y=element_text(family="Times",size=12,face="bold"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
          axis.title.y=element_text(family="Times",size = 15,face="bold"), #设置y轴标题的字体属性
          axis.title.x=element_text(family="Times",size = 10,face="bold"),
          plot.title = element_text(hjust = 0.5,family="Times",size = 15,face="bold"),
          legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold"))+
    labs(title="Tumor Purity", x="", y="Tumor Purity")
  p4
  p5 <- cowplot::plot_grid(p1,p2,p3,p4,
                           nrow = 2, 
                           align = 'h', 
                           vjust = -0.3)
  p5
}

## 02-4 GSEA富集分析--------
group_gsea<-data.frame(sample=rownames(rt_group),
                   group=c(rep('Cluster 1',85),rep('Cluster 2',66)))
gsea_exp<-dat_gsva[,rownames(rt_group)]
write.table(gsea_exp,
            file = "gsea_exp.xls",
            quote = F,
            sep = "\t",
            row.names = T)
#BiocManager::install('DESeq2')
library(DESeq2)
colData<-data.frame(sample=rownames(rt_group),
                    group=c(rep('Cluster1',85),rep('Cluster2',66)))

colData$group<-factor(colData$group,levels = c('Cluster1','Cluster2'))
write.table(colData,
            file = "group.xls",
            quote = F,
            row.names = F)
dds<-DESeqDataSetFromMatrix(countData = gsea_exp,colData=colData,design = ~group)
dds = dds[rownames(counts(dds)) > 1,]
dds<-DESeq(dds)
dds
res =results(dds, contrast = c("group","Cluster1","Cluster2"))
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
DEGeneSets <- DEGeneSets[order(DEGeneSets$padj),]
dim(DEGeneSets)
## GSEA KEGG-----
library(clusterProfiler)
library(enrichplot)
kegg_set<- read.gmt("c2.cp.kegg.v7.5.1.symbols.gmt")
set.seed(1)
kegg_gsea <- GSEA(geneList, TERM2GENE = kegg_set, pvalueCutoff = 0.05)
kegg_result <- kegg_gsea@result
dim(kegg_result)

gseaplot2(kegg_gsea,c(2,4,5,7,10),color = c('#7B68EE','#CD3333','#20B2AA','#FF8C00','#FF6666'),
          title = 'KEGG GSEA',
          base_size = 11,
          rel_heights = c(1.5, 0.3, 0.5))
## GSEA GO------
go_set<- read.gmt("c5.go.bp.v7.4.symbols.gmt")
set.seed(1)
go_gsea <- GSEA(geneList, TERM2GENE = go_set, eps = 0)
go_result <- go_gsea@result
dim(go_result)

gseaplot2(go_gsea,c(1:5),color = c('#7B68EE','#CD3333','#20B2AA','#FF8C00','#FF6666'),
          title = 'GO GSEA',
          base_size = 11,
          rel_heights = c(1.5, 0.3, 0.5))

# 03 DEGs------
setwd("/data/nas1/luchunlin/project/BJTC-292")
if (! dir.exists("./02_DEGs")){
  dir.create("./02_DEGs")
}
setwd("./02_DEGs")

DEG <- subset(res, padj < 0.05 & abs(log2FoldChange) >1 )
DEG<-as.data.frame(res)
DEG<-na.omit(DEG)
dim(DEG)
head(DEG)
## 添加change列
logFC_cutoff<-1
DEG$change=as.factor(
  ifelse(DEG$padj<0.05&abs(DEG$log2FoldChange)>logFC_cutoff,
         ifelse(DEG$log2FoldChange>logFC_cutoff,'UP','DOWN'),'NOT'))
table(DEG$change)
## DOWN   NOT    UP 
## 769 12352   506 
sig_diff <- subset(DEG,
                   DEG$padj < 0.05 & abs(DEG$log2FoldChange) >= logFC_cutoff)
## 1275
DEG_write <- cbind(GeneSymbol=rownames(DEG), DEG)
write.table(DEG_write, file = "DEG_all.xls",
            quote = F,
            sep = "\t",
            row.names = F)
sig_diff_write <- cbind(GeneSymbol=rownames(sig_diff), sig_diff)
write.table(sig_diff_write, file = "DEG_sig.xls",
            quote = F,
            sep = "\t",
            row.names = F)

### 火山图---------
#devtools::install_github("kongdd/Ipaper")
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(Ipaper)
library(scales)
library(ggrepel)

dat_rep<-DEG[rownames(DEG)%in%
               c(head(rownames(subset(sig_diff,sig_diff$log2FoldChange>3.95)),10),
                 head(rownames(subset(sig_diff,sig_diff$log2FoldChange< -3.8)),10)),]
volcano_plot<- ggplot(data = DEG, 
                      aes(x = log2FoldChange,
                          y = -log10(padj), 
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
    max.overlaps = 20,
    size = 4,
    box.padding = unit(0.5, "lines"),
    min.segment.length = 0,
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )+
  labs(x = "log (Fold Change)",
       y = "-log10 (adj.P.Val)")
volcano_plot
ggsave('volcano.png', volcano_plot,width = 8, height = 7)
ggsave('volcano.pdf', volcano_plot,width = 8, height = 7)
### 热图--------
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(gplots)
group_rt<-colData$group%>%as.data.frame()
rt<-gsea_exp
colnames(group_rt)<-'group'
rownames(group_rt)<-colData$sample
heat<-rt[rownames(rt)%in%
             c(head(rownames(subset(sig_diff,sig_diff$log2FoldChange>3.5)),10),
               head(rownames(subset(sig_diff,sig_diff$log2FoldChange< -3.5)),10)),]
x<-log2(heat+1)
ann_colors<-list(
  Group = c(Normal="lightblue",Tumor="darkorange"))
pheatmap(mat=x,
         annotation_col = group_rt,
         color=bluered(100),
         scale = "row",
         annotation_colors = ann_colors,
         fontsize = 10,
         show_colnames = F,
         cluster_cols = F,
         cluster_rows = T)

# 04 WGCNA------
setwd("/data/nas1/luchunlin/project/BJTC-292")
if (! dir.exists("./03_WGCNA")){
  dir.create("./03_WGCNA")
}
setwd("./03_WGCNA")
library(WGCNA)
library(reshape2)
library(stringr)
options(stringsAsFactors = F)
enableWGCNAThreads()
exprMat<-log2(gsea_exp+1)
dim(exprMat)
type = "unsigned"
corType = "pearson"
corFnc = ifelse(corType=="pearson", cor, bicor)
# 对二元变量，如样本性状信息计算相关性时，
# 或基因表达严重依赖于疾病状态时，需设置下面参数
maxPOutliers = ifelse(corType=="pearson",1,0.05)
# 关联样品性状的二元变量时，设置
robustY = ifelse(corType=="pearson",T,F)
dataExpr <- exprMat[rownames(gsea_exp),]
## 03-1 数据筛选-----
## 筛选中位绝对偏差（MAD）前75%的基因，至少MAD大于0.01
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
# [1] 151 10220

## 03-2 软阈值筛选----
## 样本聚类，查看是否有利群样本
sampleTree = hclust(dist(dataExpr), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
abline(h = 200, col = "red")

# Determine cluster under the line
# 剪枝算法，cutHeight修剪树枝的高度，minSize集群最小数
clust = cutreeStatic(sampleTree, cutHeight = 200, minSize = 10)
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
## 5
## 03-3一步法网络构建---------
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
# 0 (grey)表示未分入任何模块的基因。
table(net$colors)
##    0    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18   19   20   21   22 
## 1326 3013  660  659  475  438  411  397  309  280  271  269  258  249  236  209  201  150  139  104   95   40   31

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


## 03-5 绘制模块间的相关性热图-------
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


## 03-6 关联表型数据------
group_traits<-colData
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


## 03-7 筛选强相关性模块------------
## turquoise
module='turquoise'
probes=colnames(dataExpr)
inModule=(moduleColors==module)
modProbes=probes[inModule]
modGenes<-as.data.frame(modProbes)
colnames(modGenes)<-'modgene'
## 3013

# 05 MOD-DEGs-----
setwd("/data/nas1/luchunlin/project/BJTC-292")
if (! dir.exists("./04_MOD_DEGs")){
  dir.create("./04_MOD_DEGs")
}
setwd("./04_MOD_DEGs")
## 05-1 MOD_DEGs-----
## 取交集
diff_mod<-sig_diff[rownames(sig_diff)%in%modGenes$modgene,]
dim(diff_mod)
## 829
library(ggvenn)
mydata<-list(DEGs=rownames(sig_diff),modGenes=modGenes$modgene)
ggvenn(mydata,c('DEGs','modGenes'),
       fill_color = c("#B0E0E6", "#FF8888"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       set_name_color = c("#B0E0E6", "#FF8888"),
       text_color = 'black')
## 05-2 GO/KEGG------
#BiocManager::install('GOplot')
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(GOplot)
diff_gene_names <- rownames(diff_mod)
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
## GO 圈图
go_result<-read.table('GO.xls',header = T,sep = '\t',check.names = F)
go2=data.frame(Category='ALL',ID=go_result$ID,Term=go_result$Description,Genes=gsub("/", ", ", go_result$geneID), adj_pval = go_result$p.adjust)
idfc2<-diff_mod
genelist<-data.frame(ID=rownames(idfc2),logFC=idfc2$log2FoldChange)
rownames(genelist)=genelist[,1]
circ<-circle_dat(go2,genelist)
go_cir<-GOCircle(circ,rad1 = 2.5,rad2 = 3.5,label.size = 4,nsub = 10)
ggsave('GO_cir.png',go_cir,width = 12,height = 7)
ggsave('GO_cir.pdf',go_cir,width = 10,height = 6)
## GO 弦图
go_cir2<-cnetplot(ego,circular = TRUE, colorEdge = TRUE,showCategory = 5,
                  color_category = "#E5C494",color_gene = "#B3B3B3",
)
go_cir2
ggsave('GO_cir2.png',go_cir2,width = 20,height = 18)
ggsave('GO_cir2.pdf',go_cir2,width = 10,height = 6)
## KEGG富集分析
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "hsa",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk@result,file = "KEGG.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk, showCategory=25)
kk_dot
## 结果使用Cytoscape进行可视化

# 06 单因素Cox------
setwd("/data/nas1/luchunlin/project/BJTC-292")
if (! dir.exists("./05_univariate_cox")){
  dir.create("./05_univariate_cox")
}
setwd("./05_univariate_cox")
## 匹配生存数据
survival<-read_tsv(file = 'TCGA-LAML.survival.tsv')
survival<-survival[survival$sample%in%colnames(gsea_exp),]
## 132个匹配
train_data<-t(dataFilt.AML.final[,survival$sample])%>%as.data.frame()
train_data<-log2(train_data+1)
train_data<-train_data[,rownames(diff_mod)]
train_data$sample<-rownames(train_data)
train_data<-merge(survival,train_data,by='sample')
rownames(train_data)<-train_data$sample
train_data<-train_data[,-c(1,3)]

### 单因素cox
library(survival)
library(survminer)
colnames_sum <- colnames(train_data)
covariates <- colnames_sum[-which(colnames_sum %in% c("OS", "OS.time"))]
#Surv()函数产生一个生存对象  生存时间对生存的影响 对每一个变量构建生存分析公式
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste("Surv(OS.time, OS)~", x)))  #as.formula(). 将字符串转换成公式。构建formula对象
# coxph函数用于计算cox模型 循环对每一个特征做cox回归分析
univ_models <- lapply(univ_formulas,
                      function(x) {coxph(x, data = train_data)})

univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         #获取p值
                         p.value<-signif(x$wald["pvalue"], digits=3)
                         #获取HR
                         HR <-signif(x$coef[2], digits=3);
                         #获取95%置信区间
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"],3)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],3)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", 
                                      HR.confint.upper, ")")
                         res<-c(p.value,HR)
                         names(res)<-c("p.value","HR (95% CI for HR)")
                         return(res)
                       })

## coef是公式中的回归系数b（有时候也叫beta值）。exp(coef)是cox模型中的风险比（HR）
## z代表wald统计量，是coef除以其标准误se(coef)。ower .95 upper .95则是exp(coef)的95%置信区间，可信区间越窄，可信度越高，你的实验越精确，越是真理。
res_mod <- t(as.data.frame(univ_results, check.names = FALSE))
res_mod <- as.data.frame(res_mod)
res_results_0.05 <- res_mod[which(as.numeric(res_mod$p.value) < 0.05),]
res_results_0.05 <- na.omit(res_results_0.05)
write.table(res_results_0.05,
            file = "univariate_cox_result_0.05.xls",
            quote = F,
            sep = '\t',
            row.names = T)
dim(res_results_0.05)
## 278 2
library(tidyr)
res_results_0.05_2 <- separate(res_results_0.05, "HR (95% CI for HR)",
                               into = c("HR", "HR.95L", "HR.95H"),
                               sep = " ")
res_results_0.05_2 <- separate(res_results_0.05_2, "HR.95L",
                               into = c("HR.95L", "HR.95H"),
                               sep = "\\-")
res_results_0.05_2$HR.95L <- gsub("\\(", "", res_results_0.05_2$HR.95L)
res_results_0.05_2$HR.95H <- gsub("\\)", "", res_results_0.05_2$HR.95H)

res_results_0.05_2[,1:ncol(res_results_0.05_2)] <- as.numeric(unlist(res_results_0.05_2[,1:ncol(res_results_0.05_2)]))
res_results_0.05_2 <- res_results_0.05_2[order(res_results_0.05_2$p.value),]
res_results_0.05_2<-res_results_0.05_2[c(1:20),]
res_results_0.05_2 <- res_results_0.05_2[order(res_results_0.05_2$HR),]
hz <- paste(round(res_results_0.05_2$HR,3),
            "(",round(res_results_0.05_2$HR.95L,3),
            "-",round(res_results_0.05_2$HR.95H,3),")",sep = "")
tabletext <- cbind(c(NA,"Gene",rownames(res_results_0.05_2)),
                   c(NA,"P value",ifelse(res_results_0.05_2$p.value<0.001,
                                         "< 0.001",
                                         round(res_results_0.05_2$p.value,3))),
                   c(NA,"Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1
library(forestplot)
pdf(file = "univariate_cox_forest.pdf", height = 10, width = 10, onefile = F)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE, TRUE,rep(FALSE, 300)),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,NA,res_results_0.05_2$HR),
           lower=c(NA,NA,res_results_0.05_2$HR.95L), #95%置信区间下限
           upper=c(NA,NA,res_results_0.05_2$HR.95H), #95%置信区间上限
           boxsize=0.2,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(0, 0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(1,"cm"), #固定行高
           graphwidth = unit(.5,"npc"), #图在表中的宽度比例
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines = list("3" = gpar(col = "black", lty = 1, lwd = 2)),
           # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
           #                 "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
           #                 "59" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
           # mar=unit(rep(0.01, times = 4), "cm"),#图形页边距
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=1, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Hazard Ratio",
           grid = T) # 垂直于x轴的网格线，对应每个刻度
dev.off()
png(filename = "univariate_cox_forest.png", height = 700, width = 800)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE, TRUE,rep(FALSE, 70)),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,NA,res_results_0.05_2$HR),
           lower=c(NA,NA,res_results_0.05_2$HR.95L), #95%置信区间下限
           upper=c(NA,NA,res_results_0.05_2$HR.95H), #95%置信区间上限
           boxsize=0.2,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(0, 0.2,0.4,0.6,0.8,1,1.2,1.4,1.6,1.8,2,2.2), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(1,"cm"), #固定行高
           graphwidth = unit(.5,"npc"), #图在表中的宽度比例
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines = list("3" = gpar(col = "black", lty = 1, lwd = 2)),
           # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
           #                 "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
           #                 "59" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
           # mar=unit(rep(0.01, times = 4), "cm"),#图形页边距
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=1, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Hazard Ratio",
           grid = T) # 垂直于x轴的网格线，对应每个刻度

dev.off()
# 07 Lasso 回归------
setwd("/data/nas1/luchunlin/project/BJTC-292")
if (! dir.exists("./06_Lasso")){
  dir.create("./06_Lasso")
}
setwd("./06_Lasso")
library(glmnet)
x_all <- subset(train_data, select = -c(OS, OS.time))
x_all <- x_all[,rownames(res_results_0.05)]
y_all <- subset(train_data, select = c(OS, OS.time))

# 拟合模型
fit <- glmnet(as.matrix(x_all), Surv(y_all$OS.time,y_all$OS), 
              family = "cox") 
#dev.new()
png(filename = "lasso_model.png", height = 400, width = 500)
plot(fit, xvar = "lambda",label = TRUE, las=1)
dev.off()
pdf(file = "lasso_model.pdf", height = 5)
plot(fit, xvar = "lambda",label = TRUE, las=1)
dev.off()
# 交叉验证拟合模型
set.seed(43)
cvfit = cv.glmnet(as.matrix(x_all),
                  Surv(y_all$OS.time,y_all$OS),nfold=10,
                  family = "cox") 

png(filename = "lasso_verify.png", height = 400, width = 500)
plot(cvfit, las =1)
dev.off()
pdf(file = "lasso_verify.pdf", height = 5)
plot(cvfit, las =1)
dev.off()
cvfit
# 提取指定lambda时特征的系数
coef.min = coef(cvfit, s = "lambda.min")  ## lambda.min & lambda.1se 取一个
cvfit$lambda.min
# [1] 0.07122514
# 找出那些回归系数没有被惩罚为0的
active.min = which(coef.min@i != 0)
# 提取基因名称
lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i+1]
lasso_geneids
## 29
##  [1] "TNFSF15"   "GNGT2"     "TMEM217"   "CCL3"      "FMN1"      "LGALS1"    "IGLON5"    "DNAH3"     "ACSM1"     "HTR7"      "PLXNB1"   
## [12] "TMPRSS11D" "CCNJL"     "NRGN"      "MAP7D2"    "HHLA2"     "KIF26A"    "CASKIN1"   "AVPR1B"    "HRH1"      "EPS8"      "GPR56"    
## [23] "SLC5A5"    "C14orf37"  "ABLIM3"    "BVES"      "CCL23"     "GABRD"     "TRNP1"    
write(lasso_geneids, "lasso_genes.csv")
write.csv(x_all,file = "Lasso_x.csv",quote = F)
write.csv(y_all,file = "Lasso_y.csv",quote = F)

# 08 Risk Model------
setwd("/data/nas1/luchunlin/project/BJTC-292")
if (! dir.exists("./07_risk")){
  dir.create("./07_risk")
}
setwd("./07_risk")
## 08-1 riskscore-----

riskScore=predict(cvfit,newx = as.matrix(x_all),s=cvfit$lambda.min)
riskScore<-as.numeric(riskScore)
coxGene=lasso_geneids
outCol=c("OS","OS.time",coxGene)
risk=as.vector(ifelse(riskScore>median(riskScore),0,1))
risk <- as.data.frame(c(cbind(id=rownames(cbind(train_data[,outCol],
                                                riskScore,
                                                risk)),
                              cbind(train_data[,outCol],
                                    riskScore,
                                    risk))))
table(risk$risk)
##  0  1 
## 66 66 
risk_dis <- ggplot(risk, aes(x=reorder(id, riskScore), 
                             y=riskScore, 
                             color = factor(risk, 
                                            levels = c(0, 1), 
                                            labels = c("High Risk", "Low Risk")))) +
  geom_point() +
  scale_color_manual(values = c("#A73030FF", "#0073C2FF")) + 
  scale_x_discrete(breaks = risk[order(risk$riskScore),]$id[c(1,100,200,300,400,500,600,700,800,900)],
                   labels = c(1,100,200,300,400,500,600,700,800,900),
                   expand = c(0.02,0)) +
  geom_vline(xintercept = nrow(risk[which(risk$risk==1),]) + 0.5,
             lty = 2) +
  geom_hline(yintercept = median(riskScore),
             lty =2) +
  labs(x = "Patients(increasing risk score)",
       y = "Risk Score",
       title = "Train Risk Score Distribution") + 
  theme_base() +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.justification = c(0,1),
        legend.position = c(0,1),
        #        legend.margin = margin(c(-5,4,4,3)),
        legend.background = element_rect(color = "black", size = .3),
        plot.title = element_text(size = 15, hjust = 0.5))
risk_dis
surv_stat <- ggplot(risk, aes(x=reorder(id, riskScore),
                              y=OS.time/365,
                              color = factor(OS,
                                             levels = c(0,1),
                                             labels = c("Alive", "Dead")))) +
  geom_point() +
  scale_color_manual(values = c("#0073C2FF","#A73030FF")) +
  scale_x_discrete(breaks = risk[order(risk$riskScore),]$id[c(1,100,200,300,400,500,600,700,800,900)],
                   labels = c(1,100,200,300,400,500,600,700,800,900),
                   expand = c(0.02,0)) +
  ylim(c(0,15))+
  geom_vline(xintercept = nrow(risk[which(risk$risk==1),]) + 0.5,
             lty = 2) +
  labs(x = "Patients(increasing risk score)",
       y = "Survival time (years)",
       title = "Train survival state distribution") + 
  theme_base() +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.justification = c(0,1),
        legend.position = c(0,1),
        #        legend.margin = margin(c(-5,4,4,3)),
        legend.background = element_rect(color = "black", size = .3),
        plot.title = element_text(size = 15, hjust = 0.5))

surv_stat

## 08-2 KM曲线----------
kmfit<-survfit(Surv(OS.time, OS) ~ risk, data =  risk)
train_survival_median <- ggsurvplot(kmfit,
                                    pval = TRUE, 
                                    conf.int = F,
                                    legend.labs=c("High risk","Low risk" ),
                                    legend.title="Risk score",
                                    title="Train KM",
                                    font.main = c(15,"bold"),
                                    risk.table = TRUE, 
                                    risk.table.col = "strata", 
                                    linetype = "strata", 
                                    surv.median.line = "hv", 
                                    ggtheme = theme_bw(), 
                                    palette = c("#A73030FF", "#0073C2FF"))
train_survival_median

## 08-3 ROC曲线-------- 
# BiocManager::install('survivalROC')
library(survivalROC)
multi_ROC <- function(time_vector, risk_score_table){
  library(survivalROC)
  single_ROC <- function(single_time){
    for_ROC <- survivalROC(Stime=risk_score_table$OS.time,
                           status=risk_score_table$OS,
                           marker=risk_score_table$riskScore,
                           predict.time=single_time,
                           method = 'KM')
    data.frame('True_positive'=for_ROC$TP, 'False_positive'=for_ROC$FP,
               'Cut_values'=for_ROC$cut.values, 'Time'=rep(single_time, length(for_ROC$TP)),
               'AUC'=rep(for_ROC$AUC, length(for_ROC$TP)))
  }
  multi_ROC_list <- lapply(time_vector, single_ROC)
  do.call(rbind, multi_ROC_list)
}
for_multi_ROC <- multi_ROC(time_vector = c(365*seq(1,5,2)), 
                           risk_score_table = risk)
for_multi_ROC$Time <- factor(for_multi_ROC$Time)
#devtools::install_github('yikeshu0611/geomROC')
library(scales)
library(geomROC)
library(plotROC)
library(ggthemes)
auc_y1 <- round(for_multi_ROC[which(for_multi_ROC$Time==365),5][1],2)
auc_y3 <- round(for_multi_ROC[which(for_multi_ROC$Time==1095),5][1],2)
auc_y5 <- round(for_multi_ROC[which(for_multi_ROC$Time==1825),5][1],2)

ROC <- ggplot(for_multi_ROC, aes(x=False_positive,
                                 y=True_positive, 
                                 label=Cut_values, 
                                 color=Time)) + 
  scale_color_manual(breaks = c("365", "1095", "1825"),
                     labels = c("1 years", "3 years", "5 years"),
                     values = c("#4682B4", "#FF4040", "#20B2AA")) +
  geom_roc(labels = F, stat = 'identity') + 
  style_roc() + 
  geom_abline(slope = 1, intercept = 0, color = 'gray', linetype=2) +
  theme_bw() +
  labs(title = "Train ROC") +
  theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
        panel.grid = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold")) +
  # annotate('text', x=.75, y=.25, label=paste('AUC of 1 years =', round(auc_y1,2))) + 
  # annotate('text', x=.75, y=.15, label=paste('AUC of 3 years =', round(auc_y3,2))) + 
  # annotate('text', x=.75, y=.05, label=paste('AUC of 5 years =', round(auc_y5,2))) +
  annotate("text", x=0.75, y=c(0.25, 0.15, 0.05),
           label = c(paste('AUC of 1 years =', format(auc_y1,nsmall=2)),
                     paste('AUC of 3 years =', format(auc_y3,nsmall=2)),
                     paste('AUC of 5 years =', format(auc_y5,nsmall=2))))
ROC
ggsave('Train ROC.png', ROC,width = 5, height = 4)
ggsave('Train ROC.pdf', ROC,width = 5, height = 4)
# 09 外部数据库验证------
setwd("/data/nas1/luchunlin/project/BJTC-292")
if (! dir.exists("./08_External_va")){
  dir.create("./08_External_va")
}
setwd("./08_External_va")
## 2个验证集 TARGET & GSE12417
## 09-1 TRAGET-AML-----
test_data<-t(dat_va1)%>%as.data.frame()
survival_va<-read_tsv(file = 'TARGET-AML.survival.tsv')
test_data$sample<-rownames(test_data)
test_data<-merge(survival_va,test_data,by='sample')
rownames(test_data)<-test_data$sample

test_data<-test_data[,-c(1,3)]
test_data<-test_data[,colnames(test_data)%in%colnames(train_data)]
# 开始验证
x_all_out <- subset(test_data, select = -c(OS, OS.time))
x_all_out <- x_all_out[,rownames(res_results_0.05)]
riskScore_out=predict(cvfit,newx = as.matrix(x_all_out),s=cvfit$lambda.min)
riskScore_out<-as.numeric(riskScore_out)

risk_out=as.vector(ifelse(riskScore_out>median(riskScore_out),0,1))
risk_out <- as.data.frame(c(cbind(id=rownames(cbind(test_data[,outCol],
                                                    riskScore_out,
                                                    risk_out)),
                                  cbind(test_data[,outCol],
                                        riskScore_out,
                                        risk_out))))

library(ggplot2)
library(ggthemes)
median(riskScore_out)

risk_dis_out <- ggplot(risk_out, aes(x=reorder(id, riskScore_out), 
                                     y=riskScore_out, 
                                     color = factor(risk_out, 
                                                    levels = c(0, 1), 
                                                    labels = c("High Risk", "Low Risk")))) +
  geom_point() +
  scale_color_manual(values = c("#A73030FF", "#0073C2FF")) + 
  scale_x_discrete(breaks = risk_out[order(risk_out$riskScore_out),]$id[c(1,25,50,75,100,125,150,175,200)],
                   labels = c(1,25,50,75,100,125,150,175,200),
                   expand = c(0.02,0)) +
  geom_vline(xintercept = nrow(risk_out[which(risk_out$risk_out==1),]) + 0.5,
             lty = 2) +
  geom_hline(yintercept = median(riskScore_out),
             lty =2) +
  labs(x = "Patients(increasing risk score)",
       y = "Risk Score",
       title = "Validation Risk Score Distribution") + 
  theme_base() +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.justification = c(0,1),
        legend.position = c(0,1),
        #        legend.margin = margin(c(-5,4,4,3)),
        legend.background = element_rect(color = "black", size = .3),
        plot.title = element_text(size = 15, hjust = 0.5))
risk_dis_out
surv_stat_out <- ggplot(risk_out, aes(x=reorder(id, riskScore_out),
                                      y=OS.time/365,
                                      color = factor(OS,
                                                     levels = c(0,1),
                                                     labels = c("Alive", "Dead")))) +
  geom_point() +
  scale_color_manual(values = c("#0073C2FF", "#A73030FF")) +
  scale_x_discrete(breaks = risk_out[order(risk_out$riskScore_out),]$id[c(1,25,50,75,100,125,150,175,200)],
                   labels = c(1,25,50,75,100,125,150,175,200),
                   expand = c(0.02,0)) +
  ylim(x=c(0,15)) +
  geom_vline(xintercept = nrow(risk_out[which(risk_out$risk_out==1),]) + 0.5,
             lty = 2) +
  labs(x = "Patients(increasing risk score)",
       y = "Survival time (years)",
       title = "Validation Survival State Distribution") + 
  theme_base() +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.justification = c(0,1),
        legend.position = c(0,1),
        #        legend.margin = margin(c(-5,4,4,3)),
        legend.background = element_rect(color = "black", size = .3),
        plot.title = element_text(size = 15, hjust = 0.5))
surv_stat_out
kmfit_out <- survfit(Surv(OS.time, OS) ~ risk_out, data =  risk_out)
verify_survival_median <- ggsurvplot(kmfit_out,
                                     pval = TRUE, 
                                     conf.int = F,
                                     legend.labs=c("High risk","Low risk" ),
                                     legend.title="Risk score",
                                     title="Validition KM",
                                     font.main = c(15,"bold"),
                                     risk.table = TRUE, 
                                     risk.table.col = "strata", 
                                     linetype = "strata", 
                                     surv.median.line = "hv", 
                                     ggtheme = theme_bw(), 
                                     palette = c("#A73030FF", "#0073C2FF"))
verify_survival_median

multi_ROC <- function(time_vector, risk_score_table){
  library(survivalROC)   
  single_ROC <- function(single_time){
    for_ROC <- survivalROC(Stime=risk_score_table$OS.time,
                           status=risk_score_table$OS,
                           marker=risk_score_table$riskScore_out,
                           predict.time=single_time,method = 'KM')
    data.frame('True_positive'=for_ROC$TP, 'False_positive'=for_ROC$FP,
               'Cut_values'=for_ROC$cut.values, 'Time'=rep(single_time, length(for_ROC$TP)),
               'AUC'=rep(for_ROC$AUC, length(for_ROC$TP)))
  }
  multi_ROC_list <- lapply(time_vector, single_ROC)
  do.call(rbind, multi_ROC_list)
}

for_multi_ROC <- multi_ROC(time_vector = c(365*seq(1,5,2)), 
                           risk_score_table = risk_out)
for_multi_ROC$Time <- factor(for_multi_ROC$Time)

# 画ROC曲线 
library(scales)
library(geomROC)
library(plotROC)
auc_y1 <- round(for_multi_ROC[which(for_multi_ROC$Time==365),5][1],2)
auc_y3 <- round(for_multi_ROC[which(for_multi_ROC$Time==1095),5][1],2)
auc_y5 <- round(for_multi_ROC[which(for_multi_ROC$Time==1825),5][1],2)


ROC <- ggplot(for_multi_ROC, aes(x=False_positive,
                                 y=True_positive, 
                                 label=Cut_values, 
                                 color=Time)) + 
  scale_color_manual(breaks = c("365", "1095", "1825"),
                     labels = c("1 years", "3 years", "5 years"),
                     values = c("#4682B4", "#FF4040", "#20B2AA")) +
  geom_roc(labels = F, stat = 'identity') + 
  style_roc() + 
  geom_abline(slope = 1, intercept = 0, color = 'gray', linetype=2) +
  theme_bw() +
  labs(title = "Train ROC") +
  theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
        panel.grid = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold")) +
  # annotate('text', x=.75, y=.25, label=paste('AUC of 1 years =', round(auc_y1,2))) + 
  # annotate('text', x=.75, y=.15, label=paste('AUC of 3 years =', round(auc_y3,2))) + 
  # annotate('text', x=.75, y=.05, label=paste('AUC of 5 years =', round(auc_y5,2))) +
  annotate("text", x=0.75, y=c(0.25, 0.15, 0.05),
           label = c(paste('AUC of 1 years =', format(auc_y1,nsmall=2)),
                     paste('AUC of 3 years =', format(auc_y3,nsmall=2)),
                     paste('AUC of 5 years =', format(auc_y5,nsmall=2))))
ROC
ggsave('Test ROC1.png', ROC,width = 5, height = 4)
ggsave('Test ROC1.pdf', ROC,width = 5, height = 4)
## 09-2 GSE12417-----
survival_va2<-data.frame(sample=pd_va$geo_accession,
                         OS=pd_va$`event (1:ch1`,
                         OS.time=as.numeric(pd_va$`overall survival (months):ch1`)*30)
survival_va2$OS<-gsub('dead, 0:alive): ','',survival_va2$OS,fixed = T)
survival_va2$OS<-as.numeric(survival_va2$OS)
test_data2<-t(dat_va2)%>%as.data.frame()
lasso_geneids<-as.data.frame(lasso_geneids)
unmapgene<-lasso_geneids[!lasso_geneids$lasso_geneids%in%colnames(test_data2),]
unmapgene
## ADGRG1  ARMH4 分别为GPR56和C14orf37的别称。需要加回去
unmap_exp<-test_data2[,c("ADGRG1","ARMH4")]
colnames(unmap_exp)<-c('GPR56','C14orf37')
test_data2<-cbind(unmap_exp,test_data2)
test_data2$sample<-rownames(test_data2)
test_data2<-merge(survival_va2,test_data2,by='sample')
rownames(test_data2)<-test_data2$sample
test_data2<-test_data2[,-1]
# 开始验证
##手动算
test_data3<-test_data2[,lasso_geneids$lasso_geneids]
risk_out2<-data.frame(test_data3)
risk_out2$risk_out<-NA
risk_out2$riskScore_out<-NA
cnt<-1
coef.min<-coef.min[lasso_geneids$lasso_geneids,]

while (cnt < 105) {
  risk_out2$riskScore_out[cnt]<-sum(coef.min*test_data3[cnt,])
  cnt = cnt + 1
}

dim(test_data3)
cnt<-1
while (cnt < 105) {
  risk_out2$risk_out[cnt]=as.vector(ifelse(risk_out2$riskScore_out[cnt]>median(risk_out2$riskScore_out),0,1))
  cnt = cnt + 1
}
riskScore_out2<-as.numeric(risk_out2$riskScore_out)
risk_out2=as.vector(ifelse(riskScore_out2>median(riskScore_out2),0,1))
risk_out2 <- as.data.frame(c(cbind(id=rownames(cbind(test_data2[,outCol],
                                                    riskScore_out2,
                                                    risk_out2)),
                                  cbind(test_data2[,outCol],
                                        riskScore_out2,
                                        risk_out2))))

library(ggplot2)
library(ggthemes)
median(riskScore_out2)

risk_dis_out <- ggplot(risk_out2, aes(x=reorder(id, riskScore_out2), 
                                     y=riskScore_out2, 
                                     color = factor(risk_out2, 
                                                    levels = c(0, 1), 
                                                    labels = c("High Risk", "Low Risk")))) +
  geom_point() +
  scale_color_manual(values = c("#A73030FF", "#0073C2FF")) + 
  scale_x_discrete(breaks = risk_out2[order(risk_out2$riskScore_out2),]$id[c(1,10,20,30,40,50,60,70,80,90,100,110)],
                   labels = c(1,10,20,30,40,50,60,70,80,90,100,110),
                   expand = c(0.02,0)) +
  geom_vline(xintercept = nrow(risk_out2[which(risk_out2$risk_out2==1),]) + 0.5,
             lty = 2) +
  geom_hline(yintercept = median(riskScore_out2),
             lty =2) +
  labs(x = "Patients(increasing risk score)",
       y = "Risk Score",
       title = "Validation Risk Score Distribution") + 
  theme_base() +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.justification = c(0,1),
        legend.position = c(0,1),
        #        legend.margin = margin(c(-5,4,4,3)),
        legend.background = element_rect(color = "black", size = .3),
        plot.title = element_text(size = 15, hjust = 0.5))
risk_dis_out
surv_stat_out <- ggplot(risk_out2, aes(x=reorder(id, riskScore_out2),
                                      y=OS.time/365,
                                      color = factor(OS,
                                                     levels = c(0,1),
                                                     labels = c("Alive", "Dead")))) +
  geom_point() +
  scale_color_manual(values = c("#0073C2FF", "#A73030FF")) +
  scale_x_discrete(breaks = risk_out2[order(risk_out2$riskScore_out2),]$id[c(1,10,20,30,40,50,60,70,80,90,100,110)],
                   labels = c(1,10,20,30,40,50,60,70,80,90,100,110),
                   expand = c(0.02,0)) +
  ylim(x=c(0,20)) +
  geom_vline(xintercept = nrow(risk_out2[which(risk_out2$risk_out2==1),]) + 0.5,
             lty = 2) +
  labs(x = "Patients(increasing risk score)",
       y = "Survival time (years)",
       title = "Validation Survival State Distribution") + 
  theme_base() +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.justification = c(0,1),
        legend.position = c(0,1),
        #        legend.margin = margin(c(-5,4,4,3)),
        legend.background = element_rect(color = "black", size = .3),
        plot.title = element_text(size = 15, hjust = 0.5))
surv_stat_out
kmfit_out <- survfit(Surv(OS.time, OS) ~ risk_out2, data =  risk_out2)
verify_survival_median <- ggsurvplot(kmfit_out,
                                     pval = TRUE, 
                                     conf.int = F,
                                     legend.labs=c("High risk","Low risk" ),
                                     legend.title="Risk score",
                                     title="Validition KM",
                                     font.main = c(15,"bold"),
                                     risk.table = TRUE, 
                                     risk.table.col = "strata", 
                                     linetype = "strata", 
                                     surv.median.line = "hv", 
                                     ggtheme = theme_bw(), 
                                     palette = c("#A73030FF", "#0073C2FF"))
verify_survival_median

multi_ROC <- function(time_vector, risk_score_table){
  library(survivalROC)   
  single_ROC <- function(single_time){
    for_ROC <- survivalROC(Stime=risk_score_table$OS.time,
                           status=risk_score_table$OS,
                           marker=risk_score_table$riskScore_out,
                           predict.time=single_time,method = 'KM')
    data.frame('True_positive'=for_ROC$TP, 'False_positive'=for_ROC$FP,
               'Cut_values'=for_ROC$cut.values, 'Time'=rep(single_time, length(for_ROC$TP)),
               'AUC'=rep(for_ROC$AUC, length(for_ROC$TP)))
  }
  multi_ROC_list <- lapply(time_vector, single_ROC)
  do.call(rbind, multi_ROC_list)
}

for_multi_ROC <- multi_ROC(time_vector = c(365*seq(1,5,2)), 
                           risk_score_table = risk_out2)
for_multi_ROC$Time <- factor(for_multi_ROC$Time)

# 画ROC曲线 
auc_y1 <- round(for_multi_ROC[which(for_multi_ROC$Time==365),5][1],2)
auc_y3 <- round(for_multi_ROC[which(for_multi_ROC$Time==1095),5][1],2)
auc_y5 <- round(for_multi_ROC[which(for_multi_ROC$Time==1825),5][1],2)


ROC <- ggplot(for_multi_ROC, aes(x=False_positive,
                                 y=True_positive, 
                                 label=Cut_values, 
                                 color=Time)) + 
  scale_color_manual(breaks = c("365", "1095", "1825"),
                     labels = c("1 years", "3 years", "5 years"),
                     values = c("#4682B4", "#FF4040", "#20B2AA")) +
  geom_roc(labels = F, stat = 'identity') + 
  style_roc() + 
  geom_abline(slope = 1, intercept = 0, color = 'gray', linetype=2) +
  theme_bw() +
  labs(title = "Test ROC") +
  theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
        panel.grid = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold")) +
  # annotate('text', x=.75, y=.25, label=paste('AUC of 1 years =', round(auc_y1,2))) + 
  # annotate('text', x=.75, y=.15, label=paste('AUC of 3 years =', round(auc_y3,2))) + 
  # annotate('text', x=.75, y=.05, label=paste('AUC of 5 years =', round(auc_y5,2))) +
  annotate("text", x=0.75, y=c(0.25, 0.15, 0.05),
           label = c(paste('AUC of 1 years =', format(auc_y1,nsmall=2)),
                     paste('AUC of 3 years =', format(auc_y3,nsmall=2)),
                     paste('AUC of 5 years =', format(auc_y5,nsmall=2))))
ROC
ggsave('Test ROC2.png', ROC,width = 5, height = 4)
ggsave('Test ROC2.pdf', ROC,width = 5, height = 4)
# 10 风险评分与临床相关性分析-------
setwd("/data/nas1/luchunlin/project/BJTC-292")
if (! dir.exists("./09_clinical_index")){
  dir.create("./09_clinical_index")
}
setwd("./09_clinical_index")
## 要求M0-M7必须有
phenotype<-read_tsv(file = 'TCGA-LAML.GDC_phenotype.tsv')
colnames(phenotype)
train_phenotype<-data.frame(sample=phenotype$submitter_id.samples,
                             gender=phenotype$gender.demographic,
                             FAB_classification=phenotype$leukemia_french_american_british_morphology_code)
train_phenotype<-merge(train_phenotype,survival,by='sample')
train_phenotype<-train_phenotype[,-5]
train_phenotype$FAB_classification<-gsub(' Undifferentiated','',train_phenotype$FAB_classification,fixed = T)
train_phenotype$FAB_classification<-gsub('Not Classified',NA,train_phenotype$FAB_classification)

colnames(train_phenotype)<-c('id','gender','FAB_classification','OS','OS.time')
sub_risk <- subset(risk, select = c(id, riskScore))
train_phenotype2 <- merge(train_phenotype,
                          sub_risk,
                          by = "id")
write.table(train_phenotype2,
            file = "clinical_risk.csv",
            row.names = T,
            sep = "\t",
            quote = F)
library(ggpubr)
library(Ipaper)
library(ggthemes)
## 10-1 gender-----
my_comparisons <- list(c("male", "female"))
gender<-ggplot(train_phenotype2,aes(x = gender, y = riskScore, fill = gender)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "riskScore",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("Gender") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              map_signif_level = T,
              y_position = c(7))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 20))+
  guides(fill='none')
gender
## 10-2 FAB_classification-----
my_comparisons <- list(c("M0","M1"),c("M1","M2"),c("M2","M3"),c("M3","M4"),c("M4","M5"),c("M5","M6"),c("M6","M7"))
stage_data <- data.frame(riskScore = train_phenotype2$riskScore,
                         stage = factor(train_phenotype2$FAB_classification,
                                        levels = c("M0","M1","M2","M3","M4","M5","M6","M7")))
stage_data <- na.omit(stage_data)
stage<-ggplot(stage_data,aes(x = stage, y = riskScore, fill = stage)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "riskScore",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("FAB classification") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              map_signif_level = T,
              y_position = c(7,8,9,7,8,9,7))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 20))+
  guides(fill='none')
stage
library(patchwork)
all_clinical_index <- stage + gender + 
  plot_layout(ncol = 2) & 
  theme_bw() & 
  theme(legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold",size = 12))
all_clinical_index
# 11 GSEA分析---------
setwd("/data/nas1/luchunlin/project/BJTC-292")
if (! dir.exists("./10_GSEA")){
  dir.create("./10_GSEA")
}
setwd("./10_GSEA")
risk2 <- risk
risk2$risk_label <- ifelse(risk$risk == 0, "High", "Low")
gsea_exp<-dataFilt.AML.final[,risk2$id]
all(colnames(gsea_exp) == risk2$id)
dim(gsea_exp)
write.table(gsea_exp,
            file = "gsea_exp.xls",
            quote = F,
            sep = "\t",
            row.names = T)
colData<-data.frame(sample=risk2$id,
                    group=risk2$risk_label)

colData$group<-factor(colData$group,levels = c('Low','High'))
write.table(colData,
            file = "group_gsea.xls",
            quote = F,
            row.names = F)

dds<-DESeqDataSetFromMatrix(countData = gsea_exp,colData=colData,design = ~group)
dds = dds[rownames(counts(dds)) > 1,]
dds<-DESeq(dds)
dds
res =results(dds, contrast = c("group",'Low','High'))
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
DEGeneSets <- DEGeneSets[order(DEGeneSets$padj),]
dim(DEGeneSets)
set.seed(1)
KEGG_GSEA <- GSEA(geneList, TERM2GENE = kegg_set, pvalueCutoff = 0.05)
KEGG_result <- KEGG_GSEA@result
dim(KEGG_result)
gseaplot2(KEGG_GSEA,c(1:5),color = c('#7B68EE','#CD3333','#20B2AA','#FF8C00','#FF6666'),
          title = 'KEGG GSEA',
          base_size = 11,
          rel_heights = c(1.5, 0.3, 0.5))

# 12 免疫评分、基质评分----------
setwd("/data/nas1/luchunlin/project/BJTC-292")
if (! dir.exists("./11_estimate")){
  dir.create("./11_estimate")
}
setwd("./11_estimate")
library(estimate)
expr_train <- log2(gsea_exp+1)
write.table(expr_train, 
            'expr_log2.txt', 
            col.names = T, 
            row.names = T, 
            quote = F, sep="\t")
# 生成expr_train.gct
filterCommonGenes(input.f = './expr_log2.txt', 
                  output.f = 'expr_train.gct', 
                  id = 'GeneSymbol')
# 生成train_purity.gct
estimateScore('expr_train.gct', 'train_purity.gct', platform="affymetrix")
es_score <- read.table('train_purity.gct', skip = 2, header = T)
immu_score <- es_score[,3:length(es_score)]
rownames(immu_score) <- es_score$NAME
immu_score<-t(immu_score)
immu_score<-as.data.frame(immu_score)
immu_score$sample<-colnames(expr_train)
rownames(immu_score)<-immu_score$sample
write.table(es_score,
            file = "es_score.xls",
            sep = "\t",
            quote = F,
            row.names = F)
violin_dat <- immu_score
high_sample <- colData$sample[which(colData$group == 'High')]
low_sample <- colData$sample[which(colData$group == 'Low')]
length(high_sample)
# [1] 66
length(low_sample)
# [1] 66
violin_dat$group <- ifelse(violin_dat$sample %in% high_sample,
                           "High", "Low")
library(tidyr)
library(ggplot2)
library(ggpubr)
library(Ipaper)

{
  p1 <- ggplot(violin_dat, aes(x=group,y=StromalScore, fill=group))+
    geom_violin() + #绘制小提琴图
    stat_boxplot(geom="errorbar",
                 width=0.1,
                 position = position_dodge(0.9)) +
    geom_boxplot(width=0.5,
                 position=position_dodge(0.9),
                 outlier.shape = NA,
                 fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
    geom_point(aes(fill = group),
               size = 0.5,
               position = position_dodge(0.9))+
    scale_fill_manual(values = c("#FF6A6A", "#20B2AA"), name = "Group") + 
    stat_compare_means(aes(group = group),method = 'wilcox')+ 
    theme_bw()+ #背景变为白色
    ylim(-3000,2000) +
    theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
          panel.grid = element_blank()) + 
    theme(axis.text.x=element_text(family="Times",size=12,face="bold"), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
          axis.text.y=element_text(family="Times",size=12,face="bold"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
          axis.title.y=element_text(family="Times",size = 15,face="bold"), #设置y轴标题的字体属性
          axis.title.x=element_text(family="Times",size = 10,face="bold"),
          plot.title = element_text(hjust = 0.5,family="Times",size = 15,face="bold"),
          legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold"))+
    labs(title="Stromal Score", x="", y="Stromal Score")
  p1
  p2 <- ggplot(violin_dat, aes(x=group,y=ImmuneScore, fill=group))+
    geom_violin() + #绘制小提琴图
    stat_boxplot(geom="errorbar",
                 width=0.1,
                 position = position_dodge(0.9)) +
    geom_boxplot(width=0.5,
                 position=position_dodge(0.9),
                 outlier.shape = NA,
                 fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
    geom_point(aes(fill = group),
               size = 0.5,
               position = position_dodge(0.9))+
    scale_fill_manual(values = c("#FF6A6A", "#20B2AA"), name = "Group") +
    stat_compare_means(aes(group = group),method = 'wilcox')+ 
    theme_bw()+ #背景变为白色
    ylim(0,4000) +
    theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
          panel.grid = element_blank()) + 
    theme(axis.text.x=element_text(family="Times",size=12,face="bold"), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
          axis.text.y=element_text(family="Times",size=12,face="bold"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
          axis.title.y=element_text(family="Times",size = 15,face="bold"), #设置y轴标题的字体属性
          axis.title.x=element_text(family="Times",size = 10,face="bold"),
          plot.title = element_text(hjust = 0.5,family="Times",size = 15,face="bold"),     
          legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold"))+
    labs(title="Immune Score", x="", y="Immune Score")
  p2
  p3 <- ggplot(violin_dat, aes(x=group, y=ESTIMATEScore, fill=group))+
    geom_violin() + #绘制小提琴图
    stat_boxplot(geom="errorbar",
                 width=0.1,
                 position = position_dodge(0.9)) +
    geom_boxplot(width=0.5,
                 position=position_dodge(0.9),
                 outlier.shape = NA,
                 fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
    geom_point(aes(fill = group),
               size = 0.5,
               position = position_dodge(0.9))+
    scale_fill_manual(values = c("#FF6A6A", "#20B2AA"), name = "Group") +
    stat_compare_means(aes(group = group),method = 'wilcox')+ 
    theme_bw()+ #背景变为白色
    theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
          panel.grid = element_blank()) + 
    ylim(-3000, 6000) +
    theme(axis.text.x=element_text(family="Times",size=12,face="bold"), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
          axis.text.y=element_text(family="Times",size=12,face="bold"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
          axis.title.y=element_text(family="Times",size = 15,face="bold"), #设置y轴标题的字体属性
          axis.title.x=element_text(family="Times",size = 10,face="bold"),
          plot.title = element_text(hjust = 0.5,family="Times",size = 15,face="bold"),
          legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold"))+
    labs(title="ESTIMATE Score", x="", y="ESTIMATE Score")
  p3
  p4 <- ggplot(violin_dat, aes(x=group, y=TumorPurity, fill=group))+
    geom_violin() + #绘制小提琴图
    stat_boxplot(geom="errorbar",
                 width=0.1,
                 position = position_dodge(0.9)) +
    geom_boxplot(width=0.5,
                 position=position_dodge(0.9),
                 outlier.shape = NA,
                 fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
    geom_point(aes(fill = group),
               size = 0.5,
               position = position_dodge(0.9))+
    scale_fill_manual(values = c("#FF6A6A", "#20B2AA"), name = "Group") +
    stat_compare_means(aes(group = group),method = 'wilcox')+ 
    theme_bw()+ #背景变为白色
    theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
          panel.grid = element_blank()) + 
    ylim(0, 1.5) +
    theme(axis.text.x=element_text(family="Times",size=12,face="bold"), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
          axis.text.y=element_text(family="Times",size=12,face="bold"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
          axis.title.y=element_text(family="Times",size = 15,face="bold"), #设置y轴标题的字体属性
          axis.title.x=element_text(family="Times",size = 10,face="bold"),
          plot.title = element_text(hjust = 0.5,family="Times",size = 15,face="bold"),
          legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold"))+
    labs(title="Tumor Purity", x="", y="Tumor Purity")
  p4
  p5 <- cowplot::plot_grid(p1,p2,p3,p4,
                           nrow = 2, 
                           align = 'h', 
                           vjust = -0.3)
  p5
}

# 13 CIBERSORT----------
setwd("/data/nas1/luchunlin/project/BJTC-292")
if (! dir.exists("./12_CIBERSORT")){
  dir.create("./12_CIBERSORT")
}
setwd("./12_CIBERSORT")

LM22gene<-read_xlsx('F:/luchunlin/project/BJTC-292/12_CIBERSORT//LM22gene.xlsx')%>%as.data.frame()
CIBERSORT_exp<-gsea_exp[rownames(gsea_exp)%in%LM22gene$Gene,]

setwd('F:/luchunlin/pipeline/CIBERSORT')

write.table(CIBERSORT_exp,
            file = "CIBERSORT_exp.txt",
            quote = F,
            sep = "\t",
            row.names = T)
{ 
  source("Cibersort.R")
  result <- CIBERSORT('/data/nas1/luchunlin/pipeline/CIBERSORT/LM22.txt',
                      'CIBERSORT_exp.txt', 
                      perm = 1000, ##Permutations for significance analysis是用来计算单个样本估算免疫浸润的p值，大多数文章会采用1000次。数值越大，运行时间越久，
                      QN = F)
  cibersort_raw <- read.table("CIBERSORT-Results.txt",
                              header = T,
                              sep = "\t",
                              row.names = 1,
                              check.names = F)
  cibersort_result <- t(cibersort_raw[,-c(23,24,25)])
}
{
  tiics_result <- cibersort_result
  pvalue = padj = log2FoldChange <- matrix(0, nrow(tiics_result), 1)
  for (i in 1:nrow(tiics_result)){
    pvalue[i, 1] = p.value = wilcox.test(tiics_result[i, high_sample],
                                         tiics_result[i, low_sample])$p.value
    log2FoldChange[i, 1] = mean(tiics_result[i, high_sample]) - 
      mean(tiics_result[i, low_sample])
  }
  padj <- p.adjust(as.vector(pvalue), "fdr", n = length(pvalue))
  rTable <- data.frame(log2FoldChange, 
                       pvalue, 
                       padj,
                       row.names = rownames(tiics_result))
  high <- signif(apply(tiics_result[rownames(rTable), high_sample], 
                       1,
                       mean), 4)
  low <- signif(apply(tiics_result[rownames(rTable), low_sample], 
                      1, 
                      mean), 4)
  rTable <- data.frame(high, 
                       low,
                       rTable[, c("padj", "pvalue", "log2FoldChange")])
  rTable$immune_cell <- rownames(rTable)
  rTable$sig <- ifelse(rTable$padj < 0.05,
                       ifelse(rTable$padj < 0.01, 
                              ifelse(rTable$padj < 0.001,
                                     ifelse(rTable$padj < 0.0001,
                                            paste(rTable$immune_cell, "****",  sep = ""),
                                            paste(rTable$immune_cell, "***", sep = "")),
                                     paste(rTable$immune_cell, "**", sep = "")),
                              paste(rTable$immune_cell, "*",  sep = "")), 
                       rTable$immune_cell)
  
  write.table(rTable,
              file = "cibersort_tiics_wilcox_test.xls",
              quote = F,
              row.names = F,
              sep = '\t')
}
diff_cibersort_Table<-rTable[which(rTable$pvalue<0.05),]
write.table(diff_cibersort_Table,
            file = 'diff_cibersort_Table.xls',
            quote = F,
            sep = '\t')
## 13-1 柱状堆叠图--------
cibersort_result<-as.data.frame(cibersort_result)
cibersort_result$cell<-rownames(cibersort_result)
box_dat <- gather(cibersort_result, key=sample, value='score', -c("cell"))
box_dat$group<-ifelse(box_dat$sample%in%high_sample,'High','Low')

library(RColorBrewer)
mypalette <- colorRampPalette(brewer.pal(7,"Paired"))

box_plot <- ggplot(box_dat, aes(x=sample, 
                                y=100*score,
                                fill=cell)) +
  geom_bar(position = 'stack',stat = 'identity')+
  scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  labs(x='',
       y='Relative Percent',
       fill='')+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'top') +
  scale_fill_manual(values = mypalette(22))+
  facet_grid(~box_dat$group,scales= "free",space= "free")
box_plot
## 13-2 小提琴图--------
violin_dat<-cibersort_result[rownames(cibersort_result)%in%rownames(diff_cibersort_Table),]%>%as.data.frame()
#violin_dat<-cibersort_result%>%as.data.frame()
#rTable<-rTable[order(rTable$pvalue,decreasing = F),]
#violin_dat<-violin_dat[rownames(rTable),]
violin_dat$cell<-rownames(violin_dat)
violin_dat <- gather(violin_dat, key=sample, value=score, -c("cell"))
violin_dat$group <- ifelse(violin_dat$sample%in%high_sample,'High','Low')
#library(rstatix)
#stat.test<-violin_dat%>%
#  group_by(cell)%>%
#  wilcox_test(score~group)%>%
#  adjust_pvalue(method = 'fdr')%>%
#  add_significance("p.adj")

head(violin_dat)

violin_plot <- ggplot(violin_dat, aes(x=cell,
                                      y=score,
                                      fill=group)) +
  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  #stat_boxplot(geom="errorbar", 
  #             width=0.1,
  #             position = position_dodge(0.9)) +
  #geom_boxplot(width=0.7,
  #             position=position_dodge(0.9),
  #             outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#FF6A6A", "#20B2AA"), name = "Group")+
  labs(title="Immune Cell", x="", y = "Fraction",size=20) +
  stat_compare_means(data = violin_dat,
                     mapping = aes(group = group),
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
        panel.grid.minor = element_blank())+facet_wrap(~cell,scales = "free",nrow = 2) 
violin_plot

# 14 checkpoint----------
setwd("/data/nas1/luchunlin/project/BJTC-292")
if (! dir.exists("./13_checkpoint")){
  dir.create("./13_checkpoint")
}
setwd("./13_checkpoint")
checkpoint <- read.table("checkpoint.txt",
                         header = F)
checkpoint <- checkpoint$V1
length(checkpoint)
checkpoint_DEG <- allGeneSets[which(rownames(allGeneSets)%in%checkpoint),]
dim(checkpoint_DEG)
# 40 7
logFC_cutoff <- 0
checkpoint_DEG$change = as.factor(
  ifelse(checkpoint_DEG$padj < 0.05 & abs(checkpoint_DEG$log2FoldChange) > logFC_cutoff,
         ifelse(checkpoint_DEG$log2FoldChange > logFC_cutoff ,'UP','DOWN'),'NOT')
)

sig_checkpoint <- rownames(checkpoint_DEG[which(checkpoint_DEG$padj< 0.05),])
length(sig_checkpoint)
#14
sig_checkpoint
# "C10orf54" "CD80"     "CD244"    "CD274"    "CD200R1"  "TMIGD2"   "CD86"     "TNFRSF9"  "TNFRSF8"  "CD160"    "CTLA4"    "TNFSF14"  "PDCD1"    "ICOSLG" 
write.table(sig_checkpoint, 
            'sig_checkpoint.xls',
            quote = F, sep="\t")
sig_expr <- as.data.frame(gsea_exp[sig_checkpoint,])
sig_expr <- log2(sig_expr + 1)
sig_expr$gene <- rownames(sig_expr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(Ipaper)
library(rstatix)
head(sig_expr[,1:3])
violin_dat <- gather(sig_expr, key=sample, value='log2(expr+1)', -c("gene"))
head(violin_dat)
violin_dat$group <- ifelse(violin_dat$sample %in% high_sample,
                           "High", "Low") 
head(violin_dat)

violin_plot <- ggplot(violin_dat, aes(x=gene, 
                                      y=`log2(expr+1)`,
                                      fill=group)) +
  #  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#FF6A6A", "#20B2AA"), name = "Group")+
  labs(title="Immune Checkpoint", x="", y = "log2(expr+1)",size=20) +
  stat_compare_means(data = violin_dat,
                     mapping = aes(group = group),
                     label ="p.signif") +
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
        panel.grid.minor = element_blank())
violin_plot
# 15 化疗药物敏感性---------
setwd("/data/nas1/luchunlin/project/BJTC-292")
if (! dir.exists("./14_Medicinal_Sensity")){
  dir.create("./14_Medicinal_Sensity")
}
setwd("./14_Medicinal_Sensity")

#install.packages("pRRophetic_0.5.tar.gz", repos = NULL, dependencies = TRUE)
library(pRRophetic)
library(ggplot2)
set.seed(12345)
model_expr<-dataFilt.AML.final[lasso_geneids$lasso_geneids, risk$id]
riskscore<-data.frame(risk$id,risk$risk)
colnames(riskscore)<-c('sample','risk')
riskscore$risk[which(riskscore$risk==1)] <-'Low risk'
riskscore$risk[which(riskscore$risk==0)] <-'High risk'
head(riskscore)
drug<-read.table(file = 'drugs.txt',sep='\t',header=F)
ic50<-data.frame(riskscore$sample)

a<-data.frame(row.names=riskscore$sample,riskscore$risk)

colnames(a)<-'risk'

cnt<-1

while (cnt < 139) {
  
  predictedPtype <- pRRopheticPredict(as.matrix(model_expr), drug[cnt,],selection=1)
  
  Tipifarnib<-data.frame(predictedPtype)
  
  colnames(Tipifarnib)<-drug[cnt,]
  
  a<-cbind(a,Tipifarnib)
  
  cnt = cnt + 1
}

write.table(a,'IC50.xls',sep='\t',quote=F)

b<-a
b[b<0]<-NA
# 先写成函数的形式，方便调用
removeRowsAllNa  <- function(x){x[apply(x, 1, function(y) any(!is.na(y))),]}
removeColsAllNa  <- function(x){x[, apply(x, 2, function(y) any(!is.na(y)))]}
c<-removeColsAllNa(b)
na_flag <- apply(is.na(c), 2, sum)
x <- c[, which(na_flag == 0)]
View(x)
dim(x)
# [1] 132  98
medicinal_result <- t(subset(x, select = -risk)) 
high_group <- risk$id[which(risk$risk==0)]
low_group <- risk$id[which(risk$risk==1)]
pvalue = padj = log2FoldChange <- matrix(0, nrow(medicinal_result), 1)
for (i in 1:nrow(medicinal_result)){
  pvalue[i, 1] = p.value = wilcox.test(medicinal_result[i, high_group],
                                       medicinal_result[i, low_group])$p.value
  log2FoldChange[i, 1] = mean(medicinal_result[i, high_group]) - 
    mean(medicinal_result[i, low_group])
}
padj <- p.adjust(as.vector(pvalue), "fdr", n = length(pvalue))
rTable <- data.frame(log2FoldChange, 
                     pvalue, 
                     padj,
                     row.names = rownames(medicinal_result))
high_group_res <- signif(apply(medicinal_result[rownames(rTable), high_group], 
                               1,
                               median), 4)
low_group_res <- signif(apply(medicinal_result[rownames(rTable), low_group], 
                              1, 
                              median), 4)
rTable <- data.frame(high_group_res, 
                     low_group_res,
                     rTable[, c("padj", "pvalue", "log2FoldChange")])
rTable$drugs <- rownames(rTable)
rTable$sig <- ifelse(rTable$padj < 0.05,
                     ifelse(rTable$padj < 0.01, 
                            ifelse(rTable$padj < 0.001,
                                   ifelse(rTable$padj < 0.0001,
                                          paste(rTable$drugs, "****",  sep = ""),
                                          paste(rTable$drugs, "***", sep = "")),
                                   paste(rTable$drugs, "**", sep = "")),
                            paste(rTable$drugs, "*",  sep = "")), 
                     rTable$drugs)

write.table(rTable,
            file = "drugs_wilcox_test.xls",
            quote = F,
            row.names = F)

### 发散条形图绘制

#install.packages('ggprism')
library(ggprism)
## 横坐标药物，纵坐标：IC50(H)/IC50(L)-1
dat_plot<-data.frame(drug=rownames(rTable),
                     'IC50(H)/IC50(L)-1'=(rTable$high_group_res/rTable$low_group_res-1),
                     pvalue=rTable$pvalue,
                     padj=rTable$padj)
dat_plot$threshold=factor(ifelse(dat_plot$padj<0.05&dat_plot$pvalue<0.05,'P<0.05 & FDR<0.05',ifelse(dat_plot$pvalue<0.05&dat_plot$padj>0.05,'P<0.05 & FDR>0.05','P>0.05 & FDR>0.05')))
dat_plot<-dat_plot%>%arrange(desc(dat_plot$IC50.H..IC50.L..1))

dat_plot$drug<-factor(dat_plot$drug,levels=dat_plot$drug)

p <- ggplot(data = dat_plot,aes(x = drug,y = IC50.H..IC50.L..1,fill = threshold)) +
  geom_col()+
  scale_fill_manual(values = c('P<0.05 & FDR<0.05'= '#5F9EA0','P>0.05 & FDR>0.05'='#cccccc','P<0.05 & FDR>0.05'='#FFD700')) +
  xlab('') + 
  ylab('Exp(Median IC50(H))/Exp(Median IC50(L))-1') + 
  theme_prism(border = T) +
  theme(
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 65,size = 7,
                               hjust = 1,vjust = 1),
    axis.text.y = element_text(size = 13),
    legend.position = c(0.85,0.85),
    legend.text = element_text(size = 8,face = 'bold')
  )

# 16 免疫治疗反应预测---------
setwd("/data/nas1/luchunlin/project/BJTC-292")
if (! dir.exists("./15_immunotherapy")){
  dir.create("./15_immunotherapy")
}
setwd("./15_immunotherapy")
## 16-1 TIDE-------
tide_dat <- dataFilt.AML.final
# fpkm转TPM
FPKM2TPM <- function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}

# tide_dat[which(tide_dat<0)] <- 0
tide_dat <- apply(tide_dat,2,FPKM2TPM)
tide_dat <- log2(tide_dat + 1)
rownmean <- apply(tide_dat,1,mean)
tide_dat2 <- sweep(tide_dat, 1, rownmean)
dim(tide_dat2)
write.table(tide_dat2,
            file ="tide_dat.txt",
            sep = "\t",
            quote = F,
            row.names = T)

tide_result <- read.csv("AML_tide_result.csv",header = T)
View(tide_result)
tide_result2 <- subset(tide_result, select = c("Patient", "TIDE"))
tide_plot_dat <- data.frame(Patient=risk$id,
                            riskScore=risk$riskScore,
                            risk_group=risk$risk)
tide_plot_dat <- merge(tide_plot_dat, tide_result2, by = "Patient")
tide_plot_dat$risk_group <- factor(tide_plot_dat$risk_group,
                                   level = c(0, 1),
                                   labels = c("High Risk", "Low Risk"))

dim(tide_plot_dat)
library(ggplot2)
library(ggpubr)
#install.packages('ggside')
#install.packages('ggstatsplot')
{
  library(ggstatsplot)
  tide_cor <- ggscatterstats(data = tide_plot_dat,
                             y = riskScore,
                             x = TIDE,
                             centrality.para = "mean",
                             margins = "both",
                             xfill = "#A73030FF",
                             yfill = "#0073C2FF",
                             type = "pearson",
                             xlab = "TIDE Prediction Score",
                             marginal.type = "histogram",
                             title = "Relationship between TIDE Prediction Score and riskScore"
  )+theme(axis.title.x = element_text(size = 15),
          axis.title.y = element_text(size = 15),
          plot.title = element_text(size = 12))
  tide_cor
  ggsave(filename = "tide_riskscore_cor.png", height = 4, width = 6,tide_cor)
  ggsave(filename = "tide_riskscore_cor.pdf", height = 4, width = 6,tide_cor)
  
  tide_box <- ggboxplot(tide_plot_dat,
                        x = "risk_group",
                        y = "TIDE",
                        fill = "risk_group",
                        palette =c("#FF6A6A", "#20B2AA")) +
    stat_compare_means(label.y = 3.2) +
    theme_bw() +
    labs(title = "", x = "", y = "TIDE Prediction Score") +
    theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18),
          axis.text.x=element_text(colour="black",face="bold",size=15), 
          axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12),
          axis.title.y=element_text(size=18,face="bold"),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  tide_box
  ggsave(filename = "tide_riskscore_boxplot.png", height = 5, width = 4,tide_box)
  ggsave(filename = "tide_riskscore_boxplot.pdf", height = 5, width = 4,tide_box)
}
## 16-2 IPS-----
## IOBR包

#if (!requireNamespace("IOBR", quietly = TRUE))
#  devtools::install_github("IOBR/IOBR")
library(IOBR)
ips_dat<-tide_dat
ips<-deconvo_tme(eset = ips_dat, method = "ips", plot= FALSE)
# ips<-IPS_calculation(eset = ips_dat,plot = F)
head(ips)
ips_result<-subset(ips,select=c('ID','IPS_IPS'))
colnames(ips_result)<-c('ID','IPS')
ips_plot_dat<-data.frame(ID=risk$id,
                         riskScore=risk$riskScore,
                         risk_group=risk$risk)
ips_plot_dat<-merge(ips_plot_dat,ips_result,by='ID')
ips_plot_dat$risk_group<-factor(ips_plot_dat$risk_group,
                                levels = c(0,1),
                                labels = c("High Risk", "Low Risk"))
# BiocManager::install('ggridges')
my_comparisons<-list(c('High Risk','Low Risk'))
ggplot(ips_plot_dat, aes(x = IPS, y = risk_group)) +
  geom_density_ridges_gradient(aes(fill = risk_group),
                               scale = 1, size = 0.3,
                               position = 'identity') +
  geom_signif(comparisons = my_comparisons,
              map_signif_level = F,
              y_position = c(11))+
  theme(legend.position = "none")+
  labs(y = "")+
  theme(axis.text.y=element_text(hjust=0.5,colour="black",size=12),
        axis.text.x=element_text(hjust=0.5,colour="black",size=12))+
  scale_fill_manual(values= c("#FF6A6A", "#20B2AA"))

# 17 预后模型的构建与评价--------
setwd("/data/nas1/luchunlin/project/BJTC-292")
if (! dir.exists("./16_prog_model")){
  dir.create("./16_prog_model")
}
setwd("./16_prog_model")
## 17-1 单因素cox--------
train_phenotype3<-data.frame(sample=phenotype$submitter_id.samples,
                            age=phenotype$age_at_initial_pathologic_diagnosis,
                            gender=phenotype$gender.demographic,
                            FAB=phenotype$leukemia_french_american_british_morphology_code)
train_phenotype3<-merge(train_phenotype3,survival,by='sample')
train_phenotype3<-train_phenotype3[,-6]
train_phenotype3$gender<-ifelse(train_phenotype3$gender=='male',1,2)
train_phenotype3$FAB<-gsub(' Undifferentiated','',train_phenotype3$FAB,fixed = T)
train_phenotype3$FAB<-gsub('Not Classified',NA,train_phenotype3$FAB)
train_phenotype3$FAB<-gsub('M','',train_phenotype3$FAB,fixed = T)
colnames(train_phenotype3)<-c('id','age','gender','FAB','OS','OS.time')

train_risk_clinical <- merge(train_phenotype3,
                             sub_risk,
                             by = "id")
rownames(train_risk_clinical) <- train_risk_clinical$id
train_risk_clinical = subset(train_risk_clinical, select = -c(id))
dim(train_risk_clinical)
colnames_train <- colnames(train_risk_clinical)
covariates_train <- colnames_train[-which(colnames_train %in% c("OS", "OS.time"))]

train_risk_clinical$FAB<-factor(train_risk_clinical$FAB)
library(survival)
res.risk = coxph(Surv(time = OS.time, event = OS) ~ riskScore, data = train_risk_clinical) %>% summary
res.risk = c(res.risk$conf.int[-2], res.risk$coefficients[5])
res.age = coxph(Surv(time = OS.time, event = OS) ~ age, data = train_risk_clinical) %>% summary
res.age = c(res.age$conf.int[-2], res.age$coefficients[5])
res.gender = coxph(Surv(time = OS.time, event = OS) ~ gender, data = train_risk_clinical) %>% summary
res.gender = c(res.gender$conf.int[-2], res.gender$coefficients[5])
res.FAB = coxph(Surv(time = OS.time, event = OS) ~ FAB, data = train_risk_clinical) %>% summary
res.FAB = cbind(res.FAB$conf.int[,-2], res.FAB$coefficients[,5])

res.ref = c(1,1,1,NA)
res = rbind(res.risk, res.age, res.gender,res.ref,res.FAB) %>% as.data.frame()
rownames(res)
res$Indicators = c("riskScore","Age","Gender","M0(Reference)","M1","M2","M3","M4","M5","M6","M7")
colnames(res) = c("hr","low","up","pv","Indicator")
res$p = signif(res$pv, 2) %>% paste0("p = ", .)
res$p[is.na(res$pv)] = NA
res$Indicator = factor(res$Indicator, levels = rev(res$Indicator))
rownames(res) <- res$Indicator
res2 <- data.frame(p.value=res$pv,
                   HR=res$hr,
                   HR.95L=res$low,
                   HR.95H=res$up,
                   Indicator=res$Indicator)
rownames(res2) <- res2$Indicator
write.table(res2, file = "univariate_cox_prog_forest.xls", sep = "\t", quote = F)
res2 <- subset(res2, select = -c(Indicator))
library(tidyr)
hz <- paste(round(res2$HR,3),
            "(",round(res2$HR.95L,3),
            "-",round(res2$HR.95H,3),")",sep = "")
hz
hz[4] <- ""

tabletext <- cbind(c(NA,rownames(res2)),
                   c("P value",ifelse(res2$p.value<0.0001,
                                      "< 0.0001",
                                      round(res2$p.value,4))),
                   c("Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1
library(forestplot)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,rep(FALSE, 7)),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,res2$HR),
           lower=c(NA,res2$HR.95L), #95%置信区间下限
           upper=c(NA,res2$HR.95H), #95%置信区间上限
           boxsize=0.1,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(0,5,10,15,20,25,30), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(1.2,"cm"), #固定行高
           graphwidth = unit(.6,"npc"), #图在表中的宽度比例
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),
           # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
           #                 "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
           #                 "59" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
           # mar=unit(rep(0.01, times = 4), "cm"),#图形页边距
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Hazard Ratio",
           grid = T,
           title = "Univariate") # 垂直于x轴的网格线，对应每个刻度

## 17-2 多因素cox--------
res.mul = coxph(Surv(time = OS.time, event = OS) ~ riskScore + age, data = train_risk_clinical)%>% summary
res.mul = cbind(res.mul$conf.int[,-2], res.mul$coefficients[,5]) %>% as.data.frame()

res.mul$Indicators = c("riskScore","Age")
colnames(res.mul) = c("hr","low","up","pv","Indicator")
res.mul$p = signif(res.mul$pv, 2) %>% paste0("p = ", .)
res.mul$p[is.na(res.mul$pv)] = NA
res.mul$Indicator = factor(res.mul$Indicator, levels = rev(res.mul$Indicator))
rownames(res.mul) <- res.mul$Indicator

multi_res <- data.frame(p.value=res.mul$pv,
                        HR=res.mul$hr,
                        HR.95L=res.mul$low,
                        HR.95H=res.mul$up,
                        Indicator=res.mul$Indicator)
rownames(multi_res) <- multi_res$Indicator
multi_res
write.csv(multi_res,
          file = "multivariate_cox_prog_result.csv",
          quote = F,
          row.names = T)
multi_res <- subset(multi_res, select = -c(Indicator))

library(tidyr)
hz <- paste(round(multi_res$HR,3),
            "(",round(multi_res$HR.95L,3),
            "-",round(multi_res$HR.95H,3),")",sep = "")
hz

tabletext <- cbind(c(NA,rownames(multi_res)),
                   c("P value",ifelse(multi_res$p.value<0.0001,
                                      "< 0.0001",
                                      round(multi_res$p.value,4))),
                   c("Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,multi_res$HR),
           lower=c(NA,multi_res$HR.95L), #95%置信区间下限
           upper=c(NA,multi_res$HR.95H), #95%置信区间上限
           boxsize=0.1,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(0,1,2,3,4,5,6,7,8), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(1.2,"cm"), #固定行高
           graphwidth = unit(.6,"npc"), #图在表中的宽度比例
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),
           # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
           #                 "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
           #                 "59" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
           # mar=unit(rep(0.01, times = 4), "cm"),#图形页边距
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Hazard Ratio",
           grid = T,
           title = "Multivariate") # 垂直于x轴的网格线，对应每个刻度
## 17-3 构建COX模型，绘制列线图------
multi_cov<-c('riskScore','age')
cox_data_prog <- as.formula(paste0('Surv(OS.time, OS)~',
                                   paste(multi_cov,
                                         sep = '',
                                         collapse = '+')))

# Nomogram
library(rms)
ddist <- datadist(train_risk_clinical)
options(datadist='ddist')

# 构建COX模型，绘制列线图

res.cox <- psm(cox_data_prog,
               data = train_risk_clinical, dist = 'lognormal')
surv <- Survival(res.cox) # 构建生存概率函数
function(x) surv(365, x) # 1年事件发生概率
function(x) surv(1095, x) # 3年事件发生概率
function(x) surv(1825, x) # 5年事件发生概率

nom.cox <- nomogram(res.cox,
                    fun = list(function(x) surv(365, x),
                               function(x) surv(1095, x),
                               function(x) surv(1825, x)),
                    funlabel=c("1-year Survival Probability", "3-year Survival Probability", "5-year Survival Probability"),
                    maxscale = 10,
                    fun.at = c(0.01,seq(0.1,0.9,by=0.2),0.95,0.99),
                    lp=F)


plot(nom.cox, cex.axis  = 1.5, cex.var = 1.7)
png(filename = "nomogram_line_points.png", height = 600, width = 1200)
plot(nom.cox, cex.axis  = 1.5, cex.var = 1.7)
dev.off()
pdf(file = "nomogram_line_points.pdf", height = 8, width = 17)
plot(nom.cox, cex.axis  = 1.5, cex.var = 1.7)
dev.off()

## 17-4 构建校准曲线--------
##绘制1年生存期校准曲线
coxm_1 <- cph(cox_data_prog,
              data=train_risk_clinical,
              surv=T,
              x=T,
              y=T,
              time.inc = 365)
cal_1<-calibrate(coxm_1,u=365,cmethod='KM',m=100,B=1000)

##绘制3年生存期校曲线
##time.in 和 u 要是一样的，都是要评价的时间节点
coxm_3 <- cph(cox_data_prog,
              data=train_risk_clinical,
              surv=T,
              x=T,y=T,
              time.inc = 3*365)
cal_3 <-calibrate(coxm_3,u=3*365,cmethod='KM',m=100,B=1000)

coxm_5 <- cph(cox_data_prog,
              data=train_risk_clinical,
              surv=T,
              x=T,y=T,
              time.inc = 5*365)
cal_5 <-calibrate(coxm_5,u=5*365,cmethod='KM',m=100)

par(mar=c(7,4,4,3),cex=1.5)
plot(cal_1,
     subtitles = F,
     lwd=2,lty=1, ##设置线条形状和尺寸
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 1-5 year Progression-free Interval',#便签
     ylab='Actual 1-5 year Progression-free Interval(Proportion)',#标签
     col="#00468b",#设置一个颜色
     xlim = c(0,1),ylim = c(0,1)) ##x轴和y轴范围
plot(cal_3,
     add = T,
     subtitles = F,
     lwd=2,lty=1,  ##设置线条宽度和线条类型
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 1-5 year Progression-free Interval',#便签
     ylab='Actual 1-5 year Progression-free Interval(Proportion)',#标签
     col="#ed0000",#设置一个颜色
     xlim = c(0,1),ylim = c(0,1)) ##x轴和y轴范围
plot(cal_5,
     add = T,
     subtitles = F,
     lwd=2,lty=1, ##设置线条形状和尺寸
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 1-5year Progression-free Interval',#便签
     ylab='Actual 1-year Progression-free Interval(Proportion)',#标签
     col="#42b540",#设置一个颜色
     xlim = c(0,1),ylim = c(0,1)) ##x轴和y轴范围



#加上图例
legend("bottomright", legend=c("1-year", "3-year", "5-year"), 
       col=c("#00468b", "#ed0000", "#42b540"), 
       lwd=2)
#调整对角线
abline(0,1,lty=5,lwd=2,col="grey")
## 17-5 KM曲线--------

res.mul = coxph(Surv(time = OS.time, event = OS)~ riskScore + age, data = train_risk_clinical)
train_risk_clinical$riskScore_prog = predict(res.mul, newdata = train_risk_clinical, type = "lp")
train_risk_clinical$risk_prog = ifelse(train_risk_clinical$riskScore_prog > median(train_risk_clinical$riskScore_prog, na.rm = T), "high", "low")
train_risk_clinical$risk_prog = factor(train_risk_clinical$risk_prog, levels = c("high", "low"), labels = c("High risk", "Low risk"))
surv.fit = survfit(Surv(time = OS.time,event = OS) ~ risk_prog, data = train_risk_clinical)
prog_survival_median <- ggsurvplot(surv.fit,
                                   pval = TRUE, 
                                   conf.int = F,
                                   legend.labs=c("High risk","Low risk" ),
                                   legend.title="Risk score",
                                   title="Prognosis Model KM",
                                   font.main = c(15,"bold"),
                                   risk.table = TRUE, 
                                   risk.table.col = "strata", 
                                   linetype = "strata", 
                                   surv.median.line = "hv",
                                   ggtheme = theme_bw(), 
                                   palette = c("#A73030FF", "#0073C2FF"))
prog_survival_median
## 17-6 ROC曲线-------
# 开始验证
train_risk_clinical2 <- train_risk_clinical
library(survival)
library(survminer)

multi_ROC_out <- function(time_vector, risk_score_table){
  library(survivalROC)
  single_ROC <- function(single_time){
    for_ROC <- survivalROC(Stime=risk_score_table$OS.time,
                           status=risk_score_table$OS,
                           marker=risk_score_table$riskScore_prog,
                           predict.time=single_time,method = 'KM')
    data.frame('True_positive'=for_ROC$TP, 'False_positive'=for_ROC$FP,
               'Cut_values'=for_ROC$cut.values, 'Time'=rep(single_time, length(for_ROC$TP)),
               'AUC'=rep(for_ROC$AUC, length(for_ROC$TP)))
  }
  multi_ROC_list <- lapply(time_vector, single_ROC)
  do.call(rbind, multi_ROC_list)
}

for_multi_ROC <- multi_ROC_out(time_vector = c(365*seq(1,5,2)), 
                               risk_score_table = train_risk_clinical2)
for_multi_ROC$Time <- factor(for_multi_ROC$Time)

# 画ROC曲线 
library(scales)
library(geomROC)
library(plotROC)
auc_y1 <- round(for_multi_ROC[which(for_multi_ROC$Time==365),5][1],2)
auc_y3 <- round(for_multi_ROC[which(for_multi_ROC$Time==1095),5][1],2)
auc_y5 <- round(for_multi_ROC[which(for_multi_ROC$Time==1825),5][1],2)

ROC <- ggplot(for_multi_ROC, aes(x=False_positive,
                                 y=True_positive, 
                                 label=Cut_values, 
                                 color=Time)) + 
  scale_color_manual(breaks = c("365", "1095", "1825"),
                     labels = c("1 years", "3 years", "5 years"),
                     values = c("#4682B4", "#FF4040", "#20B2AA")) +
  geom_roc(labels = F, stat = 'identity') + 
  style_roc() + 
  geom_abline(slope = 1, intercept = 0, color = 'gray', linetype=2) +
  theme_bw() +
  labs(title = "ROC for Prognosis Model") +
  theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
        panel.grid = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold")) +
  # annotate('text', x=.75, y=.25, label=paste('AUC of 1 years =', round(auc_y1,2))) + 
  # annotate('text', x=.75, y=.15, label=paste('AUC of 2 years =', round(auc_y3,2))) + 
  # annotate('text', x=.75, y=.05, label=paste('AUC of 3 years =', round(auc_y5,2))) +
  annotate("text", x=0.75, y=c(0.25, 0.15, 0.05),
           label = c(paste('AUC of 1 years =', format(auc_y1,nsmall=2)),
                     paste('AUC of 3 years =', format(auc_y3,nsmall=2)),
                     paste('AUC of 5 years =', format(auc_y5,nsmall=2))))
ROC
ggsave('ROC for Prognosis Model.png', ROC,width = 5, height = 4)
ggsave('ROC for Prognosis Model.pdf', ROC,width = 5, height = 4)
