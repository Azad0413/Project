rm(list = ls())
setwd("/data/nas1/luchunlin/project/HZ0301-3/")
if (! dir.exists("./11_GSVA")){
  dir.create("./11_GSVA")
}
setwd("./11_GSVA")

dat<-read.csv('../00_rawdata/dat.fpkm.xls',sep = '\t',row.names = 1)
colnames(dat)<-gsub('.','-',colnames(dat),fixed = T)
library(GSVA)
library(GSEABase)
library(limma)
group<-read.delim2('../07_risk/risk.xls')
high.sample<-group$id[which(group$risk==0)]
group2 <- group%>%dplyr::select(c('id','risk'))
group2$risk<-ifelse(group2$risk==0,'High','Low')
colnames(group2)
colnames(group2)<-c('id','label')
gsva_exp<-log2(dat[,group2$id]+1)
all(colnames(gsva_exp) == group2$id)
dim(gsva_exp)
# 分组
group_score <- group2$label %>% as.factor()
design_score <- model.matrix(~0 + group_score)
rownames(design_score) <- colnames(gsva_exp)
colnames(design_score) <- levels(group_score)
compare_score <- makeContrasts("High-Low", levels = design_score)
KEGG_ref <- getGmt("/data/nas1/luchunlin/pipeline/GSVA/c2.cp.kegg.v7.4.symbols.gmt")
es_KEGG <- gsva(as.matrix(gsva_exp), KEGG_ref,
                min.sz=10, max.sz=1000, verbose=TRUE)
es_KEGG <- as.data.frame(es_KEGG)

fit <- lmFit(es_KEGG, design_score)
fit2 <- contrasts.fit(fit ,compare_score)
fit3 <- eBayes(fit2)
allGeneSets <- topTable(fit3, coef = 1, number = Inf)
logFCcutoff <- 0
allGeneSets$change = as.factor(
  ifelse(allGeneSets$P.Value < 0.05 & abs(allGeneSets$logFC) > logFCcutoff,
         ifelse(allGeneSets$logFC > logFCcutoff,'UP','DOWN'),'NOT')
)

DEGeneSets <- subset(allGeneSets,
                     allGeneSets$P.Value < 0.05 & abs(allGeneSets$t) > logFCcutoff )

write.table(allGeneSets, 
            file = "KEGG.xls",
            quote = F,
            sep = "\t",
            row.names = T)
DEGeneSets <- DEGeneSets[order(DEGeneSets$adj.P.Val),]
dim(DEGeneSets)
# [1]146  6
write.table(DEGeneSets,
            file = 'diff_KEGG.xls',
            quote = F,
            sep = "\t",
            row.names = T)
##热图
group.gsva<-group2[order(group2$label),]

## 按照分组排一下序
heat<-es_KEGG[rownames(es_KEGG)%in%
                rownames(rbind(head(DEGeneSets[order(DEGeneSets$logFC,decreasing = T),],10),
                               head(DEGeneSets[order(DEGeneSets$logFC,decreasing = F),],10))),]
colnames(group.gsva)<-c('sample','group')
heat <- heat[,group.gsva$sample]
annotation_col<-as.data.frame(group.gsva$group)
colnames(annotation_col)='Group'
rownames(annotation_col)=colnames(heat)
color.key<-c("#3300CC", "#3399FF", "white", "#FF3333", "#CC0000")
ann_colors<-list(
  Group = c(Low="#00FFFF",High="#FFAEB9"))
library(pheatmap)
pdf(file = '01.heatmap(KEGG).pdf',w=12,h=6)
pheatmap(heat,
         color = colorRampPalette(color.key)(50),
         border_color = NA,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         labels_row = NULL,
         clustering_method = 'ward.D2',
         show_rownames = T,
         show_colnames = F,
         fontsize_col = 5,
         cluster_cols = F,
         cluster_rows = T,
         width = 10,
         main = 'GSVA KEGG')
dev.off()

png(file = '01.heatmap(KEGG).png',w=900,h=500)
pheatmap(heat,
         color = colorRampPalette(color.key)(50),
         border_color = NA,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         labels_row = NULL,
         clustering_method = 'ward.D2',
         show_rownames = T,
         show_colnames = F,
         fontsize_col = 5,
         cluster_cols = F,
         cluster_rows = T,
         width = 10,
         main = 'GSVA KEGG')
dev.off()

GO_ref <- getGmt("/data/nas1/luchunlin/pipeline/GSVA/c5.go.bp.v7.4.symbols.gmt")
es_GO <- gsva(as.matrix(gsva_exp), GO_ref,
                min.sz=10, max.sz=1000, verbose=TRUE)
es_GO <- as.data.frame(es_GO)

fit <- lmFit(es_GO, design_score)
fit2 <- contrasts.fit(fit ,compare_score)
fit3 <- eBayes(fit2)
allGeneSets <- topTable(fit3, coef = 1, number = Inf)
logFCcutoff <- 0
allGeneSets$change = as.factor(
  ifelse(allGeneSets$P.Value < 0.05 & abs(allGeneSets$logFC) > logFCcutoff,
         ifelse(allGeneSets$logFC > logFCcutoff,'UP','DOWN'),'NOT')
)

DEGeneSets <- subset(allGeneSets,
                     allGeneSets$P.Value < 0.05 & abs(allGeneSets$t) > logFCcutoff )

write.table(allGeneSets, 
            file = "GO.xls",
            quote = F,
            sep = "\t",
            row.names = T)
DEGeneSets <- DEGeneSets[order(DEGeneSets$adj.P.Val),]
dim(DEGeneSets)
# [1]146  6
write.table(DEGeneSets,
            file = 'diff_GO.xls',
            quote = F,
            sep = "\t",
            row.names = T)
##热图
group.gsva<-group2[order(group2$label),]

## 按照分组排一下序
heat<-es_GO[rownames(es_GO)%in%
                rownames(rbind(head(DEGeneSets[order(DEGeneSets$logFC,decreasing = T),],10),
                               head(DEGeneSets[order(DEGeneSets$logFC,decreasing = F),],10))),]
colnames(group.gsva)<-c('sample','group')
heat <- heat[,group.gsva$sample]
annotation_col<-as.data.frame(group.gsva$group)
colnames(annotation_col)='Group'
rownames(annotation_col)=colnames(heat)
color.key<-c("#3300CC", "#3399FF", "white", "#FF3333", "#CC0000")
ann_colors<-list(
  Group = c(Low="#00FFFF",High="#FFAEB9"))
library(pheatmap)
pdf(file = '02.heatmap(GO).pdf',w=12,h=6)
pheatmap(heat,
         color = colorRampPalette(color.key)(50),
         border_color = NA,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         labels_row = NULL,
         clustering_method = 'ward.D2',
         show_rownames = T,
         show_colnames = F,
         fontsize_col = 5,
         cluster_cols = F,
         cluster_rows = T,
         width = 10,
         main = 'GSVA GO')
dev.off()

png(file = '02.heatmap(GO).png',w=900,h=500)
pheatmap(heat,
         color = colorRampPalette(color.key)(50),
         border_color = NA,
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         labels_row = NULL,
         clustering_method = 'ward.D2',
         show_rownames = T,
         show_colnames = F,
         fontsize_col = 5,
         cluster_cols = F,
         cluster_rows = T,
         width = 10,
         main = 'GSVA GO')
dev.off()

