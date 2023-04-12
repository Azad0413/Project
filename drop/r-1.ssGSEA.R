rm(list = ls())
setwd("/data/nas1/luchunlin/project/JNZK-214-8/")
if (! dir.exists("./01_ssGSEA")){
  dir.create("./01_ssGSEA")
}
setwd("./01_ssGSEA")
library(lance)
dat<-read.delim2("/data/nas1/luchunlin/project/JNZK-214-8/00_rawdata/dat.fpkm.xls", row.names = 1)%>% lc.tableToNum
colnames(dat)<-gsub('.','-',colnames(dat),fixed = T)
dat<-log2(dat+1)
group<-data.frame(sample=colnames(dat),group=c(rep('control',41),rep('Tumor',448)))
write.table(group,file = 'group.xls',sep = '\t',row.names = F,quote = F)
# tumor.sample<-group$sample[which(group$group=='Tumor')]
# dat<-dat[,tumor.sample]

#setwd('../00_rawdata/')
#write.table(dat,file = 'tumor.fpkm.xls',sep = '\t',row.names = T,quote = F)
#setwd('../01_ssGSEA/')
##ssGSEA--------
library(GSVA)
gene_set <- read.table("/data/nas1/luchunlin/pipeline/ssGSEA/mmc3.txt",
                       header = T,
                       sep ="\t")
dat.final <- as.matrix(dat)
gene_list <- split(as.matrix(gene_set)[,1],
                   gene_set[,2])
ssgsea_score = gsva(dat.final, gene_list, 
                    method = "ssgsea", 
                    ssgsea.norm = TRUE, 
                    verbose = TRUE)
write.table(ssgsea_score,
            file = "ssgsea_result.xls",row.names = T,
            sep = "\t",
            quote = F)
## 富集分数画热图

group<-group[order(group$group),]
ssgsea_score<-ssgsea_score[,group$sample]
annotation_col<-as.data.frame(group$group)
colnames(annotation_col)='Group'
rownames(annotation_col)=colnames(ssgsea_score)
library(pheatmap)
color.key<-c("#3300CC", "#3399FF", "white", "#FF3333", "#CC0000")
ann_colors<-list(
  Group = c(control="#00FFFF",Tumor="#FFAEB9"))
pdf('01.heatmap.pdf',w=9,h=7)
pheatmap(
  ssgsea_score,
  color = colorRampPalette(color.key)(50),
  border_color = 'darkgrey',
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  labels_row = NULL,
  clustering_method = 'ward.D2',
  show_rownames = T,
  show_colnames = F,
  fontsize_col = 5,
  cluster_cols = F,
  cluster_rows = T)
dev.off()
png('01.heatmap.png',w=900,h=700)
pheatmap(
  ssgsea_score,
  color = colorRampPalette(color.key)(50),
  border_color = 'darkgrey',
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  labels_row = NULL,
  clustering_method = 'ward.D2',
  show_rownames = T,
  show_colnames = F,
  fontsize_col = 5,
  cluster_cols = F,
  cluster_rows = T)
dev.off()

## 03-2 差异免疫细胞鉴定-------
tiics_result <- ssgsea_score
pvalue = padj = log2FoldChange <- matrix(0, nrow(tiics_result), 1)
group_High<-group[group$group=='Tumor',]
group_High<-as.character(group_High$sample)
group_Low<-group[group$group=='control',]
group_Low<-as.character(group_Low$sample)

for (i in 1:nrow(tiics_result)){
  pvalue[i, 1] = p.value = wilcox.test(tiics_result[i, group_High],
                                       tiics_result[i, group_Low])$p.value
  log2FoldChange[i, 1] = mean(tiics_result[i, group_High]) - 
    mean(tiics_result[i, group_Low])
}
padj <- p.adjust(as.vector(pvalue), "fdr", n = length(pvalue))
rTable <- data.frame(log2FoldChange, 
                     pvalue, 
                     padj,
                     row.names = rownames(tiics_result))
Low <- signif(apply(tiics_result[rownames(rTable), group_Low], 
                    1,
                    mean), 4)
High <- signif(apply(tiics_result[rownames(rTable), group_High], 
                     1, 
                     mean), 4)
rTable <- data.frame(Low, 
                     High,
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

diff_Table<-rTable[which(rTable$padj<0.05),]
## 26
write.table(rTable,
            file = "tiics_wilcox_test.xls",
            quote = F,
            row.names = F,
            sep = '\t')
write.table(diff_Table,
            file = "diff_tiics_wilcox_test.xls",
            quote = F,
            row.names = F,
            sep = '\t')
### 箱线图----
library(tidyr)
library(ggplot2)
library(ggpubr)
library(Ipaper)

xCell2 <- data.frame(Immune_Cell=rownames(tiics_result), 
                     tiics_result, 
                     pvalue=rTable$padj)
xCell3 <- xCell2[which(xCell2$pvalue<0.05),]
#xCell3<-xCell2
diff_tiics <- rownames(xCell3)
violin_dat <- gather(xCell3, key=Group, value=score, -c("Immune_Cell","pvalue"))

violin_dat$Group <- ifelse(gsub("\\.","-",violin_dat$Group) %in% group_Low,
                           "Low", "High") 
head(violin_dat)
boxplot_diff_TIICs <- ggplot(violin_dat, aes(x=Immune_Cell, 
                                             y=score,
                                             fill=Group)) +
  # geom_violin(trim=T,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)#"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar",
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  geom_point(aes(fill = Group),
             size = 0.05,
             position = position_dodge(0.9))+
  scale_fill_manual(values= c("#FF7256","#48D1CC"))+ #设置填充的颜色
  labs(title="", x="", y = "Score",size=20) +
  stat_compare_means(data = violin_dat,
                     mapping = aes(group = Group),
                     label ="p.signif",
                     hide.ns = F) +
  theme_bw()+#把背景设置为白底
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18), # 将图表标题居中
        axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=10), #设置x轴刻度标签的字体显示倾斜角度为45度，并向下调整1(hjust = 1)，字体大小为14
        axis.text.y=element_text(hjust=0.5,colour="black",size=12), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
        axis.title.x=element_text(size=16,face="bold"),#设置x轴标题的字体属性
        axis.title.y=element_text(size=14,face="bold"), #设置y轴标题的字体属性
        legend.text=element_text(face="bold", hjust = 0.5,colour="black", size=11), #设置图例的子标题的字体属性
        legend.title=element_text(face="bold", colour="black", size=11),#设置图例的总标题的字体属性
        #legend.justification=c(-0.1,1.2), #可调整图例的位置。##(1,1)第一个1是调整图例在图的内外(左右移动)，第二个1是在图中上下移动。
        #legend.position=c(0, 1.04), #legend.position=c(0,1)左上角，(1,1)是在右上角。
        panel.grid.major = element_blank(), #不显示网格线
        panel.grid.minor = element_blank()) #不显示网格线
boxplot_diff_TIICs
ggsave('02.DEcells.pdf',boxplot_diff_TIICs,w=10,h=6)
ggsave('02.DEcells.png',boxplot_diff_TIICs,w=10,h=6)
