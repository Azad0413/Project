rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-370-8/")
if (! dir.exists("./10_ssGSEA")){
  dir.create("./10_ssGSEA")
}
setwd("./10_ssGSEA")

dat<-read.delim2("../03_DEGs/dat_final.xls", row.names = 1)  %>% lc.tableToNum
group <- read.delim2("../03_DEGs/group.xls")
RB.sample<-group$sample[which(group$group=='RB')]
control.sample<-group$sample[which(group$group=='control')]
## 03-1 ssGSEA--------
library(GSVA)
gene_set <- read.table("/data/nas1/luchunlin/pipeline/ssGSEA/mmc3.txt",
                       header = T,
                       sep ="\t")
dat.final2 <- as.matrix(dat)
gene_list <- split(as.matrix(gene_set)[,1],
                   gene_set[,2])

ssgsea_score = gsva(dat.final2, gene_list, 
                    method = "ssgsea", 
                    ssgsea.norm = TRUE, 
                    verbose = TRUE)
write.table(ssgsea_score,
            file = "ssgsea_result.xls",
            sep = "\t",
            quote = F)

## 03-2 差异免疫细胞鉴定-------
tiics_result <- ssgsea_score
pvalue = padj = log2FoldChange <- matrix(0, nrow(tiics_result), 1)
group_RB<-group[group$group=='RB',]
group_RB<-as.character(group_RB$sample)
group_control<-group[group$group=='control',]
group_control<-as.character(group_control$sample)

for (i in 1:nrow(tiics_result)){
  pvalue[i, 1] = p.value = wilcox.test(tiics_result[i, group_RB],
                                       tiics_result[i, group_control])$p.value
  log2FoldChange[i, 1] = mean(tiics_result[i, group_RB]) - 
    mean(tiics_result[i, group_control])
}
padj <- p.adjust(as.vector(pvalue), "fdr", n = length(pvalue))
rTable <- data.frame(log2FoldChange, 
                     pvalue, 
                     padj,
                     row.names = rownames(tiics_result))
control <- signif(apply(tiics_result[rownames(rTable), group_control], 
                        1,
                        mean), 4)
RB <- signif(apply(tiics_result[rownames(rTable), group_RB], 
                   1, 
                   mean), 4)
rTable <- data.frame(control, 
                     RB,
                     rTable[, c("padj", "pvalue", "log2FoldChange")])
rTable$immune_cell <- rownames(rTable)
rTable$sig <- ifelse(rTable$pvalue < 0.05,
                     ifelse(rTable$pvalue < 0.01, 
                            ifelse(rTable$pvalue < 0.001,
                                   ifelse(rTable$pvalue < 0.0001,
                                          paste(rTable$immune_cell, "****",  sep = ""),
                                          paste(rTable$immune_cell, "***", sep = "")),
                                   paste(rTable$immune_cell, "**", sep = "")),
                            paste(rTable$immune_cell, "*",  sep = "")), 
                     rTable$immune_cell)

diff_Table<-rTable[which(rTable$pvalue<0.05),]
## 11
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
#xCell3 <- xCell2[which(xCell2$pvalue<0.05),]
xCell3<-xCell2
diff_tiics <- rownames(xCell3)
violin_dat <- gather(xCell3, key=Group, value=score, -c("Immune_Cell","pvalue"))

violin_dat$Group <- ifelse(gsub("\\.","-",violin_dat$Group) %in% group_control,
                           "control", "RB") 
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
ggsave('01.TIICs.pdf',boxplot_diff_TIICs,w=10,h=6)
ggsave('01.TIICs.png',boxplot_diff_TIICs,w=10,h=6)

library(Hmisc)
rownames(ssgsea_score)<-rTable$immune_cell
diff_score<-ssgsea_score[diff_Table$immune_cell,]
hubgene <- read.delim2('../06_XGBoost/hubgene.xls')
hub_exp <- dat[hubgene$symbol,]
nc<-t(rbind(diff_score,hub_exp))
m=rcorr(nc)$r[1:nrow(diff_score),(ncol(nc)-length(hubgene$symbol)+1):ncol(nc)]
m<-t(m)
p=rcorr(nc)$P[1:nrow(diff_score),(ncol(nc)-length(hubgene$symbol)+1):ncol(nc)]
p<-t(p)
library(dplyr)
library(dplyr)
tmp = matrix(case_when(p<0.01~"**",
                       p<0.05~"*",
                       T~""),nrow = nrow(p))
source("modified_pheatmap.R")
pdf('02.heatmap.pdf',w=15,h=6)
pheatmap(m,
         display_numbers =tmp,
         angle_col =45,
         color = colorRampPalette(c("#92b7d1", "white", "#d71e22"))(100),
         border_color = "white",
         treeheight_col = 0,
         treeheight_row = 0,
         fontsize_row = 12,
         fontsize_col = 12,
         fontsize = 10)

dev.off()

png('02.heatmap.png',w=1400,h=500)
pheatmap(m,
         display_numbers =tmp,
         angle_col =45,
         color = colorRampPalette(c("#92b7d1", "white", "#d71e22"))(100),
         border_color = "white",
         treeheight_col = 0,
         treeheight_row = 0,
         fontsize_row = 12,
         fontsize_col = 12,
         fontsize = 10)

dev.off()
write.table(m,'corr.xls',
            sep = '\t',
            row.names = T)
write.table(p,'pvalue.xls',
            sep = '\t',
            row.names = T)


## 富集分数画热图
# group<-group[order(group$group),]
# ssgsea_score<-ssgsea_score[,group$sample]
# heat <- ssgsea_score[diff_Table$immune_cell,]
# annotation_col<-as.data.frame(group$group)
# colnames(annotation_col)='Group'
# rownames(annotation_col)=colnames(heat)
# library(pheatmap)
# color.key<-c("#3300CC", "#3399FF", "white", "#FF3333", "#CC0000")
# ann_colors<-list(
#   Group = c('control'="#00FFFF",'RB'="#FFAEB9"))
# pdf('02.heatmap.pdf',w=8,h=6)
# pheatmap(
#   heat,
#   color = colorRampPalette(color.key)(50),
#   border_color = 'darkgrey',
#   annotation_col = annotation_col,
#   annotation_colors = ann_colors,
#   labels_row = NULL,
#   clustering_method = 'ward.D2',
#   show_rownames = T,
#   show_colnames = F,
#   fontsize_col = 5,
#   cluster_cols = F,
#   cluster_rows = T)
# dev.off()
# png('02.heatmap.png',w=650,h=500)
# pheatmap(
#   ssgsea_score,
#   color = colorRampPalette(color.key)(50),
#   border_color = 'darkgrey',
#   annotation_col = annotation_col,
#   annotation_colors = ann_colors,
#   labels_row = NULL,
#   clustering_method = 'ward.D2',
#   show_rownames = T,
#   show_colnames = F,
#   fontsize_col = 5,
#   cluster_cols = F,
#   cluster_rows = T)
# dev.off()
# 
