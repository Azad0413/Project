rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/XA-0214-1/")
if (! dir.exists("./10_ssGSEA")){
  dir.create("./10_ssGSEA")
}
setwd("./10_ssGSEA")


dat <- read.delim2('../00_rawdata/dat(GSE97537).xls',row.names = 1)%>%lc.tableToNum()
group <- read.delim2("../00_rawdata/group(GSE97537).xls")
table(group$group)

## 03-1 ssGSEA--------
library(GSVA)
gene_set <- read.table("/data/nas1/luchunlin/pipeline/ssGSEA/mmc3.txt",
                       header = T,
                       sep ="\t")

library(homologene)
homologene::taxData
homologene <- homologene(gene_set$Metagene,inTax = 9606,outTax = 10116)
dim(homologene)
homologene <- homologene[,c(1,2)]
colnames(gene_set)
colnames(homologene) <- c('Metagene','Gene')
gene_set <- merge(gene_set,homologene,by='Metagene')
gene_set$Metagene <- gene_set$Gene
gene_set <- gene_set[,-4]

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
## 富集分数画热图
group<-group[order(group$group),]
ssgsea_score<-ssgsea_score[,group$sample]
annotation_col<-as.data.frame(group$group)
colnames(annotation_col)='Group'
rownames(annotation_col)=colnames(ssgsea_score)

color.key<-c("#3300CC", "#3399FF", "white", "#FF3333", "#CC0000")
ann_colors<-list(
  Group = c('control'="#00FFFF",'CIRI'="#FFAEB9"))
pdf(file = '01.heatmap.pdf',w=7.5,h=5)
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

png(file = '01.heatmap.png',w=580,h=400)
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

## 差异-------
dat.ssgsea <- ssgsea_score %>% 
  t %>% as.data.frame() %>% 
  tibble::rownames_to_column(var = "sample")
dat.ssgsea <- merge(group, dat.ssgsea, by = "sample")
dat.ssgsea2 <- tidyr::gather(dat.ssgsea, ImmuneCell, Score, -c("sample", "group"))
library(rstatix)
library(ggplot2)
library(ggpubr)
colnames(dat.ssgsea2)
stat_ssgsea <- dat.ssgsea2 %>% 
  group_by(ImmuneCell) %>%
  wilcox_test(Score ~ group) %>% 
  adjust_pvalue(method = "BH") %>%  # method BH == fdr
  add_significance("p")
write.table(stat_ssgsea,file = 'stat.ssgsea.xls',sep = '\t',row.names = F,quote = F)
DE.ssgsea<-stat_ssgsea[which(stat_ssgsea$p<0.05),]
write.table(DE.ssgsea,file = 'DE.ssgsea.xls',sep = '\t',row.names = F,quote = F)
colnames(dat.ssgsea2)
# violin.ssgsea<-dat.ssgsea2[dat.ssgsea2$ImmuneCell%in%stat_ssgsea$ImmuneCell[which(stat_ssgsea$p<0.05)],]
violin.ssgsea<-dat.ssgsea2
ssgsea_plot <- ggplot(violin.ssgsea, aes(x=ImmuneCell,
                                         y=Score,
                                         fill=group)) +
  # geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar",
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#FF6A6A", "#20B2AA"), name = "Group")+
  labs(title="Immune Cell", x="", y = "Score",size=20) +
  stat_compare_means(data = violin.ssgsea,
                     mapping = aes(group = group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F,label.x = 1.45) +
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
        panel.grid.minor = element_blank())#+facet_wrap(~ImmuneCell,scales = "free",nrow = 1) 
ssgsea_plot
ggsave(filename = '02.ssgsea.plot.pdf',ssgsea_plot,w=10,h=6)
ggsave(filename = '02.ssgsea.plot.png',ssgsea_plot,w=10,h=6)

### 相关性-------
library(Hmisc)
# diff_score <- ssgsea_score
diff_score<-ssgsea_score[DE.ssgsea$ImmuneCell,]
hubgene <- read.delim2('../04_model/hubgene.xls')
hub_exp <- dat[hubgene$symbol,]
nc<-t(rbind(diff_score,hub_exp))
m=rcorr(nc,type = 'spearman')$r[1:nrow(diff_score),(ncol(nc)-length(hubgene$symbol)+1):ncol(nc)]
m<-t(m)
p=rcorr(nc,type = 'spearman')$P[1:nrow(diff_score),(ncol(nc)-length(hubgene$symbol)+1):ncol(nc)]
p<-t(p)
library(dplyr)
library(dplyr)

tmp = matrix(case_when(p<0.0001~"****",
                       p<0.0001~"***",
                       p<0.01~"**",
                       p<0.05~"*",
                       T~""),nrow = nrow(p))
cor <- m
cor <- round(cor,digits = 3)
cor[abs(cor)<0.56] <- ''

textMatrix = paste(cor,"\n",
                   tmp, sep = "")
dim(textMatrix) = dim(m)
colnames(textMatrix)<-colnames(m)
rownames(textMatrix)<-rownames(m)
write.table(textMatrix,file = 'correlation.xls',sep = '\t',row.names = T,quote = F)
library(WGCNA)
pdf(file = '03.correlation.pdf',w=13,h=6)
par(mar = c(9, 8, 3, 1));
labeledHeatmap(Matrix = m, 
               xLabels = colnames(m), 
               yLabels = rownames(m), 
               cex.lab = 1, 
               # ySymbols = colnames(MEs_col)[-20], 
               colorLabels = FALSE, 
               colors = blueWhiteRed(50, 0.8),
               textMatrix = textMatrix, 
               setStdMargins = FALSE, 
               cex.text = 0.8, 
               zlim = c(-1,1),
               main = paste("Biomarker-TIICs correlation"))
dev.off()
png(file = '03.correlation.png',w=900,h=400)
par(mar = c(9, 8, 3, 1));
labeledHeatmap(Matrix = m, 
               xLabels = colnames(m), 
               yLabels = rownames(m), 
               cex.lab = 1, 
               # ySymbols = colnames(MEs_col)[-20], 
               colorLabels = FALSE, 
               colors = blueWhiteRed(50, 0.8),
               textMatrix = textMatrix, 
               setStdMargins = FALSE, 
               cex.text = 0.8, 
               zlim = c(-1,1),
               main = paste("Biomarker-TIICs correlation"))
dev.off()



