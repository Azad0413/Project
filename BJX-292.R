#rm(list = ls())
# 01 差异分析---------
setwd("/data/nas1/luchunlin/project/BJX-292")
### C VS F 差异基因（6个样品）小鼠 
dat<-read_xlsx('/data/nas1/luchunlin/project/BJX-292/dat.xlsx')
dat<-as.data.frame(dat)
dat<-dat[!duplicated(dat$GeneName),]
rownames(dat)<-dat$GeneName
dat<-dat[,-1]
dat<-log2(dat+1)
group<-data.frame(sample = colnames(dat),
                  group=c(rep('C',3),rep('F',3)))
write.table(group,file = "group.xls",
            quote = F,
            sep = "\t",
            row.names = T)
## limma包 患者vs正常（p value<0.05，|logFC|>1）
## 分组矩阵
type<-group[,2]
design <- model.matrix(~ -1+factor(type,levels=c('C','F'))) 
colnames(design)<-c('C','F')
rownames(design)<-group$sample
library(limma)
# 对每一个基因进行线性模型构建
fit=lmFit(dat,design)
# 构建比较矩阵
contrast.matrix=makeContrasts(ControlVSMG=C-F,levels = design)
# 构建芯片数据线性模型，计算估计的相关系数和标准差
fit2=contrasts.fit(fit,contrast.matrix)
# 基于贝叶斯计算T值，F值和log-odds
fit2<-eBayes(fit2,0.01)
tempOutput = topTable(fit2, coef=1, n=Inf)
DEG= na.omit(tempOutput)
write.table(DEG,file = "DEG.xls",
            quote = F,
            sep = "\t",
            row.names = T)
# 筛选差异基因
logFC_cutoff <- 1
DEG$change = as.factor(
  ifelse(DEG$adj.P.Val < 0.05 & abs(DEG$logFC) > logFC_cutoff,
         ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
sig_diff <- subset(DEG,
                   DEG$adj.P.Val < 0.05 & abs(DEG$logFC) > logFC_cutoff)

dim(DEG)
dim(sig_diff)
# 1226    7
summary(sig_diff$change)
# DOWN  NOT   UP 
# 669    0  557 
write.table(DEG,file = "DEG_all.xls",
            quote = F,
            sep = "\t",
            row.names = T)
write.table(sig_diff,file = "DEG_sig.xls",
            quote = F,
            sep = "\t",
            row.names = T)


##  绘制火山图-----
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(Ipaper)
library(scales)
library(ggrepel)
library(EnhancedVolcano)

dat_rep<-DEG[rownames(DEG)%in%
               c(head(rownames(subset(sig_diff,sig_diff$logFC>6)),0),
                 head(rownames(subset(sig_diff,sig_diff$logFC< -2)),0)),]

volcano_plot<-EnhancedVolcano(DEG,lab = rownames(DEG),
                              x='logFC',
                              y='adj.P.Val',
                              pCutoff = 0.05,
                              maxoverlapsConnectors = Inf,
                              legendLabels = c('NS','Log2FC','adj.P.Value',
                                               'adj.P.Value & Log2FC'),
                              selectLab = c(rownames(dat_rep)),
                              drawConnectors = T,
                              ylim = c(0,7),
                              xlim = c(-6,10)
)+
  labs(y = "-log10 (adj.P.Value)")
volcano_plot
## 绘制热图-----
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(gplots)
group_rt<-group$group
group_rt<-as.data.frame(group_rt)
rt<-dat
colnames(group_rt)<-'group'
rownames(group_rt)<-group$sample
heat<-rt[rownames(rt)%in%
           c(head(rownames(subset(sig_diff,sig_diff$logFC>6)),10),head(rownames(subset(sig_diff,sig_diff$logFC< -2)),10)),]
#x<-log2(heat+1)
x<-t(scale(t(heat)))
ann_colors<-list(
  Group = c(Normal="lightblue",Tumor="darkorange"))
pheatmap(mat=x,
         annotation_col = group_rt,
         color=bluered(50),
         scale = "row",
         annotation_colors = ann_colors,
         fontsize = 10,
         show_colnames = F,
         cluster_cols = T,
         cluster_rows = T)


# 02 糖酵解相关基因---------
glycolysis<-read_xlsx('/data/nas1/luchunlin/project/BJX-292/Glycolysis_MSigDB.xlsx')
glycolysis<-as.data.frame(glycolysis)

## 将基因转换为小鼠的同源基因
# install.packages('homologene')
library(homologene)
homologene <- human2mouse(glycolysis$`Glycolysis-related_gene_set`)
dim(homologene)
# 288   4 
mouse_glycolysis_gene <- data.frame(homologene$mouseGene)
write.table(mouse_glycolysis_gene,file = "mouse_glycolysis_gene.xls",
            quote = F,
            sep = "\t",
            row.names = T)
## 取交集
gly_diff<-sig_diff[rownames(sig_diff)%in%mouse_glycolysis_gene$homologene.mouseGene,]
## 36个基因
write.table(gly_diff,file = "gly_diff.xls",
            quote = F,
            sep = "\t",
            row.names = T)
## 绘制热图-----
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(gplots)
group_rt<-group$group
group_rt<-as.data.frame(group_rt)
rt<-dat
colnames(group_rt)<-'group'
rownames(group_rt)<-group$sample
heat<-rt[rownames(rt)%in%rownames(gly_diff),]
#x<-log2(heat+1)
x<-t(scale(t(heat)))
ann_colors<-list(
  Group = c(Normal="lightblue",Tumor="darkorange"))
pheatmap(mat=x,
         annotation_col = group_rt,
         color=bluered(50),
         scale = "row",
         annotation_colors = ann_colors,
         fontsize = 10,
         show_colnames = F,
         cluster_cols = T,
         cluster_rows = T)





## GO/KEGG富集-----

library(clusterProfiler)
library(org.Mm.eg.db)
library(ggplot2)
library(enrichplot)
library(GOplot)
diff_gene_names <- rownames(gly_diff)
gene_transform <- bitr(diff_gene_names,
                       fromType = "SYMBOL",
                       toType = c("ENTREZID", "ENSEMBL", "REFSEQ"),
                       OrgDb = "org.Mm.eg.db")
ego <- enrichGO(gene = gene_transform$ENTREZID,
                OrgDb = org.Mm.eg.db,
                keyType = "ENTREZID",
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
write.table(ego,file = "GO.xls",sep = "\t",quote = F,row.names = F)
# 展示富集最显著的 GO term
go_bar <- barplot(ego, showCategory=7, split="ONTOLOGY") +
  facet_grid(ONTOLOGY ~ ., scales = "free",space = 'free')
go_bar
## KEGG富集分析（气泡图）
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "mmu",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)
kk <- setReadable(kk, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
write.table(kk@result,file = "KEGG.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk,showCategory=25)
kk_dot
##  hub gene  Pfkl Pygl Eno1 Tpi1 Gpi1
network<-read_xlsx('/data/nas1/luchunlin/project/BJX-292/network.xlsx')
homologene2 <- human2mouse(network$Target)
dim(homologene)
# 288   4 
mouse_glycolysis_gene <- data.frame(homologene$mouseGene)
## 03 GSVA------------
## 36个差异基因GSVA
library(GSVA)
library(GSEABase)
library(limma)
##
library(KEGGREST, quietly = TRUE)
library(tidyverse, quietly = TRUE)

# 返回信息很长，只取基因symbol.根据自己需要调整
symbolOnly <- function(x){
  items <- strsplit(x, ";", fixed = TRUE) %>% unlist()
  return(items[1])
}

# keggGet(x)[[1]]$GENE 数据基因名是个向量，其中奇数位置是 entrezgene_id 偶数位置是 symbol 
geneEntrez <- function(x){
  geneList <- keggGet(x)[[1]]$GENE
  if(!is.null(geneList)){
    listLength <- length(geneList)
    entrezList <- geneList[seq.int(from = 1, by = 2, length.out = listLength/2)]
    entrez <- stringr::str_c(entrezList, collapse = ",")
    return(entrez)
  }else{
    return(NA)
  }
  
}

# keggGet(x)[[1]]$GENE 数据基因名是个向量，其中奇数位置是 entrezgene_id 偶数位置是 symbol 
geneSymbol <- function(x){
  geneList <- keggGet(x)[[1]]$GENE
  if(!is.null(geneList)){
    listLength <- length(geneList)
    symbolList <- geneList[seq.int(from = 2, by = 2, length.out = listLength/2)] %>% map_chr(symbolOnly)
    symbol <- stringr::str_c(symbolList, collapse = ",")
    return(symbol)
  }else{
    return("")
  }
  
}

# 取得 hsaxxxxx 这种通路ID
pathwayID <- function(x){
  items <- strsplit(x, ":", fixed = TRUE) %>% unlist()
  return(items[2])
}

# 建议从这里开始读脚本。建议自己在交互模式下试一下用到的KEGGREST函数，看看返回数据的结构。
# 这是第一步，取得所有的KEGG通路列表
mouseList <- keggList("pathway", "mmu")
IDList <- names(keggList) %>% map_chr(pathwayID)

# 将通路ID和通路名放在一个表格(tibble)里
mousePathway <- tibble::tibble(pathway_id=IDList, pathway_name=mouseList)
head(mousePathway, n=3) %>% print()

# 用前面定义函数，获得每个通路的基因，然后也放在表格里
pathwayFull <- mousePathway %>% dplyr::mutate(entrezgene_id=map_chr(pathway_id, geneEntrez), hgnc_symbol=map_chr(pathway_id, geneSymbol))

# 保存数据
# write_tsv(pathwayFull, path="KEGGREST.tsv")
dim(pathwayFull) %>% print()

# 会有通路没有基因，我的话只需要有基因的，所以把无基因的移除
pathwayWithGene <- dplyr::filter(pathwayFull, !is.na(entrezgene_id) & hgnc_symbol != "")
write_tsv(pathwayWithGene, path="KEGGREST_WithGene.tsv")
dim(pathwayWithGene) %>% print()
pathwayWithGene<-pathwayWithGene[,-3]
class(pathwayWithGene)
pathwayWithGene<-as.data.frame(pathwayWithGene)
pathwayWithGene$pathway_name<-gsub(' - Mus musculus (mouse)','',pathwayWithGene$pathway_name,fixed = T)

result <- data.frame(col1="ID", col2="Pathway", col3 ="Gene")
n <- 1
for(i in 1:nrow(pathwayWithGene)){
  col1 <- pathwayWithGene[i,1]
  col2 <- pathwayWithGene[i,2]
  col3 <- pathwayWithGene[i,3]
  genelist <- unlist(strsplit(col3, ","))
  for(gene in genelist){
    result[n, 1] <- col1
    result[n, 2] <- col2
    result[n, 3] <- gene
    n <- n + 1
  }
}
mouse_KEGG<-result[,c(2:3)]
colnames(mouse_KEGG)<-c('','gene_symbol')


mouse_KEGG$gs_name<-as.character(mouse_KEGG$gs_name)
mouse_KEGG$gene_symbol<-as.character(mouse_KEGG$gene_symbol)

mouse_KEGGSet<-mouse_KEGG%>%split(x = .$gene_symbol, f = .$gs_name)

#library(org.Mm.eg.db)
#library(tidyverse)
#install.packages('msigdbr')
#library(msigdbr)
#msigdbr_species()
#mouse<-msigdbr(species = "Mus musculus")
#mouse[1:5,1:5]
#table(mouse$gs_cat) ## 查看目录，与MSigDB一样包含9个数据集
## 查看需要的数据集在哪里
#table(mouse$gs_subcat)
#mouse_KEGG<-msigdbr(species = "Mus musculus",
#                    category = 'C2',
#                    subcategory = "CP:KEGG")%>%
#  dplyr::select(gs_name,gene_symbol)
#head(mouse_KEGG)
#class(mouse_KEGG)
#mouse_KEGGSet<-mouse_KEGG%>%split(x = .$gene_symbol, f = .$gs_name)

## 将表达矩阵提取出来
gsva_exp<-dat[rownames(gly_diff),]

write.table(gsva_exp,
            file = "gsva_exp.xls",
            quote = F,
            sep = "\t",
            row.names = T)

group2<-group
group2<-group2$group%>%as.factor()
design2<-model.matrix(~0+group2)
rownames(design2)<-colnames(gsva_exp)
colnames(design2)<-levels(group2)
compare<-makeContrasts('C-F',levels = design2)
# KEGG

#KEGG_ref <- getGmt("/data/nas1/luchunlin/pipeline/GSVA/c2.cp.kegg.v7.4.symbols.gmt")
es_KEGG <- gsva(as.matrix(gsva_exp), mouse_KEGGSet,
                min.sz=2, max.sz=500)
es_KEGG <- as.data.frame(es_KEGG)
fit <- lmFit(es_KEGG, design2)
fit2 <- contrasts.fit(fit ,compare)
fit3 <- eBayes(fit2)
allGeneSets <- topTable(fit3, coef = 1, number = Inf)

logFCcutoff <- 0
allGeneSets$change = as.factor(
  ifelse(allGeneSets$P.Value < 0.05 & abs(allGeneSets$logFC) > logFCcutoff,
         ifelse(allGeneSets$logFC > logFCcutoff,'UP','DOWN'),'NOT')
)

DEGeneSets <- subset(allGeneSets,
                     allGeneSets$P.Value < 0.05 & abs(allGeneSets$logFC) > 0)


write.table(allGeneSets,
            file = "GSVA_KEGG(C_vs_F).xls",
            quote = F,
            sep = "\t",
            row.names = T)
write.table(DEGeneSets,
            file = 'diff_KEGG.xls',
            quote = F,
            sep = "\t",
            row.names = T)
### 发散条形图绘制

#install.packages('ggprism')
library(ggprism)
## barplot
dat_plot<-data.frame(id=rownames(allGeneSets),
                     t=allGeneSets$t)
dat_plot$threshold = factor(ifelse(dat_plot$t  >-1, ifelse(dat_plot$t >= 1 ,'Up','NoSignifi'),'Down'),levels=c('Up','Down','NoSignifi'))

dat_plot<-dat_plot%>%arrange(t)

dat_plot$id<-factor(dat_plot$id,levels=dat_plot$id)

p <- ggplot(data = dat_plot,aes(x = id,y = t,fill = threshold)) +
  geom_col()+
  coord_flip() +
  scale_fill_manual(values = c('Up'= '#36638a','NoSignifi'='#cccccc','Down'='#7bcd7b')) +
  geom_hline(yintercept = c(-1,1),color = 'white',size = 0.5,lty='dashed') +
  xlab('') + 
  ylab('t value of GSVA score, tumour versus non-malignant') + #注意坐标轴旋转了
  guides(fill=F)+ # 不显示图例
  theme_prism(border = T) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
p
# 小于-2的数量
low1 <- dat_plot %>% filter(t < -1) %>% nrow()
# 小于0总数量
low0 <- dat_plot %>% filter( t < 0) %>% nrow()
# 小于2总数量
high0 <- dat_plot %>% filter(t < 1) %>% nrow()
# 总的柱子数量
high1 <- nrow(dat_plot)

# 依次从下到上添加标签
p <- p + geom_text(data = dat_plot[1:low1,],aes(x = id,y = 0.1,label = id),
                   hjust = 0,color = 'black') + # 小于-1的为黑色标签
  geom_text(data = dat_plot[(low1 +1):low0,],aes(x = id,y = 0.1,label = id),
            hjust = 0,color = 'grey') + # 灰色标签
  geom_text(data = dat_plot[(low0 + 1):high0,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'grey') + # 灰色标签
  geom_text(data = dat_plot[(high0 +1):high1,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'black') # 大于1的为黑色标签
p
