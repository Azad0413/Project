# 一、训练集(pre) --------------------------------------------------------
### TCGA和GTEx数据集合并
# 01-TCGA数据的获取 ------------------------------------------------------------
library(rtracklayer)
library(TCGAbiolinks)
library(dplyr)
library(DT)
library(SummarizedExperiment)
library(readr)
library(tidyr)
# 01-1下载原始Counts数据 --------------------------------------------------------------
#TCGAbiolinks:::getProjectSummary("TARGET-OS")
#query <- GDCquery(project = "TARGET-OS",
#                  data.category = "Transcriptome Profiling",
#                  data.type = "Gene Expression Quantification",
#                  workflow.type = "STAR - Counts") ## 使用Count数据
#GDCdownload(query)
#TCGA_dat <- GDCprepare(query)
#TCGA_dat <- as.data.frame(assay(TCGA_dat)) ## TCGA:TARGET-OS表达矩阵
TCGA_dat<-read_tsv(file = 'TARGET-OS.htseq_counts.tsv.gz')
library(tidyverse)
colnames(TCGA_dat)
TCGA_dat<-column_to_rownames(TCGA_dat,var = 'Ensembl_ID')

## xena下载的数据经过了log2+1转化，需要将其还原
TCGA_dat<-2^TCGA_dat-1
# 01-2原始数据预处理 -------------------------------------------------------------
### 查看原始count数据是否存在基因表达量全部为0的行
View(TCGA_dat[which(rowSums(TCGA_dat)==0),]) 
### 表达量不为0的矩阵(去除基因表达量为0的行)
TCGA_dat <- TCGA_dat[which(rowSums(TCGA_dat)>0),] 
### 将行名变成第一列,命名为Ensembl_ID
TCGA_dat$Ensembl_ID <- rownames(TCGA_dat)
rownames(TCGA_dat) <- NULL
TCGA_dat <- TCGA_dat[,c(89,1:88)]
### 调整格式:去除TCGA_dat文件Ensembl_ID列Ensembl号后面的点号
TCGA_dat <- separate(TCGA_dat,
                     Ensembl_ID,
                     into=c('Ensembl_ID'),
                     sep='\\.')
### 将TCGA_dat按照Ensembl_ID列去重
aggr_TCGA_dat <- aggregate(TCGA_dat[,2:89],
                           by=list(Ensembl_ID=TCGA_dat$Ensembl_ID),
                           mean) 
### 注:aggr_TCGA_dat是Ensembl去0再去重之后的矩阵

# 02-GTEx数据的获取 -------------------------------------------------------------
# 02-1上传、读取原始Counts数据 -------------------------------------------------------------
## 1.本地下载Counts数据，上传至服务器
### 文件名为'GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz'
GTEx <- read.table('GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz',
                   header = T,
                   sep = '\t',
                   skip = 2)

# 02-2上传、读取注释文件 -----------------------------------------------------------
## 1.本地下载注释文件，上传至服务器
GTEx_anno <- read.table('GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt',
                        header = T,sep = '\t',quote = '')
table(GTEx_anno$SMTS) ## 查看各个组织部位有多少个样本，Muscle共1132个样本
## 2.从注释文件b中提取肌肉组织的样本ID，获得表达矩阵
### 从注释文件GTEx_anno中提取肌肉组织的样本ID
GTEx_muscle <- GTEx_anno[which(GTEx_anno$SMTS == "Muscle"),] 
### 将GTEx_muscle文件中SAMPID列的'-'替换成'.'
a <- gsub('-', '.', GTEx_muscle$SAMPID)
GTEx_muscle <- cbind(a, GTEx_muscle)
GTEx_muscle <- GTEx_muscle[,-2]
colnames(GTEx_muscle)[1] <- 'SAMPID'
### 取交集
cols <- intersect(GTEx_muscle$SAMPID, colnames(GTEx))
GTEx_dat <- subset(GTEx, select = cols)
Name <- GTEx$Name
GTEx_dat <- cbind(Name, GTEx_dat)
### 调整格式:去除GTEx_dat文件Name列Ensembl号后面的点号
GTEx_dat <- separate(GTEx_dat,
                     Name,
                     into=c('Name'),
                     sep='\\.')

### 查看原始count数据是否存在基因表达量全部为0的行
View(GTEx_dat[which(rowSums(GTEx_dat[,2:804]) == 0),]) 
### 共有2162个Ensembl在全部803个肌肉组织样本中的表达量都是0
### 原因是，从原始矩阵提取肌肉组织数据，一些基因在其它组织中表达，而在肌肉组织中不表达
### 所以要去除全部表达量数值为0的基因行
### 表达量不为0的矩阵(去除基因表达量为0的行),即:56200-2162=54038
GTEx_dat <- GTEx_dat[which(rowSums(GTEx_dat[,2:804]) > 0),] 

### 将GTEx_dat按照Name列去重,这个时候Name在第一列,不是行名
### 在行名的时候,有重复的话会自动+1,会出现问题
aggr_GTEx_dat <- aggregate(GTEx_dat[,2:804],
                           by=list(Ensembl_ID=GTEx_dat$Name),
                           mean) 
## 注:aggr_GTEx_dat是Ensembl去重且去0之后的矩阵

# 03-TCGA和GTEx数据集合并和ID转换 ----------------------------------------------------------

# 03-1数据集合并 ---------------------------------------------------------------
merge <- merge(aggr_GTEx_dat,aggr_TCGA_dat,by = 'Ensembl_ID')
### 将第一列变成行名
rownames(merge) <- merge[,1]
merge <- merge[,-1]
### 注:train_set_pre是TCGA和GTEx合并后的Ensembl_ID表达矩阵,即训练集,但是需要去除批次效应

# 03-2ID转换 --------------------------------------------------------------------
## 1.注释信息的下载和整理
library(stringr) 
library(org.Hs.eg.db)
### 提取出SYMBOL信息
g2s <- unique(toTable(org.Hs.egSYMBOL))
### 提取出ENSEMBL信息
g2e <- unique(toTable(org.Hs.egENSEMBL)) 
### 整合两个数据，得到ENSEMBL和SYMBOL的对应关系
s2e <- merge(g2e,g2s,by='gene_id')
s2e <- s2e[,-1] ## 去除gene_id列
### s2e是ENSEMBL和SYMBOL的对应关系

## 2.merge和s2e合并
### 将merge矩阵的行名转换成第一列,命名为ensembl_id,即和s2e文件保持一致
merge$ensembl_id <- rownames(merge)
merge <- merge[,c(892,1:891)]
rownames(merge) <- NULL
merge <- merge(s2e, merge, by = 'ensembl_id') ### 合并
merge <- merge[,-1] ### 删除'ensembl_id'列
### 根据symbol列去重,取均值
merge <- aggregate(merge[,2:892],
                   by=list(symbol = merge$symbol),
                   mean) 
rownames(merge) <- merge[,1] 
merge <- merge[,-1] ### (32576, 891)
### 注: merge是经过ID转换后的合并矩阵,并且经过去重,没有一行表达量全是0的基因

# 03-3去除批次效应 --------------------------------------------------------------
library(sva)
# BiocManager::install("bladderbatch")
library(bladderbatch)
## 1.样本分组
group_GTEx <- data.frame(sample = colnames(merge[1:803]),
                         group = 'normal')
group_TCGA <- data.frame(sample = colnames(merge[804:891]),
                         group = 'tumor')
group_list<- rbind(group_GTEx, group_TCGA)

## 2.设置批次
batch_GTEx <- data.frame(sample = colnames(aggr_GTEx_dat[,2:804]),
                         batch = 1)
batch_TCGA <- data.frame(sample = colnames(aggr_TCGA_dat[,2:89]),
                         batch = 2)
batch <- rbind(batch_GTEx, batch_TCGA)
batch <- merge(batch, group_list, by = 'sample')
### 注:batch含样本名、分组和批次,共3列信息

## 3.去除批次效应
### 本次分析过程中,去除批次效应不添加分组信息(协变量)
Combat_merge <- ComBat_seq(as.matrix(merge),
                           batch = batch$batch)
Combat_merge <- as.data.frame(Combat_merge)
train_set_pre <- Combat_merge
#### train_set_pre是去除批次效应之后的表达矩阵,表达量是原始count值,未经过log转换

# 二、验证集(pre) -------------------------------------------------------------------
## 1.获取GEO原始表达矩阵和临床信息
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 2)
library('GEOquery')
eSet <- getGEO('GSE21257', 
               destdir = '.', 
               getGPL = F)
exp <- exprs(eSet[[1]]) ## 获得表达矩阵exp，但行名为探针名
pheno <- pData(eSet[[1]]) ## 获得临床信息

## 2.ID转换
### 下载芯片注释信息
GPL10295 <- getGEO("GPL10295", destdir=".")
### 转换成数据框格式
GPL10295 <- Table(GPL10295) #转换成数据框形式
GPL10295 <- GPL10295[,c(1,6)] ## 提取'ID'和'Symbol'两列

## 3.表达矩阵和注释信息合并,获得validation_set_pre原始表达量矩阵
### 表达矩阵(exp)格式调整:将exp的行名变成第一列,命名为ID,和GPL10295保持一致
exp <- as.data.frame(exp)
exp$ID <- rownames(exp)
exp <- exp[,c(54,1:53)] ## 将ID列移到第一列
rownames(exp) <- NULL

### 将表达矩阵(exp)和注释信息(GPL10295)按照ID列合并
validation_set_pre <- merge(GPL10295, exp, by = 'ID')

### validation_set_pre格式调整
#### 去除ID列
validation_set_pre <- validation_set_pre[,-1] 
#### 删除Symbol为空的行
validation_set_pre <- validation_set_pre[validation_set_pre$Symbol!='', ] 
#### 查看是否存在表达量都是0的行
View(validation_set_pre[which(rowSums(validation_set_pre[,2:54])==0),]) #### 不存在
#### 查看Symbol列是否有重复
##### 存在较多重复,去重(30549 ~ 24998)
validation_set_pre <- aggregate(validation_set_pre[,2:54],
                                by=list(symbol=validation_set_pre$Symbol),
                                mean)
rownames(validation_set_pre) <- validation_set_pre[,1]
validation_set_pre <- validation_set_pre[,-1]
#### 注:validation_set_pre为原始表达量矩阵,没有经过log转换

# 三、训练集(pre)和验证集(pre)取交集 --------------------------------------------------
tv_inter <- intersect(rownames(train_set_pre), rownames(validation_set_pre)) 
## 交集基因有13696个

# 四、训练集 -------------------------------------------------------------------
## 从训练集(pre)中提取交集基因tv_inter,形成训练集矩阵
train_set <- train_set_pre[tv_inter,]
## 注:(13696, 891)

# 五、验证集 -------------------------------------------------------------------
## 从验证集(pre)中提取交集基因tv_inter,形成验证集矩阵
validation_set <- validation_set_pre[tv_inter,]
## 注:(13696, 53)

# 六、DEBRGs的鉴定 ----------------------------------------------------------------

# 01-BRGs表达矩阵的获取 ----------------------------------------------------------
## 1.本地上传、读取BMP基因集
BMP_set <- read.csv('BMP-Genecard.csv', header = T)
## 2.获取BRGs基因列表
BRGs <- BMP_set$Gene.Symbol 
### 注:N = 2770
## 3.形成BRGs表达矩阵
BRGs_dat <- train_set[BRGs,]
BRGs_dat <- na.omit(BRGs_dat)

# 02-差异分析 -----------------------------------------------------------------
library(DESeq2)
## 去除BRGs_dat表达矩阵中的小数,转换成整数
BRGs_dat1  <- as.data.frame(apply(BRGs_dat, 2, function(x) as.integer(x)))
rownames(BRGs_dat1) <- rownames(BRGs_dat)
condition <- factor(batch$group)
colData <- data.frame(row.names = colnames(BRGs_dat1), condition)

## 1.构建dds矩阵
dds <- DESeqDataSetFromMatrix(as.matrix(BRGs_dat1), 
                              colData = colData,
                              design = ~ condition) 
## 2.对dds矩阵进行标准化
dds <- DESeq(dds) 
nrDEG <- results(dds)
nrDEG <- as.data.frame(nrDEG)
write.table(nrDEG, 
            file = "nrDEG.xls",
            quote = F,
            sep = "\t",
            row.names = T)
## 3.提取差异分析结果
logFC_cutoff <- 1
nrDEG$change = as.factor(
  ifelse(nrDEG$padj < 0.05 & abs(nrDEG$log2FoldChange) > logFC_cutoff,
         ifelse(nrDEG$log2FoldChange > logFC_cutoff ,'UP','DOWN'),'NOT')
)
sig_diff <- subset(nrDEG,
                   nrDEG$padj < 0.05 & abs(nrDEG$log2FoldChange) >= 1)

write.table(sig_diff, 
            file = "sig_diff.xls",
            quote = F,
            sep = "\t",
            row.names = T)
### 注:共鉴定出46个DEBRGs

## 4.火山图
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(scales)
library(ggrepel) 
nrDEG_volcano <- nrDEG
nrDEG_volcano <- nrDEG_volcano[order(-log10(nrDEG_volcano$padj),decreasing = T) ,]
nrDEG_volcano$anno_name <- rownames(nrDEG_volcano)
nrDEG_volcano$anno_name[11:nrow(nrDEG_volcano)] <- NA
volcano_plot <- ggplot(data = nrDEG_volcano, 
                       aes(x = log2FoldChange,
                           y = -log10(padj),
                           colour = change))+
  geom_point(alpha = 0.7)+
  scale_color_manual(values=c("blue", "darkgray","red"))+
  geom_text_repel(aes(label=anno_name), show.legend = F, segment.colour = 'black', size = 3)+
  geom_vline(xintercept = c(-1,1),
             lty = 4,
             col = "darkgray",
             lwd = 0.6)+
  geom_hline(yintercept = -log10(0.05),
             lty = 4,
             col = "darkgray",
             lwd = 0.6)+
  theme_classic(base_line_size = 0.5)+
  theme(axis.title.x = element_text(size = 12, 
                                    color = "black",
                                    face = "bold"),
        axis.title.y = element_text(size = 12,
                                    color = "black",
                                    face = "bold", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        legend.title = element_blank(),
        legend.text = element_text(color="black", # 设置图例标签文字
                                   size = 8, 
                                   face = "bold"),
        axis.text.x = element_text(size = 12,  # 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0), #角度
        axis.text.y = element_text(size = 12,  
                                   color = "black",
                                   face = "bold", 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0))
volcano_plot

## 5.热图
## 5.1 数据准备
## 获取DEBRGs基因列表
DEBRGs <- rownames(sig_diff) ### 46个DEBRGs
## 从train_set中提取DEBRGs的表达量,形成DEBRGs矩阵以进行热图的绘制
DEBRGs_exp <- train_set[DEBRGs,]

## 5.2 绘制热图
library(pheatmap)
library(gplots)
library(RColorBrewer)
library(ggplot2)
# pheatmap(log2(DEBRGs_exp+1))
## 对原始基因表达量取log,因为个别基因的表达量为2万多,个别为0
DEBRGs_exp1 <- log2(DEBRGs_exp+1)
## 样本分组
## 分组1:对列(样本)进行分组
df_info <- group_list
rownames(df_info) <- df_info[,1]
df1 <- df_info[,-1,drop=FALSE]
## 分组2:对行(基因)进行分组
df2 <- sig_diff
df2$gene <-rownames(df2) 
df2 <- df2[,c(8,7)]
df2 <- df2[,-1,drop=FALSE]

bk = unique(c(seq(-2,2, length=100)))
pheatmap(DEBRGs_exp1,
         cellwidth = 0.5, cellheight = 9, ## 小方块的长和框
         breaks = bk,
         annotation_col = df1,
         annotation_row = df2,
         scale = 'row',
         fontsize_row = 8,
         cluster_row = T,
         cluster_col = F,
         show_rownames = T,
         show_colnames = F,
         border_color = NA,
         drop_levels = T,
         color = colorRampPalette(c("blue1", "white", "red"))(100),) 

# 03-富集分析 -----------------------------------------------------------------
library('BiocManager')
library('AnnotationHub')
library('readxl')
library('org.Hs.eg.db')
library('clusterProfiler')
library('GOplot')
library('ggplot2')  
library('stringr')  
## 获取DEBRGs的ENTREZID,即将SYMBOL转换成ENTREZID
DEBRGs_ENTREZID <- bitr(DEBRGs, fromType = 'SYMBOL',toType = 'ENTREZID','org.Hs.eg.db')

## 1.GO分析
# ALL <- enrichGO(gene = DEBRGs_ENTREZID$ENTREZID,
#                 keyType = 'ENTREZID',
#                 ont = "ALL",
#                 OrgDb = org.Hs.eg.db,
#                 pAdjustMethod = 'BH',
#                 pvalueCutoff = 0.05,
#                 qvalueCutoff = 0.2,
#                 readable = T)
# head(ALL)
# write.table(ALL,file = "GO.xls",sep = "\t",quote = F,row.names = T)
# 
# ### GO条形图
# barplot(ALL, split="ONTOLOGY") + 
#   facet_grid(ONTOLOGY~.,scale="free",space='free')

### 20211125客户反馈,补充弦图
# GO弦图
GO_all <- enrichGO(gene = DEBRGs_ENTREZID$ENTREZID,
                   keyType = 'ENTREZID',
                   ont = "ALL",
                   OrgDb = org.Hs.eg.db,
                   pAdjustMethod = 'BH',
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.2,
                   readable = T)
GO_res <- GO_all@result
write.table(GO_res,file = "DEBRGs_GO.xls",sep = "\t",quote = F,row.names = T)
GO <- GO_res[order(GO_res$p.adjust,decreasing = F),]
GO <- GO[1:10,c(1,2,3,9,7)]
GO$geneID <- str_replace_all(GO$geneID,"/",",") ## 替换富集条目中基因一栏的分隔符为逗号,这是GOplot默认的
names(GO) <- c("Category","ID","term","Genes","adj_pval")

## 添加差异倍数logFC值
logFC_dat <- sig_diff[DEBRGs,2,drop = F]
colnames(logFC_dat) <- "logFC"
logFC_dat$ID <- rownames(logFC_dat)
logFC_dat <- logFC_dat[DEBRGs_ENTREZID$SYMBOL,]
rownames(logFC_dat) <- NULL
logFC_dat <- logFC_dat[,2:1]
circ <- circle_dat(GO,logFC_dat)
chord <- chord_dat(data = circ, genes = logFC_dat, process = GO$term)
GOChord(chord, 
        space = 0.02,
        gene.order = 'logFC', 
        gene.space = 0.25, 
        gene.size = 5)
# 20:21

## 2.KEGG分析
KEGG<-enrichKEGG(
  gene = DEBRGs_ENTREZID$ENTREZID,
  organism = "hsa",
  keyType = "kegg",
  pvalueCutoff = 0.2,
  pAdjustMethod = "BH",
  qvalueCutoff = 0.2)
head(KEGG)
write.table(KEGG@result,file = "KEGG.csv",sep = "\t",quote = F,row.names = F)
## 气泡图
dotplot(KEGG,showCategory = 10,title = 'KEGG Pathway')

# 七、单因素Cox回归分析初步筛选预后相关BRGs ----------------------------------------------------------------
## 单因素Cox分析需要去除正常样本,使用FPKM文件
## 本次分析就是仅使用TCGA数据集,但需下载FPKM文件
library(SummarizedExperiment)
library(survival)
library(survminer)
library(TCGAbiolinks)

# 01-TCGA(TARGET-OS)数据集的FPKM文件的下载和预处理 ------------------------------------------
## 1.原始FPKM表达数据的下载
#geneexp <- assay(mydata,i = "fpkm_unstranded")
#query_FPKM <- GDCquery(project = "TARGET-OS",
#                       data.category = "Transcriptome Profiling",
#                       data.type = "Gene Expression Quantification",
#                       workflow.type = "STAR - FPKM") ## 使用FPKM数据
#GDCdownload(query_FPKM)
#TCGA_dat_FPKM <- GDCprepare(query_FPKM)
#TCGA_dat_FPKM <- as.data.frame(assay(TCGA_dat_FPKM)) 
TCGA_dat_FPKM<-read_tsv(file = 'TARGET-OS.htseq_fpkm.tsv.gz')
colnames(TCGA_dat_FPKM)
TCGA_dat_FPKM<-column_to_rownames(TCGA_dat_FPKM,var = 'Ensembl_ID')
### (56602,88) 56602个ENSGxxx,88个样本

## 2.数据预处理
## 查看原始FPKM数据是否存在基因表达量全部为0的行
View(TCGA_dat_FPKM[which(rowSums(TCGA_dat_FPKM)==0),]) ## N=1654
## 表达量不为0的矩阵(去除基因表达量为0的行)
TCGA_dat_FPKM <- TCGA_dat_FPKM[which(rowSums(TCGA_dat_FPKM)>0),] ## 56602-1654=54948
TCGA_dat_FPKM<-2^TCGA_dat_FPKM-1
## 将行名变成第一列,命名为Ensembl_ID
TCGA_dat_FPKM$Ensembl_ID <- rownames(TCGA_dat_FPKM)
rownames(TCGA_dat_FPKM) <- NULL
TCGA_dat_FPKM <- TCGA_dat_FPKM[,c(89,1:88)]
TCGA_dat_FPKM <- separate(TCGA_dat_FPKM,
                     Ensembl_ID,
                     into=c('Ensembl_ID'),
                     sep='\\.')
## ID转换:将TCGA_dat_FPKM按照Ensembl_ID列去重
aggr_TCGA_dat_FPKM <- aggregate(TCGA_dat_FPKM[,2:89],
                                by=list(Ensembl_ID=TCGA_dat_FPKM$Ensembl_ID),
                                mean) 
## 注:aggr_TCGA_dat_FPKM是Ensembl去0再去重之后的矩阵

## 将aggr_TCGA_dat_FPKM的Ensembl_ID转换成symbol
colnames(aggr_TCGA_dat_FPKM)[1] <- 'ensembl_id'
TCGA_dat_FPKM1 <- merge(s2e,aggr_TCGA_dat_FPKM,by='ensembl_id')
TCGA_dat_FPKM1 <- TCGA_dat_FPKM1[,-1] ## 去除ensembl_id列 (33878)
## symbol列去重
TCGA_dat_FPKM1 <- aggregate(TCGA_dat_FPKM1[,2:89],
                            by=list(symbol=TCGA_dat_FPKM1$symbol),
                            mean) 
rownames(TCGA_dat_FPKM1) <- TCGA_dat_FPKM1[,1]
TCGA_dat_FPKM1 <- TCGA_dat_FPKM1[,-1]
## TCGA_dat_FPKM1是去0去重之后的symbol名FPKM表达矩阵

# 02-DEBRGs_FPKM表达矩阵 ------------------------------------------------------
## 1.从TCGA_dat_FPKM1文件中提取DEBRGs的表达量,形成小矩阵,以进行后续的单因素Cox分析
sig_diff<-read_xlsx('sig_diff.xlsx')
DEBRGs<-sig_diff$Symbol
DEBRGs_FPKM_exp <- TCGA_dat_FPKM1[DEBRGs,] ## 46个基因,46行
## 2.将FPKM表达量进行log转换
DEBRGs_FPKM_exp <- log2(DEBRGs_FPKM_exp+1)
## 注:DEBRGs_FPKM_exp是显著差异表达BRGs的FPKM表达矩阵,FPKM值经过log转换

# 03-生存信息的下载和整理 -----------------------------------------------------------
## 1.从UCSC Xena下载OS患者生存信息:TARGET-OS.survival.tsv,本地上传、读取
TARGET_OS_survival <- read.table('TARGET-OS.survival.tsv',
                                 header = T,
                                 sep = '\t',
                                 quote = '')
# 04-将DEBRGs_FPKM表达矩阵和生存信息进行合并 ---------------------------------------------------
## 1.DEBRGs_FPKM表达矩阵数据框格式调整
DEBRGs_FPKM_exp <- t(DEBRGs_FPKM_exp) ## 数据框转置
DEBRGs_FPKM_exp <- as.data.frame(DEBRGs_FPKM_exp)
## 将DEBRGs_FPKM_exp行名变成第一列,命名为sample,和生存信息文件(TARGET_OS_survival)列名保持一致
DEBRGs_FPKM_exp$sample <- rownames(DEBRGs_FPKM_exp)
DEBRGs_FPKM_exp <- DEBRGs_FPKM_exp[,c(47,1:46)]
rownames(DEBRGs_FPKM_exp) <- NULL
## 将DEBRGs_FPKM_exp文件中sample列中的-01R替换为空
DEBRGs_FPKM_exp$sample <- gsub(pattern = '-01R',
                               replacement = '',
                               x = DEBRGs_FPKM_exp$sample)

## 如果基因名存在"-",会导致错误,将"-"更改为"_"
colnames_sum <- colnames(DEBRGs_FPKM_exp)
colnames_sum <- gsub("-","_",colnames_sum)
colnames_sum <- gsub(" ","",colnames_sum)
colnames(DEBRGs_FPKM_exp) <- colnames_sum

## 2.生存信息文件格式调整
TARGET_OS_survival <- TARGET_OS_survival[,-3] ## 将生存数据TARGET_OS_survival中X_PATIENT列去除

## 3.合并表达矩阵和生存信息
DEBRGs_FPKM_exp_survival <- merge(TARGET_OS_survival,
                                  DEBRGs_FPKM_exp,
                                  by='sample')
## 只有85个样本有生存信息
rownames(DEBRGs_FPKM_exp_survival) <- DEBRGs_FPKM_exp_survival[,1]
DEBRGs_FPKM_exp_survival <- DEBRGs_FPKM_exp_survival[,-1]

## DEBRGs_FPKM_exp_survival文件是合并了DEBRGs FPKM表达量(log转换后)和生存信息的矩阵

# 05-单因素Cox ---------------------------------------------------------------
library("survival")
library("survminer")
## 对46个变量(covariates)进行单因素Cox分析
colnames_sum <- colnames(DEBRGs_FPKM_exp_survival)
covariates <- colnames_sum[-which(colnames_sum %in% c("OS", "OS.time"))]
## 对每一个变量构建生存分析的公式
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(OS.time, OS)~', x)))
## 循环对每一个变量做Cox回归分析
univ_models <- lapply(univ_formulas, 
                      function(x){coxph(x, data = DEBRGs_FPKM_exp_survival)})
univ_models$DLX2
## #提取HR，95%置信区间和p值
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
## 转换成数据框，并转置
res <- t(as.data.frame(univ_results, check.names = FALSE))
res <- as.data.frame(res)
## 将结果文件导出
write.table(file="univariate_cox_result.xls",
            res,
            quote=F,
            sep="\t",
            row.names = T)
## 筛选p.value < 0.05的DEBRGs,即这些基因与生存显著相关(初步筛选),N=6
res_results_0.05 <- res[which(as.numeric(res$p.value) < 0.05),]
res_results_0.05 <- na.omit(res_results_0.05)

### 20211125客户反馈,补充森林图
# 单因素cox森林图
# 因为要绘制和后面保持一致的森林图,所以先重做一下单因素cox
uniTab <- data.frame()
for(i in colnames(DEBRGs_FPKM_exp_survival[,3:ncol(DEBRGs_FPKM_exp_survival)])){
  cox <- coxph(Surv(OS.time, OS) ~ DEBRGs_FPKM_exp_survival[,i], data = DEBRGs_FPKM_exp_survival)
  coxSummary = summary(cox)
  uniTab = rbind(uniTab,
                 cbind(id = i,
                       HR = coxSummary$conf.int[,'exp(coef)'],
                       HR.95L = coxSummary$conf.int[,'lower .95'],
                       HR.95H = coxSummary$conf.int[,'upper .95'],
                       pvalue = coxSummary$coefficients[,'Pr(>|z|)']))
}
rownames(uniTab) <- uniTab[,1]
uniTab$HR <- as.numeric(uniTab$HR)
uniTab$HR.95L <- as.numeric(uniTab$HR.95L)
uniTab$HR.95H <- as.numeric(uniTab$HR.95H)
uniTab$pvalue <- as.numeric(uniTab$pvalue)
write.table(uniTab, file = 'uniCox.txt', sep = '\t', row.names = F, quote = F)
# 绘制森林图
bioForest = function(coxFile = null, forestFile = null, forestCol = null){
  rt <- read.table(coxFile, header = T, sep = '\t', check.names = F, row.names = 1)
  Variables <-  rownames(rt)
  hr <- sprintf('%.3f', rt$'HR')
  hrLow <- sprintf('%.3f', rt$'HR.95L')
  hrHigh <- sprintf('%.3f', rt$'HR.95H')
  Hazard_ratio <- paste0(hr,'(',hrLow,'-',hrHigh,')')
  pVal <- ifelse(rt$pvalue<0.001, '<0.001', sprintf('%.3f', rt$pvalue))
  
  ## 输出图形
  
  pdf(file = forestFile, width = 8, height = 4.5)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc = 2),width = c(3,2.5))
  
  ## 绘制森林图左边的临床信息
  
  xlim = c(0,3)
  par(mar = c(4,2.5,2,1))
  plot(1, xlim = xlim, ylim = ylim, type = 'n', axes = F, xlab = '', ylab = '')
  text.cex = 0.8
  text(0, n:1, Variables, adj = 0, cex = text.cex)
  text(1.5-0.5*0.2, n:1, pVal, adj = 1, cex = text.cex); text(1.5-0.5*0.2, n+1, 'pvalue', cex = text.cex, font = 2, adj = 1)
  text(3.1, n:1, Hazard_ratio, adj = 1, cex = text.cex); text(3.1, n+1, 'Hazard ratio', cex = text.cex, font = 2, adj = 1)
  
  ## 绘制森林图
  
  par(mar = c(4,1,2,1), mgp = c(2,0.5,0))
  xlim = c(0, max(as.numeric(hrLow), as.numeric(hrHigh)))
  plot(1, xlim = xlim, ylim = ylim, type = 'n', ylab = '', xaxs = 'i', xlab = 'Hazard ratio')
  arrows(as.numeric(hrLow), n:1, as.numeric(hrHigh), n:1, angle = 90, code = 3, length = 0.05, col = 'darkblue', lwd = 2.5)
  abline(v = 1, col = 'black', lty = 2, lwd = 2)
  boxcolor = ifelse(as.numeric(hr) > 1, forestCol, forestCol)
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex = 1.5)
  axis(1)
  dev.off()}

bioForest(coxFile = 'uniCox.txt', forestFile = 'uniForest.pdf', forestCol = 'green')
# x轴过大,无法绘制。

# 八、LASSO -----------------------------------------------------------------
library('glmnet')
## x_all为基因表达矩阵，需为matrix格式，y_all为临床信息矩阵，需要包含生存时间和生存状态，行名为样本名
x_all <- subset(DEBRGs_FPKM_exp_survival,select = -c(OS, OS.time)) 
x_all <- x_all[,rownames(res_results_0.05)]
y_all <- subset(DEBRGs_FPKM_exp_survival, select = c(OS, OS.time))
cvfit <-  cv.glmnet(as.matrix(x_all),
                    Surv(y_all$OS.time,y_all$OS),nfold=10,
                    #10倍交叉验证，非必须限定条件，这篇文献有，其他文献大多没提nfold=10,
                    family = "cox") 
plot(cvfit, las =1) ## 7:7
fit <- glmnet(as.matrix(x_all), Surv(y_all$OS.time,y_all$OS), 
              family = "cox") 
cvfit$lambda.min
plot(fit, xvar = "lambda",label = TRUE, las=1)
coef.min = coef(cvfit, s = "lambda.min")  ## lambda.min & lambda.1se 取一个
active.min = which(coef.min@i != 0) ## 找出那些回归系数没有被惩罚为0的
lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i+1] ## 提取基因名称
coef.min

# 九、多因素Cox ----------------------------------------------------------------
# gene_list <- rownames(res_results_0.05)
gene_list <- lasso_geneids
cox_data <- as.formula(paste0('Surv(OS.time, OS)~',
                              paste(gene_list,
                                    sep = '',
                                    collapse = '+')))
cox_more <- coxph(cox_data,data = as.data.frame(DEBRGs_FPKM_exp_survival))

## check if variances are supported by PH hypothesis.
cox_zph <- cox.zph(cox_more)
cox_table <- cox_zph$table[-nrow(cox_zph$table),]
#remove variances not supported by ph hypothesis and perform the 2nd regression
cox_formula <- as.formula(paste0('Surv(OS.time, OS)~',
                                 paste(rownames(cox_table)[cox_table[,3]>0.05],
                                       sep = '',
                                       collapse = '+')))
cox_more_2 <- coxph(cox_formula,data = as.data.frame(DEBRGs_FPKM_exp_survival))
cox_more_2_Sum <- summary(cox_more_2)
cox_more_2_Tab <- data.frame()
cox_more_2_Tab <- cbind(
  HR = cox_more_2_Sum$conf.int[,'exp(coef)'],
  HR.95L = cox_more_2_Sum$conf.int[,'lower .95'],
  HR.95H = cox_more_2_Sum$conf.int[,'upper .95'],
  pvalue = cox_more_2_Sum$coefficients[,'Pr(>|z|)'])
cox_more_2_Tab <- cbind(id = row.names(cox_more_2_Tab), cox_more_2_Tab)
cox_more_2_Tab <- as.data.frame(cox_more_2_Tab)
cox_more_2_Tab <- cox_more_2_Tab[,-1]
cox_more_2_Tab$id <- rownames(cox_more_2_Tab)
cox_more_2_Tab <- cox_more_2_Tab[,c(5,1:4)] 
## check the co-linearity between samples检测是否共线性，vif < 2
cox_correlation <- cor(DEBRGs_FPKM_exp_survival[,rownames(cox_table)[cox_table[,3]>0.05]],
                       method = 'pearson')

library(GGally)
cox_corr <- ggpairs(DEBRGs_FPKM_exp_survival[,rownames(cox_table)[cox_table[,3]>0.05]],
                    axisLabels = 'show')+
  theme_bw()+
  theme(panel.background = element_rect(colour = 'black',
                                        size = 1,
                                        fill = 'white'),
        panel.grid = element_blank())
cox_corr
library('rms')
vif <- rms::vif(cox_more_2)
#some people said if the square root of VIF >2, they might be co-linear
vif
sqrt(vif) < 2#均未达到共线性

Hazard_ratio <- ggforest(model = cox_more_2, 
                         data = DEBRGs_FPKM_exp_survival,
                         main = 'Hazard ratio', 
                         fontsize = 1)
Hazard_ratio
## 输出多因素结果
## #提取HR，95%置信区间和p值

## 计算C-index
C_index <- cox_more_2$concordance['concordance']
if(C_index >= 0.9){
  print("High accuracy")
}else{
  if(C_index < 0.9 & C_index >= 0.7){
    print("Medium accuracy")
  }else{
    print("Low accuracy")
  }
}
sum.surv<-summary(cox_more_2)
c_index<-sum.surv$concordance
c_index

# 202111客户反馈,提供多因素cox森林图
# 进行多因素cox
gene_list <- lasso_geneids
cox_data <- as.formula(paste0('Surv(OS.time, OS)~',
                              paste(gene_list,
                                    sep = '',
                                    collapse = '+')))
multiCox <- coxph(cox_data,data = as.data.frame(DEBRGs_FPKM_exp_survival))
multiCoxSum <- summary(multiCox)
multiTab <- data.frame()
multiTab <- cbind(
  HR = multiCoxSum$conf.int[,'exp(coef)'],
  HR.95L = multiCoxSum$conf.int[,'lower .95'],
  HR.95H = multiCoxSum$conf.int[,'upper .95'],
  pvalue = multiCoxSum$coefficients[,'Pr(>|z|)'])
multiTab <- cbind(id = row.names(multiTab), multiTab)
multiTab <- as.data.frame(multiTab)
multiTab <- multiTab[,-1]
multiTab$id <- rownames(multiTab)
multiTab <- multiTab[,c(5,1:4)]
write.table(multiTab, file = 'multiCox.txt', sep = '\t', row.names = F, quote = F)
# 绘制森林图
bioForest = function(coxFile = null, forestFile = null, forestCol = null){
  rt <- read.table(coxFile, header = T, sep = '\t', check.names = F, row.names = 1)
  Variables <-  rownames(rt)
  hr <- sprintf('%.3f', rt$'HR')
  hrLow <- sprintf('%.3f', rt$'HR.95L')
  hrHigh <- sprintf('%.3f', rt$'HR.95H')
  Hazard_ratio <- paste0(hr,'(',hrLow,'-',hrHigh,')')
  pVal <- ifelse(rt$pvalue<0.001, '<0.001', sprintf('%.3f', rt$pvalue))
  
  ## 输出图形
  
  pdf(file = forestFile, width = 8, height = 4.5)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc = 2),width = c(3,2.5))
  
  ## 绘制森林图左边的临床信息
  
  xlim = c(0,3)
  par(mar = c(4,2.5,2,1))
  plot(1, xlim = xlim, ylim = ylim, type = 'n', axes = F, xlab = '', ylab = '')
  text.cex = 0.8
  text(0, n:1, Variables, adj = 0, cex = text.cex)
  text(1.5-0.5*0.2, n:1, pVal, adj = 1, cex = text.cex); text(1.5-0.5*0.2, n+1, 'pvalue', cex = text.cex, font = 2, adj = 1)
  text(3.1, n:1, Hazard_ratio, adj = 1, cex = text.cex); text(3.1, n+1, 'Hazard ratio', cex = text.cex, font = 2, adj = 1)
  
  ## 绘制森林图
  
  par(mar = c(4,1,2,1), mgp = c(2,0.5,0))
  xlim = c(0, max(as.numeric(hrLow), as.numeric(hrHigh)))
  plot(1, xlim = xlim, ylim = ylim, type = 'n', ylab = '', xaxs = 'i', xlab = 'Hazard ratio')
  arrows(as.numeric(hrLow), n:1, as.numeric(hrHigh), n:1, angle = 90, code = 3, length = 0.05, col = 'darkblue', lwd = 2.5)
  abline(v = 1, col = 'black', lty = 2, lwd = 2)
  boxcolor = ifelse(as.numeric(hr) > 1, forestCol, forestCol)
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex = 1.5)
  axis(1)
  dev.off()}
bioForest(coxFile = 'multiCox.txt', forestFile = 'multiForest.pdf', forestCol = 'red')

# 十、Kaplan–Meier生存分析 -------------------------------------------------------
riskScore = predict(cox_more_2,type="risk",newdata=DEBRGs_FPKM_exp_survival)
coxGene = rownames(as.data.frame(cox_more_2$coefficients))
coxGene = gsub("`","",coxGene)
outCol = c("OS","OS.time",coxGene)
risk_train = as.vector(ifelse(riskScore>median(riskScore),0,1))
risk_train <- as.data.frame(c(cbind(id=rownames(cbind(DEBRGs_FPKM_exp_survival[,outCol],
                                                      riskScore,
                                                      risk_train)),
                                    cbind(DEBRGs_FPKM_exp_survival[,outCol],
                                          riskScore,
                                          risk_train))))
## 将生存时间转换成年
risk_train$OS.time <- risk_train$OS.time / 365
rownames(risk_train) <- risk_train[,1]
risk_train <- risk_train[,-1]
## 生存曲线
kmfit <- survfit(Surv(OS.time, OS) ~ risk_train, data = risk_train)
Kaplan_Meier <- ggsurvplot(kmfit,
                           pval = TRUE, 
                           conf.int = F,
                           legend.labs=c("High risk","Low risk" ),
                           legend.title="Risk score",
                           risk.table = TRUE, 
                           risk.table.col = "strata", 
                           linetype = "strata", 
                           surv.median.line = "hv", 
                           ggtheme = theme_bw(), 
                           palette = c("#E7B800", "#2E9FDF"))+xlab("Time(years)")
Kaplan_Meier ## 5:6

# 十一、ROC曲线 ----------------------------------------------------------------
library(timeROC)
library(survival)
## 将生存时间除以365,由天变换成年
risk_roc <- risk_train
## 开始绘图
bioROC=function(inputFile=null, rocFile=null){
  ROC_rt=timeROC(T=risk_roc$OS.time, delta=risk_roc$OS,
                 marker=risk_roc$riskScore, cause=1,
                 weighting='aalen',
                 times=c(1,3,5), ROC=TRUE)
  pdf(file=rocFile,width=5,height=5)
  plot(ROC_rt,time=1,col='green',title=FALSE,lwd=2)
  plot(ROC_rt,time=3,col='blue',add=TRUE,title=FALSE,lwd=2)
  plot(ROC_rt,time=5,col='red',add=TRUE,title=FALSE,lwd=2)
  legend('bottomright',
         c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
           paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
           paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
         col=c("green","blue","red"),lwd=2,bty = 'n')
  dev.off()
}
bioROC(inputFile="trainRisk.txt",rocFile="train.ROC.pdf")

# 十二、风险曲线 -----------------------------------------------------------------
# 01-Train风险评分图 ----------------------------------------------------------------
## 数据准备
survStat_train <- risk_train 
## 对riskScore进行排序
risk_train <- risk_train[order(risk_train$riskScore),]
riskClass <- risk_train[,'risk_train'] ## 提取风险分类
## 绘制train组风险图
lowLength <- length(riskClass[riskClass == "1"])
highLength <- length(riskClass[riskClass == "0"])
line <- risk_train[,"riskScore"]
line[line>10] = 10
pdf(file="riskScoreTrain.pdf",width = 12,height = 5)
plot(line,
     type = "p",
     pch = 20,
     xlab = "Patients (increasing risk score)",
     ylab = "Risk score",
     col=c(rep("green",lowLength),
           rep("red",highLength)))
par(mfrow=c(1,1))
legend("topleft",
       c('High risk','Low risk'),
       pch = 20,
       ncol = 2,
       cex = 1,
       col = c("red","green"))
trainMedianScore=median(risk_train$riskScore)
abline(h=trainMedianScore,v=lowLength,lty=2)
dev.off()

## 绘制train组生存状态图
survStat_train <- risk_train
survStat_train <- survStat_train[order(survStat_train$riskScore),]
riskClass <- survStat_train[,'risk_train']
lowLength = length(riskClass[riskClass == "1"])
highLength = length(riskClass[riskClass == "0"])
color = as.vector(survStat_train$OS)
color[color == 1] = "red"
color[color == 0] = "green"
pdf(file="survStatTrain.pdf", width = 12, height = 5)
plot(survStat_train$OS.time,
     pch=19,
     xlab="Patients (increasing risk score)",
     ylab="Survival time (years)",
     col=color)
legend("topleft",
       c('Dead','Alive'),
       pch = 20,
       ncol = 2,
       cex = 1,
       col = c("red","green"))
abline(v=lowLength,lty=2)
dev.off()

## 绘制风险热图
library(pheatmap)
riskpheat_dat <- survStat_train
riskpheat_dat <- riskpheat_dat[order(riskpheat_dat$riskScore),] ## 对风险评分进行排序有小到大
## 提取基因表达矩阵
riskpheat_dat1 <- riskpheat_dat[c(3:(ncol(riskpheat_dat)-2))] 
## 对基因表达矩阵进行转置,行名是基因名,列名是样品名
riskpheat_dat1 <- t(riskpheat_dat1)
annotation <- data.frame(Type = riskpheat_dat[,ncol(riskpheat_dat)]) 
## 将1变成'Low',0变成'High'
annotation$Type <- ifelse(annotation$Type == 1, 'Low', 'High')
rownames(annotation) <- rownames(riskpheat_dat)
pdf(file="heatmapTrain.pdf",width = 15,height = 10)
pheatmap(riskpheat_dat1, 
         annotation_col = annotation,
         cluster_cols = FALSE, ## 样品不聚类
         show_rownames = T,
         show_colnames = F,
         fontsize_row = 11,
         fontsize_col=3,
         cellwidth = 2, cellheight = 70, ## 小方块的长和宽
         color = colorRampPalette(c("green", "black", "red"))(50))
dev.off()
# 十三、验证集验证 ----------------------------------------------------------------
# 01-数据预处理 ----------------------------------------------------------------
## 1.将验证集的表达量进行log转换
validation_set <- log2(validation_set+1)
## 转置
validation_set <- t(validation_set)
validation_set <- as.data.frame(validation_set)

## 将行名变成第一列,命名为'geo_accession'
validation_set$geo_accession <- rownames(validation_set)
validation_set <- validation_set[,c(13680,1:13679)]
rownames(validation_set) <- NULL

## 2.生存信息整理
### 创建一个数据框
geo_accession <- pheno$geo_accession
OS <- c(1,1,0,0,1,0,0,1,1,0,1,1,1,1,1,0,1,1,0,0,1,0,0,0,0,1,
        1,0,1,0,1,1,0,1,0,0,0,1,0,1,0,0,1,0,0,1,0,0,0,0,0,0,0)
OS.time <- c(810,630,1380,840,330,1110,1350,390,990,750,540,900,
             1050,120,810,780,540,5670,1080,3690,3300,1890,1800,
             1800,1800,300,1170,2850,2490,7380,750,1200,4290,330,
             3150,2340,2910,990,2310,1410,3600,2730,870,960,930,
             750,6570,5790,5520,5820,2820,2610,1800)
validation_set_survival <- data.frame(geo_accession,OS,OS.time)

## 3.基因表达矩阵和生存信息合并
validation_set <- merge(validation_set_survival,validation_set,by='geo_accession')
rownames(validation_set) <- validation_set[,1]
validation_set <- validation_set[,-1]
## validation_set是临床数据和表达量数据合并后的验证集

# 02-风险分析 -----------------------------------------------------------------
riskScore <- predict(cox_more_2, type = "risk", newdata = validation_set)
coxGene = rownames(as.data.frame(cox_more_2$coefficients))
coxGene = gsub("`","",coxGene)
outCol = c("OS","OS.time",coxGene)
risk_test = as.vector(ifelse(riskScore > median(riskScore),0,1)) ## 0代表高风险,1代表低风险
risk_test <- as.data.frame(c(cbind(id=rownames(cbind(validation_set[,outCol],
                                                     riskScore,
                                                     risk_test)),
                                   cbind(validation_set[,outCol],
                                         riskScore,
                                         risk_test))))
## 将risk中OS.time变成年
risk_test$OS.time <- risk_test$OS.time / 365

# 03-绘制风险曲线 ---------------------------------------------------------------
## 将第一列作为行名
rownames(risk_test) <- risk_test[,1]
risk_test <- risk_test[,-1]
write.table(risk_test,file = "risk_test.xls",sep = "\t",quote = F,row.names = T)

## (不绘制)
# 04-Kaplan-Meier生存分析 -----------------------------------------------------
kmfit <- survfit(Surv(OS.time, OS) ~ risk_test, data = risk_test)
Kaplan_Meier <- ggsurvplot(kmfit,
                           pval = TRUE, 
                           conf.int = F,
                           legend.labs=c("High risk","Low risk" ),
                           legend.title="Risk score",
                           risk.table = TRUE, 
                           risk.table.col = "strata", 
                           linetype = "strata", 
                           surv.median.line = "hv", 
                           ggtheme = theme_bw(), 
                           palette = c("#E7B800", "#2E9FDF"))+xlab("Time(years)")
Kaplan_Meier ## 5:6

# 05-ROC曲线 ----------------------------------------------------------------
library(timeROC)
library(survival)
## 将生存时间为年
risk_roc <- risk_test
## 开始绘图
bioROC=function(inputFile=null, rocFile=null){
  ROC_rt=timeROC(T=risk_roc$OS.time, delta=risk_roc$OS,
                 marker=risk_roc$riskScore, cause=1,
                 weighting='aalen',
                 times=c(1,3,5), ROC=TRUE)
  pdf(file=rocFile,width=5,height=5)
  plot(ROC_rt,time=1,col='green',title=FALSE,lwd=2)
  plot(ROC_rt,time=3,col='blue',add=TRUE,title=FALSE,lwd=2)
  plot(ROC_rt,time=5,col='red',add=TRUE,title=FALSE,lwd=2)
  legend('bottomright',
         c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
           paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
           paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
         col=c("green","blue","red"),lwd=2,bty = 'n')
  dev.off()
}
bioROC(inputFile="validationRisk.txt",rocFile="validation.ROC.pdf")

## 数据准备
survStat_test <- risk_test 
## 对riskScore进行排序
risk_test <- risk_test[order(risk_test$riskScore),]
riskClass <- risk_test[,'risk_test'] ## 提取风险分类
## 绘制test组风险图
lowLength <- length(riskClass[riskClass == "1"])
highLength <- length(riskClass[riskClass == "0"])
line <- risk_test[,"riskScore"]
#line[line>10] = 10
pdf(file="riskScoretest.pdf",width = 12,height = 5)
plot(line,
     type = "p",
     pch = 20,
     xlab = "Patients (increasing risk score)",
     ylab = "Risk score",
     col=c(rep("green",lowLength),
           rep("red",highLength)))
par(mfrow=c(1,1))
legend("topleft",
       c('High risk','Low risk'),
       pch = 20,
       ncol = 2,
       cex = 1,
       col = c("red","green"))
testMedianScore=median(risk_test$riskScore)
abline(h=testMedianScore,v=lowLength,lty=2)
dev.off()

## 绘制test组生存状态图
survStat_test <- risk_test
survStat_test <- survStat_test[order(survStat_test$riskScore),]
riskClass <- survStat_test[,'risk_test']
lowLength = length(riskClass[riskClass == "1"])
highLength = length(riskClass[riskClass == "0"])
color = as.vector(survStat_test$OS)
color[color == 1] = "red"
color[color == 0] = "green"
pdf(file="survStattest.pdf", width = 12, height = 5)
plot(survStat_test$OS.time,
     pch=19,
     xlab="Patients (increasing risk score)",
     ylab="Survival time (years)",
     col=color)
legend("topleft",
       c('Dead','Alive'),
       pch = 20,
       ncol = 2,
       cex = 1,
       col = c("red","green"))
abline(v=lowLength,lty=2)
dev.off()

## 绘制风险热图
library(pheatmap)
riskpheat_dat <- survStat_test
riskpheat_dat <- riskpheat_dat[order(riskpheat_dat$riskScore),] ## 对风险评分进行排序有小到大
## 提取基因表达矩阵
riskpheat_dat1 <- riskpheat_dat[c(3:(ncol(riskpheat_dat)-2))] 
## 对基因表达矩阵进行转置,行名是基因名,列名是样品名
riskpheat_dat1 <- t(riskpheat_dat1)
annotation <- data.frame(Type = riskpheat_dat[,ncol(riskpheat_dat)]) 
## 将1变成'Low',0变成'High'
annotation$Type <- ifelse(annotation$Type == 1, 'Low', 'High')
rownames(annotation) <- rownames(riskpheat_dat)
pdf(file="heatmaptest.pdf",width = 15,height = 10)
pheatmap(riskpheat_dat1, 
         annotation_col = annotation,
         cluster_cols = FALSE, ## 样品不聚类
         show_rownames = T,
         show_colnames = F,
         fontsize_row = 11,
         fontsize_col=3,
         cellwidth = 8, cellheight = 70, ## 小方块的长和宽
         color = colorRampPalette(c("green", "black", "red"))(50))
dev.off()

# 十三、肿瘤突变负荷差异分析 -----------------------------------------------------------
# library(maftools)
# library(TCGAmutations)
# library(TCGAbiolinks)
# (数据下载需要权限)

# 十四、免疫浸润细胞差异分析 ---------------------------------------------------------------
# 01-TCGA(TARGET-OS)_FPKM数据下载 ---------------------------------------------
## 同七、01
## 1.原始FPKM表达数据的下载
query_FPKM <- GDCquery(project = "TARGET-OS",
                       data.category = "Transcriptome Profiling",
                       data.type = "Gene Expression Quantification",
                       workflow.type = "HTSeq - FPKM") ## 使用FPKM数据
GDCdownload(query_FPKM)
TCGA_dat_FPKM <- GDCprepare(query_FPKM)
TCGA_dat_FPKM <- as.data.frame(assay(TCGA_dat_FPKM)) 
### (56602,88) 56602个ENSGxxx,88个样本

## 2.数据预处理
## 查看原始FPKM数据是否存在基因表达量全部为0的行
View(TCGA_dat_FPKM[which(rowSums(TCGA_dat_FPKM)==0),]) ## N=1654
## 表达量不为0的矩阵(去除基因表达量为0的行)
TCGA_dat_FPKM <- TCGA_dat_FPKM[which(rowSums(TCGA_dat_FPKM)>0),] ## 56602-1654=54948

## 将行名变成第一列,命名为Ensembl_ID
TCGA_dat_FPKM$Ensembl_ID <- rownames(TCGA_dat_FPKM)
rownames(TCGA_dat_FPKM) <- NULL
TCGA_dat_FPKM <- TCGA_dat_FPKM[,c(89,1:88)]

## ID转换:将TCGA_dat_FPKM按照Ensembl_ID列去重
aggr_TCGA_dat_FPKM <- aggregate(TCGA_dat_FPKM[,2:89],
                                by=list(ensembl_id=TCGA_dat_FPKM$Ensembl_ID),
                                mean) 
## 注:aggr_TCGA_dat_FPKM是Ensembl去0再去重之后的矩阵

## 将aggr_TCGA_dat_FPKM的Ensembl_ID转换成symbol
TCGA_dat_FPKM1 <- merge(s2e,aggr_TCGA_dat_FPKM,by='ensembl_id')
TCGA_dat_FPKM1 <- TCGA_dat_FPKM1[,-1] ## 去除ensembl_id列 (33878)
## symbol列去重
TCGA_dat_FPKM1 <- aggregate(TCGA_dat_FPKM1[,2:89],
                            by=list(symbol=TCGA_dat_FPKM1$symbol),
                            mean) 
rownames(TCGA_dat_FPKM1) <- TCGA_dat_FPKM1[,1]
TCGA_dat_FPKM1 <- TCGA_dat_FPKM1[,-1]
## TCGA_dat_FPKM1是去0去重之后的symbol名FPKM表达矩阵


# 02-TCGA_TIIC数据校正 --------------------------------------------------------
library('limma')
## 将表达矩阵重新命名为TCGA_TIIC
TCGA_TIIC <- TCGA_dat_FPKM1
TCGA_TIIC1 <- voom(TCGA_TIIC, plot = F, save.plot = F) ## 处理成类似芯片的格式
TCGA_TIIC1 <- as.data.frame(TCGA_TIIC1)
## 将行名作为第一列,命名为Genes
TCGA_TIIC1$Genes <- rownames(TCGA_TIIC1)
TCGA_TIIC1 <- TCGA_TIIC1[,c(89,1:88)]
rownames(TCGA_TIIC1) <- NULL
## 将文件读出
write.table(TCGA_TIIC1,
            file = "TCGA_TIIC1.xls",
            quote = F,
            sep = "\t",
            row.names = F)

# 03-免疫细胞组成(CIBERSORT) ---------------------------------------------------------------
library('parallelly')
library('e1071')
library('preprocessCore')
source("Cibersort.R") ## 将文件拷贝到当前路径下
Cibersort_result <- CIBERSORT('LM22.txt','TCGA_TIIC1.xls', 
                              perm = 1000, 
                              QN = F)
# 这个过程大约耗时1h
# 单独写一个脚本(ciber.R),在linux上运行
# load("Cibersort_result.rda")
Cibersort_result <- as.data.frame(Cibersort_result)
write.table(Cibersort_result,
            file = "CIBERSORT_result.xls",
            quote = F,
            sep = "\t",
            row.names = T)
## 去除Cibersort_result文件后3列
Cibersort_result <- Cibersort_result[,-c(23:25)]
Cibersort_result$id <- rownames(Cibersort_result)
Cibersort_result <- Cibersort_result[,c(23,1:22)]
rownames(Cibersort_result) <- NULL
## 将cibersort_result文件中样本名后的'-01R'去掉
Cibersort_result$id <- gsub(pattern = '-01R',
                            replacement = '',
                            x = Cibersort_result$id)
## 去除高低风险分组中不存在的3个样本
rownames(Cibersort_result) <- Cibersort_result[,1]
Cibersort_result <- Cibersort_result[,-1]
risk_sample <- rownames(risk_train)
Cibersort_result <- Cibersort_result[risk_sample,] ## 提取(N = 85)
## 宽数据转换成长数据
Cibersort_result1 <- t(Cibersort_result)
Cibersort_result1 <- as.data.frame(Cibersort_result1)
Cibersort_result1$TIICs <- rownames(Cibersort_result1)
Cibersort_result1 <- Cibersort_result1[,c(86,1:85)]
rownames(Cibersort_result1) <- NULL
Cibersort_result2 <- gather(Cibersort_result1,
                            key = sample,
                            value = fraction,
                            -c("TIICs"))
## 高低风险分组
high_risk <- risk_train[which(as.numeric(risk_train$risk_train) == 0),]## 提取risk=0(高风险组)的样本
high_risk_sample <- rownames(high_risk)
## CIBERSORT评分+高低风险分组
Cibersort_result2$group <- ifelse(Cibersort_result2$sample %in% high_risk_sample, "High-risk", "Low-risk")
colnames(Cibersort_result2)[4] <- "Group"
## 绘制箱线图
library(ggpubr)
library(RColorBrewer)
ggplot(Cibersort_result2,aes(TIICs, fraction, fill = Group)) + 
  geom_boxplot(outlier.shape = 21,color = "black")+
  theme_classic()+
  labs(x = "Cell Type", y = "Estimated Proportion")+
  theme(legend.position = "top")+ 
  theme(axis.text.x = element_text(angle=45,vjust = 1, hjust = 1))+
  scale_fill_manual(values = c(brewer.pal(7,'Set2')[c(4,5)]))+
  stat_compare_means(aes(group = Group,label = ..p.signif..),method = "kruskal.test")
# all_TIICs 5:8

## 提取差异显著的免疫细胞列表
sigdiff_TIICs <- c('B cells memory',
                   'Mast cells activated',
                   'Mast cells resting',
                   'Plasma cells',
                   'T cells CD4 memory activated')
Cibersort_result_sigdiff <- t(Cibersort_result)
Cibersort_result_sigdiff <- Cibersort_result_sigdiff[sigdiff_TIICs,]
Cibersort_result_sigdiff <- as.data.frame(Cibersort_result_sigdiff)
## 差异显著的免疫细胞:宽数据转换成长数据
Cibersort_result_sigdiff$TIICs <- rownames(Cibersort_result_sigdiff)
Cibersort_result_sigdiff <- Cibersort_result_sigdiff[,c(86,1:85)]
rownames(Cibersort_result_sigdiff) <- NULL
Cibersort_result_sigdiff <- gather(Cibersort_result_sigdiff,
                                   key = sample,
                                   value = fraction,
                                   -c("TIICs"))
## 添加分组
Cibersort_result_sigdiff$group <- ifelse(Cibersort_result_sigdiff$sample %in% high_risk_sample, "High-risk", "Low-risk")
colnames(Cibersort_result_sigdiff)[4] <- 'Group'
## 差异显著免疫细胞(N=5)绘图

# # 04-免疫细胞浸润热图 -------------------------------------------------------------
# ## 数据准备:
# pheat_TIICs_dat <- Ciber_vioplot_dat_Ordered ## 高风险组在前(N=42),低风险组(N=43)在后
# ## 将第一列变成行名
# rownames(pheat_TIICs_dat) <- pheat_TIICs_dat[,1]
# pheat_TIICs_dat <- pheat_TIICs_dat[,-1]
# ## 将数据框转置
# pheat_TIICs_dat <- t(pheat_TIICs_dat)
# pheat_TIICs_dat <- as.data.frame(pheat_TIICs_dat)
# df3 <- Ciber_vioplot_dat[,c(1,2)]
# ## 将0替换成High_risk, 1替换成Low_risk
# df3$risk <- ifelse(df3$risk == 0, 'High-risk', 'Low-risk')
# colnames(df3)[2] <- 'Group'
# rownames(df3) <- df3[,1]
# df3 <- df3[,-1,drop=FALSE] ## 注释分组
# 
# ## 绘制免疫细胞丰度热图
# library(pheatmap)
# library(RColorBrewer)
# bk = unique(c(seq(0,0.3, length=50)))
# pheat_TIICs_dat1 <- pheat_TIICs_dat
# TIICs_pheatplot <- pheatmap(pheat_TIICs_dat,
#                             bk = bk,
#                             cellwidth = 6, cellheight = 15, ## 小方块的长和框
#                             annotation_col = df3,
#                             cluster_row = F,
#                             cluster_col = F,
#                             show_rownames = T,
#                             show_colnames = F,
#                             border_color = NA,
#                             fontsize_row = 8,
#                             drop_levels = T,
#                             color = colorRampPalette(c("green", "black", "red"))(50),)
# TIICs_pheatplot ## 10:7


# 05-显著差异免疫浸润细胞 -----------------------------------------------------------
library(ggpubr)
library(RColorBrewer)
library(ggplot2)
sigdiff_TIICs <- ggplot(Cibersort_result_sigdiff,
                        aes(TIICs, fraction), 
                        fill = Group)+
  geom_boxplot(aes(colour = Group))+
  labs(x = NULL, y = "Fraction")+
  theme_bw()+
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust = 1))+
  scale_fill_manual(values = c(brewer.pal(7,'Set2')[c(4,5)]))+
  stat_compare_means(aes(group = Group,label = ..p.signif..), method = "kruskal.test")
sigdiff_TIICs  ## 10:6





# 20211129按照客户要求修改为小提琴图。
sigdiff_TIICs <- ggplot(Cibersort_result_sigdiff,
                        aes(x = TIICs,y = fraction,fill = Group))+
  geom_violin(alpha = 0.8,width = 1)+
  labs(x = NULL, y = "Fraction")+
  theme_bw()+
  theme(legend.position = "top") + 
  theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust = 1))+
  scale_fill_manual(values = c(brewer.pal(7,'Set2')[c(4,5)]))+
  stat_compare_means(aes(group = Group,label = ..p.signif..), method = "kruskal.test")
sigdiff_TIICs
# 10:6

# 十五、风险评分与临床指标的相关性分析 -------------------------------------------------------------
# 01-临床数据下载 ---------------------------------------------------------------
## UCSC Xena本地下载、上传,文件名:TARGET-OS.clinical.tsv
## 读取临床数据
library(tidyr)
clinical <- read.table('TARGET-OS.clinical.tsv',
                       header = T,
                       sep = '\t',
                       quote = '')
## 提取85个样本的临床数据
clinical1 <- clinical[which(clinical$sample_id %in% risk_sample),]
## 查看clinical1文件中sample_id列是否存在重复
duplic_list <- clinical1$sample_id
aa <- duplic_list[duplicated(clinical1$sample_id)]
dupli_expr <- clinical1[which(clinical1$sample_id %in% aa),]
## 发现,样本名TARGET-40-PASKZZ-01A重复,有两个
clinical1 <- clinical1[!duplicated(clinical1$sample_id),] ## 删除重复值

## 样本分组
clinical_sample <- clinical1$sample_id ## N=84
risk_group <- risk_train
risk_group <- risk_group[,c(1:2,6:7)]
risk_group <- risk_group[clinical_sample,]
risk_group$sample_id <- rownames(risk_group)
rownames(risk_group) <- NULL
risk_group <- risk_group[,c(5,1:4)]

# 02-gender ---------------------------------------------------------------
clinical_gender <- clinical1[,c(1,3)]
rownames(clinical_gender) <- NULL
## 添加分组信息(高低风险组)
clinical_gender$Group <- ifelse(clinical_gender$sample_id %in% high_risk_sample, "High-risk", "Low-risk")
risk1 <- risk_train[,c(1,6,7)]
## 将risk1第一列的样本名更改为'sample_id'
risk1$sample_id <- rownames(risk1)
## 将risk1和clinical_gender通过'sample_id'列合并
clinical_gender <- merge(risk1, clinical_gender, by = 'sample_id')
## 去除risk_traim列
clinical_gender <- clinical_gender[,-4]
## clinical_gender文件中Gender列存在空格需要去除
a <- gsub(' ', '', clinical_gender$Gender)
clinical_gender1 <- cbind(a, clinical_gender)
clinical_gender1 <- clinical_gender1[,-5]
colnames(clinical_gender1)[1] <- 'Gender'
clinical_gender1$Gender <- as.factor(clinical_gender1$Gender)
## 绘制箱线图
library(ggpubr)
library(ggsignif) 
#compaired <- list(c('Female', 'Male'))
gender_plot <- ggboxplot(clinical_gender1,
                         x = 'Gender', 
                         y = 'riskScore',
                         fill = 'Gender',
                         add = "jitter", size = 0.5)+
  stat_compare_means(method = "wilcox.test")+
  scale_fill_manual(values = c(brewer.pal(7,'Set2')[c(5,4)]))
theme_light()
gender_plot ##5:5

stat.test.gender<-clinical_gender1%>%
#  group_by(Group)%>%
  wilcox_test(riskScore ~ Gender)%>%
  adjust_pvalue(method = 'fdr')
# 20211129按照客户要求,绘制小提琴图
gender_plot <- ggplot(clinical_gender1,
                      aes(x = Gender,y = riskScore,fill = Gender))+
  geom_violin(alpha = 0.8,width = 1)+
  stat_compare_means(method = "wilcox.test")+
  scale_fill_manual(values = c(brewer.pal(7,'Set2')[c(5,4)]))
gender_plot

# 03-stage -----------------------------------------------------------------
clinical_stage <- clinical1[,c(1,27)]
rownames(clinical_stage) <- NULL
clinical_stage <- clinical_stage[-c(1:48,55:57,62,64:82),] ## 去除空值,仅剩13个样本
rownames(clinical_stage) <- NULL

## 设置分组
risk_group <- clinical_gender[,c(1,2,4)]
## 将clinical_stage1和risk1合并
clinical_stage <- merge(risk1, clinical_stage, by = 'sample_id')
##clinical_stage文件中Histologic.response列名更换为'stage'
colnames(clinical_stage)[5] <- 'Stage'
clinical_stage$Stage <- ifelse(clinical_stage$Stage == 'Stage 1/2', 'Stage I/II', 'Stage III/IV')

## 绘制箱线图
library(ggpubr)
## 将'Stage 1/2'移到前面
clinical_stage$Stage <- factor(clinical_stage$Stage,levels = c('Stage I/II', 'Stage III/IV'))
stage_plot <- ggboxplot(clinical_stage,
                        x = 'Stage', 
                        y = 'riskScore',
                        fill = 'Stage',
                        add = "jitter", size = 0.5)+
  stat_compare_means(method = 'wilcox.test')+
  scale_fill_manual(values = c(brewer.pal(7,'Set2')[c(5,4)]))
stage_plot 
# 20211129按照客户要求,绘制小提琴图
stage_plot <- ggplot(clinical_stage,
                     aes(x = Stage,y = riskScore,fill = Stage))+
  geom_violin(alpha = 0.8,width = 1)+
  stat_compare_means(method = "wilcox.test")+
  scale_fill_manual(values = c(brewer.pal(7,'Set2')[c(5,4)]))
stage_plot

# 04-age ------------------------------------------------------------------
clinical_age <- clinical1
clinical_age <- clinical1[,c(1,6)]
colnames(clinical_age)[2] <- 'Age' 
rownames(clinical_age) <- NULL
## 将年龄由days转换成years
clinical_age$Age <- clinical_age$Age / 365
clinical_age$Age <- floor(clinical_age$Age) ## 转换成整数
clinical_age <- merge(risk1, clinical_age, by = 'sample_id')
clinical_age$Age <- ifelse(clinical_age$Age >15, '>15', '≤ 15')
## 绘图
library(ggpubr)
age_plot <- ggboxplot(clinical_age,
                      x = 'Age', 
                      y = 'riskScore',
                      fill = 'Age',
                      add = "jitter", size = 0.5)+
  stat_compare_means(method = 'wilcox.test')+
  scale_fill_manual(values = c(brewer.pal(7,'Set2')[c(5,4)]))+
  xlab("Age (years)")
theme_light()
age_plot ## 5:5

# 20211129按照客户要求,绘制小提琴图
age_plot <- ggplot(clinical_age,
                   aes(x = Age,y = riskScore,fill = Age))+
  geom_violin(alpha = 0.8,width = 1)+
  stat_compare_means(method = "wilcox.test")+
  scale_fill_manual(values = c(brewer.pal(7,'Set2')[c(5,4)]))
age_plot

# 05-site -----------------------------------------------------------------
clinical_site <- clinical1[,c(1,16)]
rownames(clinical_site) <- NULL
colnames(clinical_site)[2] <- 'Site'
## 将'Arm/hand'替换成'Upper limb','Leg/Foot'和'Pelvis'替换成'Lower limb'
clinical_site$Site <- ifelse(clinical_site$Site == 'Arm/hand', 
                             'Upper limb', 
                             'Lower limb')
clinical_site <- merge(clinical_site, risk1, by = 'sample_id')
## 绘图
library(ggpubr)
library(ggsignif) 
site_plot <- ggboxplot(clinical_site,
                       x = 'Site', 
                       y = 'riskScore',
                       fill = 'Site',
                       add = "jitter", size = 0.5)+
  stat_compare_means(method = 'wilcox.test')+
  scale_fill_manual(values = c(brewer.pal(7,'Set2')[c(5,4)]))
theme_light()
site_plot ## 5:5

# 20211129按照客户要求,绘制小提琴图
site_plot <- ggplot(clinical_site,
                    aes(x = Site,y = riskScore,fill = Site))+
  geom_violin(alpha = 0.8,width = 1)+
  stat_compare_means(method = "wilcox.test")+
  scale_fill_manual(values = c(brewer.pal(7,'Set2')[c(5,4)]))
site_plot

# 十八、肿瘤纯度差异分析 -------------------------------------------------------------
library(utils)
library(estimate)

# 01-数据处理 -----------------------------------------------------------------
estimate_dat <- TCGA_dat_FPKM1 ## N = 88
estimate_dat <- t(estimate_dat)
estimate_dat <- as.data.frame(estimate_dat)
## 去除样本名中的'-01R'
estimate_dat$id <- rownames(estimate_dat)
estimate_dat <- estimate_dat[,c(33777,1:33776)]
rownames(estimate_dat) <- NULL
estimate_dat$id <- gsub(pattern = '-01R',
                        replacement = '',
                        x = estimate_dat$id)
rownames(estimate_dat) <- estimate_dat$id
estimate_dat <- estimate_dat[,-1]

## 去除高低风险组中不存在的3个样本
b <- rownames(risk_train)
estimate_dat <- estimate_dat[b,]
## 转置,将表达矩阵estimate文件调整为基因在行,样本在列
estimate_dat <- t(estimate_dat)
estimate_dat <- as.data.frame(estimate_dat)
write.table(estimate_dat, 
            'estimate_dat.txt', 
            col.names = T, 
            row.names = T, 
            quote = F, sep="\t")
# 02-ESTIMATE计算 -------------------------------------------------------------
filterCommonGenes(input.f = 'estimate_dat.txt', 
                  output.f = 'estimate_dat.gct', 
                  id = 'GeneSymbol')
estimateScore('estimate_dat.gct', 'estimate_purity.gct')
scores <- read.table('estimate_purity.gct', skip = 2, header = T)
rownames(scores) <- scores[,1]
scores <- t(scores[,3:ncol(scores)])
## 输出文件
write.table(scores, 
            'scores.xls', 
            col.names = T, 
            row.names = T, 
            quote = F, sep="\t")

# 03-差异分析 -----------------------------------------------------------------
# 03-1数据处理 ----------------------------------------------------------------
## 加入分组信息
scores1 <- scores
scores1 <- as.data.frame(scores1)
scores1$sample_id <- rownames(scores1)
scores1 <- scores1[,c(5,1:4)]
rownames(scores1) <- NULL
## 将'sample_id'列'.'替换成'-'
a <- gsub('\\.', '-', scores1$sample_id)
scores1 <- cbind(a,scores1)
scores1 <- scores1[,-2]
colnames(scores1)[1] <- 'sample_id'
a1 <- risk_train[,7,drop = F]
a1$sample_id <- rownames(a1)
rownames(a1) <- NULL
## 加入分组信息
scores1 <- merge(scores1,a1,by = 'sample_id')
colnames(scores1)[6] <- 'Group'
scores1$Group <- ifelse(scores1$Group == 1, 'Low-risk', 'High-risk')
rownames(scores1) <- scores1[,1]
scores1 <- scores1[,-1]

# 03-2 stormal+immune+estimateScore -----------------------------------------
## 数据处理
scores2 <- scores1
scores2 <- scores2[,-c(4:5)]
scores2 <- t(scores2)
scores2 <- as.data.frame(scores2)
scores2$Score <- rownames(scores2)
scores2 <- scores2[,c(86,1:85)]
rownames(scores2) <- NULL
library(tidyr)
scores3 <- gather(scores2,
                  key = sample,
                  value = fraction,
                  -c("Score"))
## 添加分组
scores3$Group <- ifelse(scores3$sample %in% high_risk_sample, "High-risk", "Low-risk")
colnames(scores3)[3] <- 'value'

## 绘制小提琴图
library(ggpubr)
library(ggplot2)
library(RColorBrewer)
p1 <- ggboxplot(scores3,
                x = 'Score', 
                y = 'value',
                fill = 'Group')+
  labs(x = NULL, y = "value") +
  theme_classic()+
  theme(legend.position = "top") + 
  stat_compare_means(aes(group = Group,label = ..p.signif..), method = 'wilcox.test')+
  scale_fill_manual(values = c(brewer.pal(7,'Set2')[c(4,5)]))
p1 ## 6:4.5

# 03-3肿瘤纯度 ----------------------------------------------------------------
scores4 <- scores1
scores4 <- scores4[,-c(1:3,5)]
scores4 <- as.data.frame(scores4)
rownames(scores4) <- rownames(scores1)
colnames(scores4)[1] <- 'TumorPurity'
scores4$sample_id <- rownames(scores4)
rownames(scores4) <- NULL
## 添加分组信息
scores4$Group <- ifelse(scores4$sample_id %in% high_risk_sample, "High-risk", "Low-risk")
rownames(scores4) <- scores4$sample_id
scores4 <- scores4[,-2]

## 绘制箱线图
TumorPurity_plot <- ggboxplot(scores4,
                              x = 'Group', 
                              y = 'TumorPurity',
                              fill = 'Group',
                              add = "jitter", size = 0.5)+
  stat_compare_means(method = 'wilcox.test')+
  scale_fill_manual(values = c(brewer.pal(7,'Set2')[c(4,5)]))+
  labs(x = NULL) 
theme_light()
TumorPurity_plot ## 5:5

# 十六、GSEA -----------------------------------------------------------------
# 01-样本分组 -----------------------------------------------------------------
group_median <- data.frame(sample = risk$id,
                           group = risk$risk)

# 一般logFC以control为对照，这样logFC>1就表示实验组表达大于对照组
# 写在前面的levels是对照
group_median[group_median==0] <- "High"
group_median[group_median==1] <- "Low"
group_median$group <- factor(group_median$group, levels = c("Low", "High"))
rownames(group_median) <- group_median$sample
# 02-表达矩阵 -----------------------------------------------------------------
GSEA_exp <- train_set
GSEA_exp <- t(GSEA_exp)
GSEA_exp <- as.data.frame(GSEA_exp)
GSEA_exp$sample_id <- rownames(GSEA_exp) ## 将行名变成一列,命名为'sample_id'
GSEA_exp <- GSEA_exp[,c(13697,1:13696)]
rownames(GSEA_exp) <- NULL
## 将GSEA_exp文件中sample_id列中的-01R替换为空
GSEA_exp$sample_id <- gsub(pattern = '-01R',
                           replacement = '',
                           x = GSEA_exp$sample_id)
rownames(GSEA_exp) <- GSEA_exp[,1]
GSEA_exp <- GSEA_exp[,-1]
## 提取高低风险组中的85个样本,形成新的表达矩阵
GSEA_exp <- GSEA_exp[risk_sample,]
GSEA_exp <- t(GSEA_exp)
GSEA_exp <- as.data.frame(GSEA_exp)
# 判断分组文件和表达矩阵中的样品是否一致
all(rownames(group_median) %in% colnames(GSEA_exp)) ## TRUE
# 判断分组文件和表达矩阵中的样品顺序是否一致
all(rownames(group_median) == colnames(GSEA_exp)) ## TRUE
## 去除小数
GSEA_exp1 <- as.data.frame(apply(GSEA_exp, 2, function(x) as.integer(x)))
rownames(GSEA_exp1) <- rownames(GSEA_exp)
countData <- GSEA_exp1 ## 最终的表达矩阵
condition <- factor(group_median$group)
colData <- data.frame(row.names=colnames(countData), condition)

# 03-计算logFC值 -------------------------------------------------------------
## DESeq2使用原始count值做差异分析
## 1.构建dds矩阵
library(DESeq2)
dds <- DESeqDataSetFromMatrix(as.matrix(countData), 
                              colData = colData,
                              design = ~ condition) 
## 2.对原始dds矩阵进行标准化
dds <- DESeq(dds)
## 3.两两比较
res <- results(dds, contrast = c("condition","Low","High"))
resOrdered <- res[order(res$pvalue),] # 按照P值排序
DEG <- as.data.frame(resOrdered)
DEG <- na.omit(DEG)

logFC_cutoff <- 1
DEG$change = as.factor(
  ifelse(DEG$padj < 0.05 & abs(DEG$log2FoldChange) > logFC_cutoff,
         ifelse(DEG$log2FoldChange > logFC_cutoff ,'UP','DOWN'),'NOT')
)
table(DEG$change)

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(patchwork)
DEG$symbol <- rownames(DEG)
df <- bitr(unique(DEG$symbol), 
           fromType = "SYMBOL",
           toType = c("ENTREZID"),
           OrgDb = org.Hs.eg.db)
DEG=merge(DEG, df, by.y="SYMBOL", by.x="symbol")
data_all_sort <- DEG %>%
  arrange(desc(log2FoldChange))
geneList <- data_all_sort$log2FoldChange
names(geneList) <- data_all_sort$ENTREZID
head(geneList)
kk2 <- gseKEGG(geneList = geneList,
               organism = 'hsa',
               minGSSize = 10,
               maxGSSize = 200,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none")

rownames(kk2@result,)[head(order(kk2@result$enrichmentScore))]
kk2 <- as.data.frame(kk2)

gsea_median <- gseaplot2(kk2,1:10,color="red",base_size = 15)
gsea_median ## 12:12

# 十七、独立预后分析 ------------------------------------------------------------------

# 01-数据准备 -----------------------------------------------------------------
indep_dat <- risk_train
## 提取有临床信息的84个样本
indep_dat <- indep_dat[clinical_sample,]
## 将行名变成第一列,命名为'sample_id'
indep_dat$sample_id <- rownames(indep_dat)
rownames(indep_dat) <- NULL
indep_dat <- indep_dat[,c(8, 1:7)]
## 去除'DLX2'、'TERT'、'EVX1'3列
indep_dat <- indep_dat[,-c(4:6)] ## 0代表高风险组
## 纳入年龄
indep_age <- clinical_age
indep_age <- indep_age[,-c(2:4)]
## indep_dat和indep_age合并
indep_dat <- merge(indep_dat, indep_age, by = 'sample_id')
## 纳入性别
indep_gender <- clinical_gender1
indep_gender <- indep_gender[,-c(3:5)]
## indep_dat和indep_gender合并
indep_dat <- merge(indep_dat, indep_gender, by = 'sample_id')
## 纳入肿瘤发生位点
indep_site <- clinical_site
indep_site <- indep_site[,-c(3:5)]
## indep_dat和indep_site合并
indep_dat <- merge(indep_dat, indep_site, by = 'sample_id')
## 将第一列变成行名
rownames(indep_dat) <- indep_dat[,1]
indep_dat <- indep_dat[,-1]
## 去除risk列(高低风险组)
indep_dat <- indep_dat[,-4]
## 将indep_dat1文件中'gender'列Male更改为1, Female更改为0
indep_dat1 <- indep_dat
indep_dat1$Gender <- ifelse(indep_dat1$Gender == 'Male', 1, 0)
## 将indep_dat1文件中'Site'列Lower limb更改为1, Upper limb更改为0
indep_dat1$Site <- ifelse(indep_dat1$Site == 'Lower limb', 1, 0)
## 将indep_dat1文件中'Age'列去除
indep_dat1 <- indep_dat1[,-4]
clinical_age1 <- clinical1[,c(1,6)]
colnames(clinical_age1)[2] <- 'Age' 
rownames(clinical_age1) <- NULL
clinical_age1$Age <- clinical_age1$Age / 365
clinical_age1$Age <- floor(clinical_age1$Age) ## 转换成整数
clinical_age1$Age <- ifelse(clinical_age1$Age >15, 1, 0) ## ＞15更改为1, ≤15更改为0
## 合并clinical_age1和indep_dat1
indep_dat1$sample_id <- rownames(indep_dat1)
indep_dat1 <- indep_dat1[,c(6,1:5)]
rownames(indep_dat1) <- NULL
indep_dat1 <- merge(indep_dat1, clinical_age1, by = 'sample_id')
## 将indep_dat1文件的第一列变成行名
rownames(indep_dat1) <- indep_dat1[,1]
indep_dat1 <- indep_dat1[,-1]
indep_dat1 <- indep_dat1[,c(1:3,6,4:5)]
# 02-单因素Cox ---------------------------------------------------------------
library("survival")
library("survminer")
uniTab <- data.frame()
for(i in colnames(indep_dat1[,3:ncol(indep_dat1)])){
  cox <- coxph(Surv(OS.time, OS) ~ indep_dat1[,i], data = indep_dat1)
  coxSummary = summary(cox)
  uniTab = rbind(uniTab,
                 cbind(id = i,
                       HR = coxSummary$conf.int[,'exp(coef)'],
                       HR.95L = coxSummary$conf.int[,'lower .95'],
                       HR.95H = coxSummary$conf.int[,'upper .95'],
                       pvalue = coxSummary$coefficients[,'Pr(>|z|)']))
}
rownames(uniTab) <- uniTab[,1]
uniTab$HR <- as.numeric(uniTab$HR)
uniTab$HR.95L <- as.numeric(uniTab$HR.95L)
uniTab$HR.95H <- as.numeric(uniTab$HR.95H)
uniTab$pvalue <- as.numeric(uniTab$pvalue)

write.table(uniTab, file = 'uniCox.txt', sep = '\t', row.names = F, quote = F)

# 03-绘制森林图(uniCox) ----------------------------------------------------------------
bioForest = function(coxFile = null, forestFile = null, forestCol = null){
  rt <- read.table(coxFile, header = T, sep = '\t', check.names = F, row.names = 1)
  Variables <-  rownames(rt)
  hr <- sprintf('%.3f', rt$'HR')
  hrLow <- sprintf('%.3f', rt$'HR.95L')
  hrHigh <- sprintf('%.3f', rt$'HR.95H')
  Hazard_ratio <- paste0(hr,'(',hrLow,'-',hrHigh,')')
  pVal <- ifelse(rt$pvalue<0.001, '<0.001', sprintf('%.3f', rt$pvalue))
  
  ## 输出图形
  
  pdf(file = forestFile, width = 8, height = 4.5)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc = 2),width = c(3,2.5))
  
  ## 绘制森林图左边的临床信息
  
  xlim = c(0,3)
  par(mar = c(4,2.5,2,1))
  plot(1, xlim = xlim, ylim = ylim, type = 'n', axes = F, xlab = '', ylab = '')
  text.cex = 0.8
  text(0, n:1, Variables, adj = 0, cex = text.cex)
  text(1.5-0.5*0.2, n:1, pVal, adj = 1, cex = text.cex); text(1.5-0.5*0.2, n+1, 'pvalue', cex = text.cex, font = 2, adj = 1)
  text(3.1, n:1, Hazard_ratio, adj = 1, cex = text.cex); text(3.1, n+1, 'Hazard ratio', cex = text.cex, font = 2, adj = 1)
  
  ## 绘制森林图
  
  par(mar = c(4,1,2,1), mgp = c(2,0.5,0))
  xlim = c(0, max(as.numeric(hrLow), as.numeric(hrHigh)))
  plot(1, xlim = xlim, ylim = ylim, type = 'n', ylab = '', xaxs = 'i', xlab = 'Hazard ratio')
  arrows(as.numeric(hrLow), n:1, as.numeric(hrHigh), n:1, angle = 90, code = 3, length = 0.05, col = 'darkblue', lwd = 2.5)
  abline(v = 1, col = 'black', lty = 2, lwd = 2)
  boxcolor = ifelse(as.numeric(hr) > 1, forestCol, forestCol)
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex = 1.5)
  axis(1)
  dev.off()}
bioForest(coxFile = 'uniCox.txt', forestFile = 'uniForest.pdf', forestCol = 'green')

# 04-多因素Cox ---------------------------------------------------------------
multiCox <- coxph(Surv(OS.time, OS) ~ riskScore + Age + Gender + Site, data = indep_dat1)
multiCoxSum <- summary(multiCox)
multiTab <- data.frame()
multiTab <- cbind(
  HR = multiCoxSum$conf.int[,'exp(coef)'],
  HR.95L = multiCoxSum$conf.int[,'lower .95'],
  HR.95H = multiCoxSum$conf.int[,'upper .95'],
  pvalue = multiCoxSum$coefficients[,'Pr(>|z|)'])
multiTab <- cbind(id = row.names(multiTab), multiTab)
multiTab <- as.data.frame(multiTab)
multiTab <- multiTab[,-1]
multiTab$id <- rownames(multiTab)
multiTab <- multiTab[,c(5,1:4)]
write.table(multiTab, file = 'multiCox.txt', sep = '\t', row.names = F, quote = F)

# 05-绘制森林图(multiCox) ------------------------------------------------------
bioForest = function(coxFile = null, forestFile = null, forestCol = null){
  rt <- read.table(coxFile, header = T, sep = '\t', check.names = F, row.names = 1)
  Variables <-  rownames(rt)
  hr <- sprintf('%.3f', rt$'HR')
  hrLow <- sprintf('%.3f', rt$'HR.95L')
  hrHigh <- sprintf('%.3f', rt$'HR.95H')
  Hazard_ratio <- paste0(hr,'(',hrLow,'-',hrHigh,')')
  pVal <- ifelse(rt$pvalue<0.001, '<0.001', sprintf('%.3f', rt$pvalue))
  
  ## 输出图形
  
  pdf(file = forestFile, width = 8, height = 4.5)
  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc = 2),width = c(3,2.5))
  
  ## 绘制森林图左边的临床信息
  
  xlim = c(0,3)
  par(mar = c(4,2.5,2,1))
  plot(1, xlim = xlim, ylim = ylim, type = 'n', axes = F, xlab = '', ylab = '')
  text.cex = 0.8
  text(0, n:1, Variables, adj = 0, cex = text.cex)
  text(1.5-0.5*0.2, n:1, pVal, adj = 1, cex = text.cex); text(1.5-0.5*0.2, n+1, 'pvalue', cex = text.cex, font = 2, adj = 1)
  text(3.1, n:1, Hazard_ratio, adj = 1, cex = text.cex); text(3.1, n+1, 'Hazard ratio', cex = text.cex, font = 2, adj = 1)
  
  ## 绘制森林图
  
  par(mar = c(4,1,2,1), mgp = c(2,0.5,0))
  xlim = c(0, max(as.numeric(hrLow), as.numeric(hrHigh)))
  plot(1, xlim = xlim, ylim = ylim, type = 'n', ylab = '', xaxs = 'i', xlab = 'Hazard ratio')
  arrows(as.numeric(hrLow), n:1, as.numeric(hrHigh), n:1, angle = 90, code = 3, length = 0.05, col = 'darkblue', lwd = 2.5)
  abline(v = 1, col = 'black', lty = 2, lwd = 2)
  boxcolor = ifelse(as.numeric(hr) > 1, forestCol, forestCol)
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex = 1.5)
  axis(1)
  dev.off()}

bioForest(coxFile = 'multiCox.txt', forestFile = 'multiForest.pdf', forestCol = 'red')

# 十八、列线图 ------------------------------------------------------------------
library(rms)
library(regplot)
library(survival)
library(survminer)

# 01-数据准备 -----------------------------------------------------------------
nomogram_dat <- indep_dat ## N=84
nomogram_dat <- nomogram_dat[,-c(4:5)]
## 1.添加高低风险分组评分
a <- clinical_gender1
a <- a[,-c(1,4)]
nomogram_dat <- merge(nomogram_dat, a, by = 'sample_id')
## 2.删除Age列,加入真实的年龄数值
nomogram_dat <- nomogram_dat[,-4] ## 删除Age列
a1 <- clinical1
a1 <- a1[,c(1,6)]
rownames(a1) <- NULL
colnames(a1)[2] <- 'Age'
## 将a1与nomogram_dat合并,即:添加'Age'列
nomogram_dat <- merge(nomogram_dat, a1, by = 'sample_id')
## 将'Age'列的天数转换成年,调整成整数
nomogram_dat$Age <- nomogram_dat$Age / 365
nomogram_dat$Age <- floor(nomogram_dat$Age) ## 转换成整数
## 3.纳入肿瘤发生位点
indep_site <- clinical_site
indep_site <- indep_site[,-c(3:6)]
## indep_dat和indep_site合并
nomogram_dat <- merge(nomogram_dat, indep_site, by = 'sample_id')
#将数值转化为真实属性后再因子化
## Gender:1 = male; 0 = female
## nomogram_dat$Gender <- as.factor(ifelse(nomogram_dat$Gender == 1,"Male","Female"))
## Site:1 = Lower limb; 0 = Upper limb
## nomogram_dat$Site <- as.factor(ifelse(nomogram_dat$Site == 1,"Lower limb","Upper limb"))
## 将第一列变成行名
rownames(nomogram_dat) <- nomogram_dat[,1]
nomogram_dat <- nomogram_dat[,-1]

# 02-构建Cox回归Nomogram预测模型 ----------------------------------------------------------------
dd <- datadist(nomogram_dat)
options(datadist="dd") 
## 构建COX比例风险模型
coxNomo_model <- psm(Surv(OS.time, OS) ~ riskScore + Age + Gender + Site, 
                     data = nomogram_dat, 
                     dist = 'lognormal')
med <- Quantile(coxNomo_model) ## 计算中位生存时间
surv <- Survival(coxNomo_model) ## 构建生存概率函数

## 绘制COX回归生存概率Nomogram图
## 注意:nomogram_dat数据的OS.time是以'天'为单位
nom <- nomogram(coxNomo_model, 
                fun=list(function(x) surv(365, x),
                         function(x) surv(1095, x),
                         function(x) surv(1825, x)),
                funlabel=c("1-year Survival Probability",
                           "3-year Survival Probability",
                           "5-year Survival Probability"))
plot(nom, xfrac=.6) ## 10:7.5

## 评价COX回归的预测效果
## 第一步 计算c-index
rcorrcens(Surv(OS.time, OS) ~ predict(coxNomo_model), data = nomogram_dat)
# 03-多指标ROC曲线 -------------------------------------------------------------
library(survivalROC)
# 03-1数据准备 ----------------------------------------------------------------
multiROC_dat <- nomogram_dat
## 生存时间由天变成年
multiROC_dat$OS.time <- multiROC_dat$OS.time / 365


# 03-2绘制风险值ROC曲线 ----------------------------------------------------------
rocCol <- rainbow(ncol(multiROC_dat)-2) ## 彩虹色带取色
aucText=c()
pdf(file="multiROC.pdf", width=6, height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc <- survivalROC(Stime = multiROC_dat$OS.time, 
                   status = multiROC_dat$OS, 
                   marker = multiROC_dat$riskScore, 
                   predict.time = 1, 
                   method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1), col = rocCol[1], 
     xlab="False positive rate", ylab="True positive rate",
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
aucText <- c(aucText, paste0("riskScore"," (AUC=",sprintf("%.3f",roc$AUC),")"))
abline(0,1)

# 03-3绘制其它指标的ROC曲线 --------------------------------------------------------
j=1
for(i in colnames(multiROC_dat[,4:(ncol(multiROC_dat))])){
  roc = survivalROC(Stime = multiROC_dat$OS.time, 
                    status = multiROC_dat$OS, 
                    marker = multiROC_dat[,i], 
                    predict.time =1, 
                    method="KM")
  j=j+1
  aucText = c(aucText, paste0(i," (AUC=",sprintf("%.3f",roc$AUC),")"))
  lines(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1), col=rocCol[j],lwd = 2)
}
legend("bottomright", aucText,lwd=2,bty="n",col=rocCol)
dev.off()

# 04-列线图ROC曲线 -----------------------------------------------------------------
library(survival)
library(survminer)
library(timeROC)
library(rms)
library(regplot)
multiROC_dat$Gender <- ifelse(multiROC_dat$Gender == 'Male', 1, 0)
multiROC_dat$Site <- ifelse(multiROC_dat$Site == 'Lower limb', 1, 0)
res.cox=coxph(Surv(OS.time, OS) ~ riskScore+Age+Gender+Site, data = multiROC_dat)
nomoRisk=predict(res.cox, data=multiROC_dat, type="risk")
multiROC_dat$nomoRisk=nomoRisk
ROC_rt=timeROC(T=multiROC_dat$OS.time, delta=multiROC_dat$OS,
               marker=multiROC_dat$nomoRisk, cause=1,
               weighting='aalen',
               times=c(1,3,5), ROC=TRUE)
pdf(file="ROC.pdf", width=5, height=5)
plot(ROC_rt,time=1,col='green',title=FALSE,lwd=2)
plot(ROC_rt,time=3,col='blue',add=TRUE,title=FALSE,lwd=2)
plot(ROC_rt,time=5,col='red',add=TRUE,title=FALSE,lwd=2)
legend('bottomright',
       c(paste0('AUC at 1 years: ',sprintf("%.03f",ROC_rt$AUC[1])),
         paste0('AUC at 3 years: ',sprintf("%.03f",ROC_rt$AUC[2])),
         paste0('AUC at 5 years: ',sprintf("%.03f",ROC_rt$AUC[3]))),
       col=c("green","blue","red"),lwd=2,bty = 'n')
dev.off()

##
pdf(file="calibration.pdf", width=5, height=5)
multiROC_dat$OS.time <- multiROC_dat$OS.time*365
## 1年校准曲线
f <- cph(Surv(OS.time, OS) ~ nomoRisk, x=T, y=T, surv=T, data=multiROC_dat, time.inc=365)
cal <- calibrate(f, cmethod="KM", method="boot", u=365, m=(nrow(multiROC_dat)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1),
     xlab="Nomogram-predicted OS (%)", ylab="Observed OS (%)", lwd=1.5, col="green", sub=F)
## 3年校准曲线
f <- cph(Surv(OS.time, OS) ~ nomoRisk, x=T, y=T, surv=T, data=multiROC_dat, time.inc=365*3)
cal <- calibrate(f, cmethod="KM", method="boot", u=365*3, m=(nrow(multiROC_dat)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", lwd=1.5, col="blue", sub=F, add=T)
## 5年校准曲线
f <- cph(Surv(OS.time, OS) ~ nomoRisk, x=T, y=T, surv=T, data=multiROC_dat, time.inc=365*5)
cal <- calibrate(f, cmethod="KM", method="boot", u=365*5, m=(nrow(multiROC_dat)/3), B=1000)
plot(cal, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="",  lwd=1.5, col="red", sub=F, add=T)
legend('bottomright', c('1-year', '3-year', '5-year'),
       col=c("green","blue","red"), lwd=1.5, bty = 'n')
dev.off()

# 十九、GSVA -----------------------------------------------------------------
library(GSVA)
library(GSEABase) ## 用于读取gmt格式基因集文件
library(limma)
# 01-表达矩阵 -----------------------------------------------------------------
GOBP_exp <- train_set
GOBP_exp <- t(GOBP_exp)
GOBP_exp <- as.data.frame(GOBP_exp)
GOBP_exp <- GOBP_exp[804:891,] ## 提取TCGA表达量数据
## 将行名作为第一列
GOBP_exp$id <- rownames(GOBP_exp)
rownames(GOBP_exp) <- NULL
GOBP_exp <- GOBP_exp[,c(13697,1:13696)]
## 将GOBP_exp文件id列中的-01R替换为空
GOBP_exp$id <- gsub(pattern = '-01R',
                    replacement = '',
                    x = GOBP_exp$id)
rownames(GOBP_exp) <- GOBP_exp[,1]
GOBP_exp <- GOBP_exp[,-1]
## 提取高低风险组的85个样本
a <- risk$id
GOBP_exp <- GOBP_exp[a,]
GOBP_exp <- t(GOBP_exp) ## 数据框转置,行为gene symbol,列为样本
GOBP_exp <- as.data.frame(GOBP_exp)
## 将低风险组的样本排到前面
GOBP_exp <- t(GOBP_exp)
GOBP_exp <- as.data.frame(GOBP_exp)
GOBP_exp$sample_id <- rownames(GOBP_exp)
GOBP_exp <- GOBP_exp[,c(13697,1:13696)]
rownames(GOBP_exp) <- NULL
GOBP_exp <- merge(a1, GOBP_exp, by = 'sample_id')
a2 <- risk
rownames(a2) <- a2$id
a2 <- a2[,-1]
a2 <- a2[,-c(1:5)]
a2$sample_id <- rownames(a2)
a2 <- a2[,-2]
rownames(a2) <- NULL
GOBP_exp <- merge(a2, GOBP_exp, by = 'sample_id')
## 按照风险评分高低进行排序,前43个为低风险,后42个为高风险
GOBP_exp <- GOBP_exp[order(GOBP_exp$riskScore),] 
rownames(GOBP_exp) <- NULL
GOBP_exp <- GOBP_exp[,-c(2:3)]
rownames(GOBP_exp) <- GOBP_exp[,1]
GOBP_exp <- GOBP_exp[,-1]
GOBP_exp <- t(GOBP_exp)
GOBP_exp <- as.data.frame(GOBP_exp)
exp <- GOBP_exp

# 01-1 GO-BP --------------------------------------------------------------
## 从GSEA | MSigDB数据库下载数据集(C5-BP:subset of GO)
## 文件名:c5.go.bp.v7.4.symbols.gmt
## 1.读取参考基因集文件
GOBP_ref <- getGmt('c5.go.bp.v7.4.symbols.gmt')
## 2.使用GSVA方法进行分析,获得GSVA得分
es_GOBP <- gsva(as.matrix(GOBP_exp), GOBP_ref,
                min.sz=10, max.sz=500, verbose=TRUE)
es_GOBP <- as.data.frame(es_GOBP)

## 得到gsva计算的数值后再用limma包做差异分析得到差异的pathway
## 3.设置分组
grouP <- c(rep("Low", 43), rep("High", 42)) %>% as.factor()
desigN <- model.matrix(~ grouP + 0)
rownames(desigN) <- colnames(GOBP_exp)
desigN
# 用desigN的列名，设定High-risk组比Low-risk组
comparE <- makeContrasts(grouPHigh-grouPLow,levels = desigN)
## 进行差异分析
fit <- lmFit(es_GOBP,desigN)
fit2 <- contrasts.fit(fit, comparE)
fit3 <- eBayes(fit2)
allGeneSets <- topTable(fit3, coef = 1, number = Inf)
## 画火山图
logFCcutoff <- 0
allGeneSets$change = as.factor(
  ifelse(allGeneSets$P.Value < 0.05 & abs(allGeneSets$logFC) > logFCcutoff,
         ifelse(allGeneSets$logFC > logFCcutoff,'UP','DOWN'),'NOT')
)

DEGeneSets <- subset(allGeneSets,
                     allGeneSets$P.Value < 0.05 & abs(allGeneSets$logFC) > 0)

write.table(allGeneSets,
            file = "GOBP.xls",
            quote = F,
            sep = "\t",
            row.names = T)
## 画图
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(scales)
library(ggrepel) 
volcano <- allGeneSets
volcano <- allGeneSets[order(volcano$P.Value,decreasing = F),]
volcano$anno_name <- rownames(volcano)
volcano$anno_name[11:nrow(volcano)] <- NA
colnames(volcano)[1] <- 'log2FoldChange'
colnames(volcano)[4] <- 'P.Value'
volcano_plot <- ggplot(data = volcano, 
                       aes(x = -log10(P.Value),
                           y = log2FoldChange,
                           colour = change))+
  geom_point(alpha = 0.7)+
  scale_color_manual(values=c('blue','darkgray','red'))+
  geom_text_repel(aes(label=anno_name), show.legend = F, segment.colour = 'black', size = 2, hjust = 10)+
  geom_vline(xintercept = -log10(0.05),
             lty = 4,
             col = "darkgray",
             lwd = 0.6)+
  xlim(0,10)+
  ylim(-0.5,0.5)+
  theme_classic(base_line_size = 0.5)+
  theme(axis.title.x = element_text(size = 12, 
                                    color = "black",
                                    face = "bold"),
        axis.title.y = element_text(size = 12,
                                    color = "black",
                                    face = "bold", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        legend.title = element_blank(),
        legend.text = element_text(color="black", # 设置图例标签文字
                                   size = 8, 
                                   face = "bold"),
        axis.text.x = element_text(size = 12,  # 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0), #角度
        axis.text.y = element_text(size = 12,  
                                   color = "black",
                                   face = "bold", 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0))
volcano_plot # 6:4.5(GOBP)

# 01-2 GO-CC ---------------------------------------------------------------
## 1.读取参考基因集文件c5.go.cc.v7.4.symbols.gmt
GOCC_ref <- getGmt('c5.go.cc.v7.4.symbols.gmt')
es_GOCC <- gsva(as.matrix(exp), GOCC_ref,
                min.sz=10, max.sz=500, verbose=TRUE)
es_GOCC <- as.data.frame(es_GOCC)
write.table(es_GOCC,
            file = "es_GOCC.xls",
            quote = F,
            sep = "\t",
            row.names = T)
fit <- lmFit(es_GOCC,desigN)
fit2 <- contrasts.fit(fit,comparE)
fit3 <- eBayes(fit2)
allGeneSets <- topTable(fit3,coef = 1, number = Inf)
## 画火山图
logFCcutoff <- 0
allGeneSets$change = as.factor(
  ifelse(allGeneSets$P.Value < 0.05 & abs(allGeneSets$logFC) > logFCcutoff,
         ifelse(allGeneSets$logFC > logFCcutoff,'UP','DOWN'),'NOT')
)

DEGeneSets <- subset(allGeneSets,
                     allGeneSets$P.Value < 0.05 & abs(allGeneSets$logFC) > 0)
write.table(allGeneSets,
            file = "GOCC.xls",
            quote = F,
            sep = "\t",
            row.names = T)
## 画图
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(scales)
library(ggrepel) 
volcano <- allGeneSets
volcano <- allGeneSets[order(volcano$P.Value,decreasing = F),]
volcano$anno_name <- rownames(volcano)
volcano$anno_name[11:nrow(volcano)] <- NA
colnames(volcano)[1] <- 'log2FoldChange'
colnames(volcano)[4] <- 'P.Value'
volcano_plot <- ggplot(data = volcano, 
                       aes(x = -log10(P.Value),
                           y = log2FoldChange,
                           colour = change))+
  geom_point(alpha = 0.7)+
  scale_color_manual(values=c('blue','darkgray','red'))+
  geom_text_repel(aes(label=anno_name), show.legend = F, segment.colour = 'black', size = 2, hjust = 10)+
  geom_vline(xintercept = -log10(0.05),
             lty = 4,
             col = "darkgray",
             lwd = 0.6)+
  xlim(0,8)+
  ylim(-0.5,0.5)+
  theme_classic(base_line_size = 0.5)+
  theme(axis.title.x = element_text(size = 12, 
                                    color = "black",
                                    face = "bold"),
        axis.title.y = element_text(size = 12,
                                    color = "black",
                                    face = "bold", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        legend.title = element_blank(),
        legend.text = element_text(color="black", # 设置图例标签文字
                                   size = 8, 
                                   face = "bold"),
        axis.text.x = element_text(size = 12,  # 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0), #角度
        axis.text.y = element_text(size = 12,  
                                   color = "black",
                                   face = "bold", 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0))
volcano_plot # 6:4.5(GOCC)

# 01-3 GO-MF --------------------------------------------------------------
GOMF_ref <- getGmt('c5.go.mf.v7.4.symbols.gmt')
es_GOMF <- gsva(as.matrix(exp), GOMF_ref,
                min.sz=10, max.sz=500, verbose=TRUE)
es_GOMF <- as.data.frame(es_GOMF)
fit <- lmFit(es_GOMF,desigN)
fit2 <- contrasts.fit(fit,comparE)
fit3 <- eBayes(fit2)
allGeneSets <- topTable(fit3,coef = 1, number = Inf)
## 画火山图
logFCcutoff <- 0
allGeneSets$change = as.factor(
  ifelse(allGeneSets$P.Value < 0.05 & abs(allGeneSets$logFC) > logFCcutoff,
         ifelse(allGeneSets$logFC > logFCcutoff,'UP','DOWN'),'NOT')
)

DEGeneSets <- subset(allGeneSets,
                     allGeneSets$P.Value < 0.05 & abs(allGeneSets$logFC) > 0)
write.table(allGeneSets,
            file = "GOMF.xls",
            quote = F,
            sep = "\t",
            row.names = T)
## 画图
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(scales)
library(ggrepel) 
volcano <- allGeneSets
volcano <- allGeneSets[order(volcano$P.Value,decreasing = F),]
volcano$anno_name <- rownames(volcano)
volcano$anno_name[11:nrow(volcano)] <- NA
colnames(volcano)[1] <- 'log2FoldChange'
colnames(volcano)[4] <- 'P.Value'
volcano_plot <- ggplot(data = volcano, 
                       aes(x = -log10(P.Value),
                           y = log2FoldChange,
                           colour = change))+
  geom_point(alpha = 0.7)+
  scale_color_manual(values=c('blue','darkgray','red'))+
  geom_text_repel(aes(label=anno_name), show.legend = F, segment.colour = 'black', size = 2, hjust = 10)+
  geom_vline(xintercept = -log10(0.05),
             lty = 4,
             col = "darkgray",
             lwd = 0.6)+
  xlim(0,8)+
  ylim(-0.5,0.5)+
  theme_classic(base_line_size = 0.5)+
  theme(axis.title.x = element_text(size = 12, 
                                    color = "black",
                                    face = "bold"),
        axis.title.y = element_text(size = 12,
                                    color = "black",
                                    face = "bold", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        legend.title = element_blank(),
        legend.text = element_text(color="black", # 设置图例标签文字
                                   size = 8, 
                                   face = "bold"),
        axis.text.x = element_text(size = 12,  # 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0), #角度
        axis.text.y = element_text(size = 12,  
                                   color = "black",
                                   face = "bold", 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0))
volcano_plot # 6:4.5(GOMF)

# 01-4 KEGG ---------------------------------------------------------------
## 1.读取参考基因集文件c5.go.cc.v7.4.symbols.gmt
KEGG_ref <- getGmt('c2.cp.kegg.v7.4.symbols.gmt')
es_KEGG <- gsva(as.matrix(exp), KEGG_ref,
                min.sz=10, max.sz=500, verbose=TRUE)
es_KEGG <- as.data.frame(es_KEGG)
fit <- lmFit(es_KEGG,desigN)
fit2 <- contrasts.fit(fit,comparE)
fit3 <- eBayes(fit2)
allGeneSets <- topTable(fit3,coef = 1, number = Inf)
## 画火山图
logFCcutoff <- 0
allGeneSets$change = as.factor(
  ifelse(allGeneSets$P.Value < 0.05 & abs(allGeneSets$logFC) > logFCcutoff,
         ifelse(allGeneSets$logFC > logFCcutoff,'UP','DOWN'),'NOT')
)

DEGeneSets <- subset(allGeneSets,
                     allGeneSets$P.Value < 0.05 & abs(allGeneSets$logFC) > 0)

write.table(allGeneSets,
            file = "KEGG.xls",
            quote = F,
            sep = "\t",
            row.names = T)
## 画图
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(scales)
library(ggrepel) 
volcano <- allGeneSets
volcano <- allGeneSets[order(volcano$P.Value,decreasing = F),]
volcano$anno_name <- rownames(volcano)
volcano$anno_name[11:nrow(volcano)] <- NA
colnames(volcano)[1] <- 'log2FoldChange'
colnames(volcano)[4] <- 'P.Value'
volcano_plot <- ggplot(data = volcano, 
                       aes(x = -log10(P.Value),
                           y = log2FoldChange,
                           colour = change))+
  geom_point(alpha = 0.7)+
  xlim(0,8)+
  ylim(-0.5,0.5)+
  scale_color_manual(values=c('blue','darkgray','red'))+
  geom_text_repel(aes(label=anno_name), show.legend = F, segment.colour = 'black', size = 2, hjust = 10)+
  geom_vline(xintercept = -log10(0.05),
             lty = 4,
             col = "darkgray",
             lwd = 0.6)+
  theme_classic(base_line_size = 0.5)+
  theme(axis.title.x = element_text(size = 12, 
                                    color = "black",
                                    face = "bold"),
        axis.title.y = element_text(size = 12,
                                    color = "black",
                                    face = "bold", 
                                    vjust = 1.9, 
                                    hjust = 0.5, 
                                    angle = 90),
        legend.title = element_blank(),
        legend.text = element_text(color="black", # 设置图例标签文字
                                   size = 8, 
                                   face = "bold"),
        axis.text.x = element_text(size = 12,  # 修改X轴上字体大小，
                                   color = "black", # 颜色
                                   face = "bold", #  face取值：plain普通，bold加粗，italic斜体，bold.italic斜体加粗
                                   vjust = 0.5, # 位置
                                   hjust = 0.5, 
                                   angle = 0), #角度
        axis.text.y = element_text(size = 12,  
                                   color = "black",
                                   face = "bold", 
                                   vjust = 0.5, 
                                   hjust = 0.5, 
                                   angle = 0))
volcano_plot # 6:4.5(KEGG)
