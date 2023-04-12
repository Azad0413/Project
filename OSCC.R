rm(list = ls())
# 01 获取数据集 
setwd("/data/nas1/luchunlin/project/BJTC-228")
if (! dir.exists("./00_rawdata")){
  dir.create("./00_rawdata")
}
setwd("./00_rawdata")
# 获取头颈鳞状细胞癌数据
library(TCGAbiolinks)
# 获取TCGA中最新的不同癌种的项目号
getGDCprojects()$project_id
query<-GDCquery(project = "TCGA-HNSC",
                data.category = "Transcriptome Profiling",
                data.type = "Gene Expression Quantification",
                workflow.type = "HTSeq - Counts")
# 546 个samples
samplesDown <- getResults(query = query,
                          cols = c("cases"))
# 500个Tumor samples
dataSmTP <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "TP")
# 44个Normal samples
dataSmNT <- TCGAquery_SampleTypes(barcode = samplesDown,
                                  typesample = "NT")
group_list<-data.frame(sample=c(dataSmNT,dataSmTP),
                       group=c(rep("Normal",44),rep("Tumor",500)))
GDCdownload(query)
## GDC prepare将前面的GDCquery的结果准备成R语言可处理的SE（summarizedExperiment）文件
TCGA_RNASeq<- GDCprepare(query = query,
                         save = TRUE,
                         save.filename = "TCGA_query.Rdata")
## 数据预处理，去除dataPrep中的异常值
dataPrep<-TCGAanalyze_Preprocessing(object = TCGA_RNASeq,
                                    cor.cut = 0.6,
                                    datatype = "HTSeq - Counts")
dataPrep <- subset(dataPrep, select = c(dataSmNT, dataSmTP))
# 将基因序号更换为基因名称
rownames(dataPrep) <- rowData(TCGA_RNASeq)[rownames(dataPrep),]$external_gene_name
dim(dataPrep)
# 下载临床信息
clinical_data <- GDCquery_clinic(project = "TCGA-HNSC",
                                 type = "clinical")
write.table(clinical_data,
            file = "clinical.xls",
            quote = F,
            sep = '\t',
            row.names = T)

# 将22个m6A调节因子的矩阵筛选出来
m6A_factor<-c("ALKBH5","FTO","CBLL1","HNRNPA2B1","HNRNPC","IGF2BP1","IGF2BP2","IGF2BP3",
              "METTL14","METTL16","METTL3","RBM15","RBM15B","VIRMA","WTAP","YTHDC1","YTHDC2","YTHDF1",
              "YTHDF2","YTHDF3","ZC3H13","ZCCHC4")
m6A_factor<-as.data.frame(m6A_factor)
colnames(m6A_factor)<-"symbol"
ids=m6A_factor[m6A_factor$symbol%in%rownames(dataPrep),]           # 判断是否匹配
dataPrep1<-dataPrep
dataPrep1=dataPrep1[m6A_factor$symbol,]                            # 取表达矩阵中匹配的部分，形成22个基因的小矩阵
write.table(dataPrep1,
            file = "dataPrep.xls",
            quote = F,
            sep = '\t',
            row.names = T)
#筛选到需要的样本 310个OSCC样本和44个正常样本,形成最终的表达矩阵
pd<-read_xlsx("/data/nas1/luchunlin/project/BJTC-228/00_rawdata/clinical_data.xlsx")
pd<-as.data.frame(pd)
rownames(pd)=pd[,1]
pd=pd[,-1]
dataPrep_final<-read_xlsx("/data/nas1/luchunlin/project/BJTC-228/00_rawdata/dataPrep_final.xlsx")
dataPrep_final<-as.data.frame(dataPrep_final)
rownames(dataPrep_final)=dataPrep_final[,1]
dataPrep_final=dataPrep_final[,-1]


# 样本过滤
###Tumor purity filtering
###vector containing all TCGA barcodes that hhave 60% tumor purity or more
# TCGAtumor_purity使用来自5种方法的5个估计值作为阈值对TCGA样本进行过滤
# 这5个值是estimate, absolute, lump, ihc, cpe，这里设置cpe=0.6，cpe是派生的共识度量，是将所有方法的标准含量归一化后的均值纯度水平。
purityDATA<-TCGAtumor_purity(colnames(TCGA_RNASeq), 0, 0, 0, 0, 0.6)  # 筛选肿瘤纯度大于60%的样本
# filtered 为被过滤的数据， pure_barcodes是我们要的数据
Purity.HNSC<-purityDATA$pure_barcodes
length(Purity.HNSC)
# 标准化处理
library(EDASeq)
## 标准化mRNA转录本和miRNA
dataNorm.OSCC <- TCGAanalyze_Normalization(tabDF = dataPrep,
                                           geneInfo = geneInfo,
                                           method = "gcContent")

# 将标准化后的数据再过滤，去除表达量较低的基因得到最终的数据

dataFilt.OSCC <- TCGAanalyze_Filtering(tabDF = dataNorm.OSCC,
                                       method = "quantile", 
                                       qnt.cut =  0.25)
dim(dataFilt.OSCC)
dataFilt.OSCC.final<-dataFilt.OSCC
dataFilt.OSCC.final=dataFilt.OSCC.final[rownames(dataFilt.OSCC.final)%in%rownames(dataPrep_final),]    # 数据过滤后剩下20个基因
dataFilt.OSCC.final=t(dataFilt.OSCC.final)
dataFilt.OSCC.final=dataFilt.OSCC.final[rownames(dataFilt.OSCC.final)%in%colnames(dataPrep_final),]    
dataFilt.OSCC.final=t(dataFilt.OSCC.final)
write.table(dataFilt.OSCC.final,
            file = "dataFilt.OSCC.final.xls",
            quote = F,
            row.names = T,
            sep = '\t')
# 02 表达差异分析
setwd("/data/nas1/luchunlin/project/BJTC-228")
if (! dir.exists("./01_DEG")){
  dir.create("./01_DEG")
}
setwd("./01_DEG")
## 02-1 OSCC和对照组的20个m6A甲基化调节因子的表达情况

## 数据清洗
dataPrep2<-dataFilt.OSCC.final
dataPrep2=log2(dataPrep2)
class(dataPrep2)
dataPrep2[dataPrep2==-Inf]<-0
dataPrep2<-as.data.frame(dataPrep2)
# dataPrep2=dataPrep2[rownames(dataPrep2)%in%rownames(dataFilt.OSCC.final),]
Symbol<-rownames(dataPrep2)
dataPrep2<-cbind(Symbol,dataPrep2)
dataPrep2<-gather(dataPrep2,
                  key = sample,
                  value = expr,
                  -c("Symbol"))

class(dataPrep2)
## 样本分组
group<-as.data.frame(group_list)
group<-group[group$sample%in%colnames(dataFilt.OSCC.final),]  #  取匹配的354个样本（310个肿瘤，44个对照）
write.table(group,file = "group.xls",
            quote = F,
            sep = "\t")
group_exp<-as.data.frame(group)
Group<-group_exp[,2]
dataPrep2<-cbind(Group,dataPrep2)

## 绘制箱线图

##分面图形
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
# 分面图形
exp_boxplot<-ggplot(dataPrep2,aes(x = Group, y = expr, fill = Group)) +
   geom_boxplot(alpha=0.7) +
   scale_y_continuous(name = "Expression",expand = c(0.1,0.1))+
   scale_x_discrete(name = "group") +
   ggtitle("Boxplot of 20 m6A regulator") +
   theme_bw() +
   geom_signif(comparisons = list(c("Normal","Tumor")),
              test = t.test,
              map_signif_level = T)+
   theme(plot.title = element_text(size = 16, face =  "bold"),
         text = element_text(size = 14),
         axis.title = element_text(face="bold"),
         axis.text.x=element_text(size = 16)) 
exp_boxplot
exp_boxplot2<-exp_boxplot+scale_fill_lancet()+facet_wrap(~Symbol,scales = "free")
exp_boxplot2
# 

exp_boxplot3<-ggplot(dataPrep2,aes(x = Symbol, y = expr, fill = Group)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "Expression")+
  scale_x_discrete(name = "Gene") +
  ggtitle("Boxplot of 20 m6A regulator") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1))+
  theme(plot.title = element_text(size = 14, face =  "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 11)) 
exp_boxplot3

## 02-2 差异分析
## DESeq2筛选差异基因
library(DESeq2)
expr<-dataFilt.OSCC.final
rownames(group)<-group$sample
group$group<-factor(group$group,levels = c("Normal","Tumor"))
dds<-DESeqDataSetFromMatrix(countData = expr,
                            colData = group,
                            design = ~group)
dds=dds[rownames(counts(dds))>1,]
## 对原始dds进行normalize
dds<-DESeq(dds)
## 显示dds信息
dds
# 提取DESeq2分析结果
## 使用Result函数提取差异分析结果。
## 将提取的差异分析结果定义为变量“res”。
## contrast：定义谁和谁比较
res =results(dds, contrast = c("group","Tumor","Normal"))
## 对结果res利用order（）函数按照pvalue值进行排序。
res =res[order(res$pvalue),]
head(res)
summary(res)
## 保存所有输出结果
write.csv(res,file="All_results.csv")
# 显示显著差异基因数目
table(res$padj<0.05)
# 对显著性差异结果进行提取和保存
## 获取padj< 0.05，表达倍数取以2为对数后大于1或者小于-1的差异表达基因。
## 使用subset（）函数过滤需要结果至变量DEG中
## Usage:subset(x,...),x为objects，...为筛选的参数或者条件
DEG <- subset(res, padj < 0.05 & abs(log2FoldChange) > 0)
DEG<-as.data.frame(res)
DEG<-na.omit(DEG)
## 使用dim查看该结果的维度、规模
dim(DEG)
head(DEG)
## 添加change列
logFC_cutoff<-0
DEG$change=as.factor(
  ifelse(DEG$padj<0.05&abs(DEG$log2FoldChange)>0,
         ifelse(DEG$log2FoldChange>logFC_cutoff,'UP','DOWN'),'NOT'))
table(DEG$change)
## down 6 up 8 

## 输出结果
write.csv(DEG,file = "DEG.csv")
sig_diff <- subset(DEG,
                   DEG$padj < 0.05 & abs(DEG$log2FoldChange) >= 0)

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
## 绘制热图
sig_diff_heatmap<-sig_diff[order(abs(sig_diff$log2FoldChange)[1:14],decreasing = T),]
sig_diff_heatmap<-sig_diff_heatmap[order(sig_diff_heatmap$change),]
mat<-assay(dds[rownames(sig_diff_heatmap),])
mat<-t(scale(t(mat)))

library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(gplots)
annotation_col <- data.frame(
  Group = factor(colData(dds)$group),
  row.names = rownames(colData(dds)))
ann_colors <- list(
  Group = c(Normal="lightblue", Tumor="darkorange"))

#绘制基因的热图，并调整参数

pheatmap(mat=mat,
         color=bluered(100),
         scale="row",
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         fontsize = 10,
         show_colnames = FALSE,
         cluster_cols = F)

# 03 GO_KEGG富集分析
setwd("/data/nas1/luchunlin/project/BJTC-228")
if (! dir.exists("./02_GO_KEGG")){
  dir.create("./02_GO_KEGG")
}
setwd("./02_GO_KEGG")
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
diff_gene_names <- rownames(sig_diff)
gene_transform <- bitr(diff_gene_names,
                       fromType = "SYMBOL",
                       toType = c("ENTREZID", "ENSEMBL", "REFSEQ"),
                       OrgDb = "org.Hs.eg.db")

ego <- enrichGO(gene = gene_transform$ENTREZID, 
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID", 
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.01, 
                qvalueCutoff = 0.05,
                readable = TRUE)
write.table(ego@result,file = "GO.xls",sep = "\t",quote = F,row.names = F)

## GO富集（条形图）
go_bar <- barplot(ego, showCategory=5, split="ONTOLOGY") + 
  facet_grid(ONTOLOGY ~ ., scales = "free")
go_bar
write_fig(go_bar,
          file = "GO_bar.pdf",
          width = 10,
          height = 9,
          devices = NULL,
          res = 600,
          show = F)
write_fig(go_bar,
          file = "GO_bar.png",
          width = 10,
          height = 9,
          devices = NULL,
          res = 600,
          show = F)

## KEGG富集（气泡图）
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "hsa",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.2,
                 qvalueCutoff = 2)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk@result,file = "KEGG.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk, showCategory=25)
kk_dot
kk_bar <- barplot(kk, showCategory=20,
                  title = "The KEGG of diff genes")
kk_bar
write_fig(kk_dot,
          file = "KEGG_dot.pdf",
          width = 10,
          height = 9,
          devices = NULL,
          res = 600,
          show = F)
write_fig(kk_dot,
          file = "KEGG_dot.png",
          width = 10,
          height = 9,
          devices = NULL,
          res = 600,
          show = F)

# 04 PPI网络构建
setwd("/data/nas1/luchunlin/project/BJTC-228")
if (! dir.exists("./03_PPI")){
  dir.create("./03_PPI")
}
setwd("./03_PPI")
## String数据库

## MCODE插件 or degree

# 05 诊断基因的筛选
setwd("/data/nas1/luchunlin/project/BJTC-228")
if (! dir.exists("./04_hub")){
  dir.create("./04_hub")
}
setwd("./04_hub")
## hub基因的ROC曲线， AUC>0.7
## 14个hub gene:IGF2BP1、IGF2BP2、IGF2BP3、RBM15B、HNRNPA2B1、YTHDC1、YTHDC2、HNRNPC、YTHDF2、RBM15、ALKBH5、METTL3、YTHDF1、ZC3H13
## 提取14个hub基因的表达矩阵，形成小矩阵
hub_gene<-c("IGF2BP1","IGF2BP2","IGF2BP3","RBM15B","HNRNPA2B1","YTHDC1","YTHDC2","HNRNPC","YTHDF2","RBM15","ALKBH5","METTL3","YTHDF1","ZC3H13")
hub_gene<-as.data.frame(hub_gene)
hub_gene_exp<-dataFilt.OSCC.final[rownames(dataFilt.OSCC.final)%in%hub_gene$hub_gene,]
## 绘制ROC曲线
library(pROC)
library(ggplot2)
library(geomROC)
hub_gene_exp<-t(hub_gene_exp)
hub_gene_exp<-cbind(group,hub_gene_exp)
hub_gene_exp<-hub_gene_exp[,-1]
type<-group[,2]
design<-model.matrix(~ -1+factor(type,levels = c("Normal","Tumor"),ordered = T))
colnames(design)<-c("Normal","Tumor")
rownames(design)=rownames(hub_gene_exp)
hub_gene_exp<-cbind(design,hub_gene_exp)

# IGF2BP1
roc_IGF2BP1<-roc(hub_gene_exp$group,hub_gene_exp$IGF2BP1,
                   levels=c("Normal","Tumor"))
plot(roc_IGF2BP1,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,      
     #auc.polygon=T,                      
     #auc.polygon.con="#fff7f7",           
     grid=c(0.5,0.2),                     
     grid.col=c("black","black"),          
     #print.thres=T,                         
     main="IGF2BP1 ROC curve",              
     col="#FF2E63",                       
     legacy.axes=T,)
# IGF2BP2
roc_IGF2BP2<-roc(hub_gene_exp$group,hub_gene_exp$IGF2BP2,
                 levels=c("Normal","Tumor"))
plot(roc_IGF2BP2,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,      
     #auc.polygon=T,                      
     #auc.polygon.con="#fff7f7",           
     grid=c(0.5,0.2),                     
     grid.col=c("black","black"),          
     #print.thres=T,                         
     main="IGF2BP2 ROC curve",              
     col="#FF2E63",                       
     legacy.axes=T,)
# IGF2BP3
roc_IGF2BP3<-roc(hub_gene_exp$group,hub_gene_exp$IGF2BP3,
                 levels=c("Normal","Tumor"))
plot(roc_IGF2BP3,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,      
     #auc.polygon=T,                      
     #auc.polygon.con="#fff7f7",           
     grid=c(0.5,0.2),                     
     grid.col=c("black","black"),          
     #print.thres=T,                         
     main="IGF2BP3 ROC curve",              
     col="#FF2E63",                       
     legacy.axes=T,)
# RBM15B
roc_RBM15B<-roc(hub_gene_exp$group,hub_gene_exp$RBM15B,
                 levels=c("Normal","Tumor"))
plot(roc_RBM15B,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,      
     #auc.polygon=T,                      
     #auc.polygon.con="#fff7f7",           
     grid=c(0.5,0.2),                     
     grid.col=c("black","black"),          
     #print.thres=T,                         
     main="RBM15B ROC curve",              
     col="#FF2E63",                       
     legacy.axes=T,)
# HNRNPA2B1
roc_HNRNPA2B1<-roc(hub_gene_exp$group,hub_gene_exp$HNRNPA2B1,
                levels=c("Normal","Tumor"))
plot(roc_HNRNPA2B1,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,      
     #auc.polygon=T,                      
     #auc.polygon.con="#fff7f7",           
     grid=c(0.5,0.2),                     
     grid.col=c("black","black"),          
     #print.thres=T,                         
     main="HNRNPA2B1 ROC curve",              
     col="#FF2E63",                       
     legacy.axes=T,)
# YTHDC1
roc_YTHDC1<-roc(hub_gene_exp$group,hub_gene_exp$YTHDC1,
                   levels=c("Normal","Tumor"))
plot(roc_YTHDC1,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,      
     #auc.polygon=T,                      
     #auc.polygon.con="#fff7f7",           
     grid=c(0.5,0.2),                     
     grid.col=c("black","black"),          
     #print.thres=T,                         
     main="YTHDC1 ROC curve",              
     col="#FF2E63",                       
     legacy.axes=T,)
# YTHDC2
roc_YTHDC2<-roc(hub_gene_exp$group,hub_gene_exp$YTHDC2,
                levels=c("Normal","Tumor"))
plot(roc_YTHDC2,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,      
     #auc.polygon=T,                      
     #auc.polygon.con="#fff7f7",           
     grid=c(0.5,0.2),                     
     grid.col=c("black","black"),          
     #print.thres=T,                         
     main="YTHDC2 ROC curve",              
     col="#FF2E63",                       
     legacy.axes=T,)
# HNRNPC
roc_HNRNPC<-roc(hub_gene_exp$group,hub_gene_exp$HNRNPC,
                levels=c("Normal","Tumor"))
plot(roc_HNRNPC,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,      
     #auc.polygon=T,                      
     #auc.polygon.con="#fff7f7",           
     grid=c(0.5,0.2),                     
     grid.col=c("black","black"),          
     #print.thres=T,                         
     main="HNRNPC ROC curve",              
     col="#FF2E63",                       
     legacy.axes=T,)
# YTHDF2
roc_YTHDF2<-roc(hub_gene_exp$group,hub_gene_exp$YTHDF2,
                levels=c("Normal","Tumor"))
plot(roc_YTHDF2,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,      
     #auc.polygon=T,                      
     #auc.polygon.con="#fff7f7",           
     grid=c(0.5,0.2),                     
     grid.col=c("black","black"),          
     #print.thres=T,                         
     main="YTHDF2 ROC curve",              
     col="#FF2E63",                       
     legacy.axes=T,)
# RBM15
roc_RBM15<-roc(hub_gene_exp$group,hub_gene_exp$RBM15,
                levels=c("Normal","Tumor"))
plot(roc_RBM15,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,      
     #auc.polygon=T,                      
     #auc.polygon.con="#fff7f7",           
     grid=c(0.5,0.2),                     
     grid.col=c("black","black"),          
     #print.thres=T,                         
     main="RBM15 ROC curve",              
     col="#FF2E63",                       
     legacy.axes=T,)
# ALKBH5
roc_ALKBH5<-roc(hub_gene_exp$group,hub_gene_exp$ALKBH5,
               levels=c("Normal","Tumor"))
plot(roc_ALKBH5,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,      
     #auc.polygon=T,                      
     #auc.polygon.con="#fff7f7",           
     grid=c(0.5,0.2),                     
     grid.col=c("black","black"),          
     #print.thres=T,                         
     main="ALKBH5 ROC curve",              
     col="#FF2E63",                       
     legacy.axes=T,)
# METTL3
roc_METTL3<-roc(hub_gene_exp$group,hub_gene_exp$METTL3,
                levels=c("Normal","Tumor"))
plot(roc_METTL3,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,      
     #auc.polygon=T,                      
     #auc.polygon.con="#fff7f7",           
     grid=c(0.5,0.2),                     
     grid.col=c("black","black"),          
     #print.thres=T,                         
     main="METTL3 ROC curve",              
     col="#FF2E63",                       
     legacy.axes=T,)
# YTHDF1
roc_YTHDF1<-roc(hub_gene_exp$group,hub_gene_exp$YTHDF1,
                levels=c("Normal","Tumor"))
plot(roc_YTHDF1,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,      
     #auc.polygon=T,                      
     #auc.polygon.con="#fff7f7",           
     grid=c(0.5,0.2),                     
     grid.col=c("black","black"),          
     #print.thres=T,                         
     main="YTHDF1 ROC curve",              
     col="#FF2E63",                       
     legacy.axes=T,)
# ZC3H13
roc_ZC3H13<-roc(hub_gene_exp$group,hub_gene_exp$ZC3H13,
                levels=c("Normal","Tumor"))
plot(roc_ZC3H13,
     print.auc=T,
     print.auc.x=0.4,print.auc.y=0.5,      
     #auc.polygon=T,                      
     #auc.polygon.con="#fff7f7",           
     grid=c(0.5,0.2),                     
     grid.col=c("black","black"),          
     #print.thres=T,                         
     main="ZC3H13 ROC curve",              
     col="#FF2E63",                       
     legacy.axes=T,)

# 06 诊断基因的表达相关性分析
setwd("/data/nas1/luchunlin/project/BJTC-228")
if (! dir.exists("./05_correlation")){
  dir.create("./05_correlation")
}
setwd("./05_correlation")
# 筛选到的9个具有诊断意义的hub基因：IGF2BP1,IGF2BP2,IGF2BP3,HNRNPA2B1,HNRNPC,RBM15,YTHDC2,YTHDF1,RBM15B
library(ggcorrplot)
library(corrplot)
hub_gene_final<-c("IGF2BP1","IGF2BP2","IGF2BP3","HNRNPA2B1","HNRNPC","RBM15","YTHDF1","YTHDC2","RBM15B")
hub_gene_final<-as.data.frame(hub_gene_final)
#提取9个hub基因的表达矩阵
hub_gene_final_exp<-dataFilt.OSCC.final[rownames(dataFilt.OSCC.final)%in%hub_gene_final$hub_gene_final,]
hub_corr<-round(cor(t(hub_gene_final_exp)),3)

## 检验基因之间的相关性p值
## 计算相关性系数并显示基因之间的相关性。相关性系数大于0为正相关，小于0为负相关。
## p小于0.05认为相关性显著
hub_p.mat<-round(cor_pmat(t(hub_gene_final_exp)),3)
## 对基因聚类，聚类算法为ward。D
col1 <- colorRampPalette(c("#67001F", "#B2182B", "#D6604D", "#F4A582",
                           "#FDDBC7", "#FFFFFF", "#D1E5F0", "#92C5DE",
                           "#4393C3", "#2166AC", "#053061"))
hub_corr_plot<-corrplot(hub_corr,
                        method = "square",
                        is.corr = T,
                        type = "lower",
                        p.mat = hub_p.mat,
                        insig = "blank",
                        outline = "white",
                        addCoef.col ="black",
                        col = col1(200))


hub_corr_plot<-ggcorrplot(hub_corr, 
                          method = "square",
                          hc.order = TRUE,
                          hc.method = "ward.D",
                          outline.color = "white",
                          ggtheme = theme_bw(),
                          type = "upper",
                          colors = c("#FF69B4", "white", "#6495ED"),
                          lab = TRUE,
                          lab_size = 5,
                          p.mat = hub_p.mat,
                          insig = "blank")
hub_corr_plot
write_fig(hub_corr_plot,
          file = "correlation_hub_p.pdf",
          width = 10,
          height = 10,
          devices = NULL,
          res = 600,
          show = F)
write_fig(hub_corr_plot,
          file = "correlation_hub_p.png",
          width = 10,
          height = 10,
          devices = NULL,
          res = 600,
          show = F)
# 07 诊断基因调控网络预测
# 选择在验证集验证成功的诊断基因，使用StarBase和PROMO数据可对诊断基因的miRNAs和TFs进行预测，将
# 诊断基因及其miRNA和TF数据整合到一个调控网络中，使用Cytoscape可视化。

save.image(file = OSCC.RData)

