rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/XA-0214-1/")
if (! dir.exists("./07_GSEA")){
  dir.create("./07_GSEA")
}
setwd("./07_GSEA")
library(lance)
library(tidyverse)
hubgene <- read.delim2('../04_model/hubgene.xls')
dat <- read.delim2('../00_rawdata/dat(GSE97537).xls',row.names = 1)%>%lc.tableToNum()

hub.exp <- dat[hubgene$symbol,]
library(data.table)
library(clusterProfiler)
library(biomaRt)
library(enrichplot)
library(limma)
## 
# kegggmt<-read.gmt("/data/nas1/luchunlin/pipeline/GSVA/c2.cp.kegg.v7.5.1.symbols.gmt") #读gmt文件

library(org.Rn.eg.db)
library(tidyverse)
library(msigdbr)
msigdbr_species()
rat<-msigdbr(species = "Rattus norvegicus")
rat[1:5,1:5]
table(rat$gs_cat) ## 查看目录，与MSigDB一样包含9个数据集
# 查看需要的数据集在哪里
table(rat$gs_subcat)
rat_KEGG<-msigdbr(species = "Rattus norvegicus",
                   category = 'C2',
                   subcategory = "CP:KEGG")%>%
 dplyr::select(gs_name,gene_symbol)
head(rat_KEGG)
class(rat_KEGG)
# rat_KEGGSet<-rat_KEGG%>%split(x = .$gene_symbol, f = .$gs_name)
colnames(rat_KEGG) <- c('term','gene')
class(kegggmt$term)
rat_KEGG <- as.data.frame(rat_KEGG)
rat_KEGG$term <- as.factor(rat_KEGG$term)
kegggmt <- rat_KEGG

hub.exp<-hub.exp[order(rownames(hub.exp),decreasing = F),]  ##按照字母进行排序

# 每个hub基因均分析一遍 ------------------------------------------------------------
##分组
for (i in c(1:5)){
  train<-t(hub.exp)
  ## 以hub基因将样本进行高低表达量分析
  group=as.vector(ifelse(train[,i]>median(train[,i]),'high','low'))  ##根据TNFRSF1A表达量的大小将样本分组
  group <- factor(group,levels = c("high","low")) 
  condition<-data.frame(group)
  rownames(condition)<-rownames(train)
  design<-model.matrix(~0+group)   
  contrast.matrix<-makeContrasts(contrasts = "grouphigh-grouplow", levels = design) 
  fit <- lmFit(dat,design)     
  fit1 <- contrasts.fit(fit, contrast.matrix)    
  fit2 <- eBayes(fit1)  
  tempOutput <- topTable(fit2, coef=1, n=nrow(fit2),adjust="fdr")  ###所有基因检验结果
  genelist<-data.frame(SYMBOL=rownames(tempOutput ),logFC=tempOutput $logFC)
  #开始ID转换
  gene=bitr(genelist$SYMBOL,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Rn.eg.db") #会有部分基因数据丢失，或者ENSEMBL
  ## 去重
  gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
  gene_df <- data.frame(ENTREZID=gene$ENTREZID ,#可以是foldchange
                        SYMBOL = gene$SYMBOL) #记住你的基因表头名字
  gene_df <- merge(gene_df,genelist,by="SYMBOL")
  geneList<-gene_df $logFC #第二列可以是folodchange，也可以是logFC
  names(geneList)=gene_df $SYMBOL #使用转换好的ID
  geneList=sort(geneList,decreasing = T) #从高到低排序
  
  KEGG<-GSEA(geneList,TERM2GENE = kegggmt,pvalueCutoff = 0.05) #GSEA分析
  sortKEGG<-data.frame(KEGG)
  sortKEGG<-sortKEGG[order(sortKEGG$p.adjust, decreasing = F),]#按照enrichment score从高到低排序
  write.table(sortKEGG,file=paste('0',i,'.',rownames(hub.exp)[i],'(KEGG)','.xls',sep=''),sep="\t",quote=F,row.names=F)
  paths <- rownames(sortKEGG[c(1:5),])#选取你需要展示的通路ID
  
  #特定通路作图
  #trace(gseaplot2,edit=T)
  p = gseaplot2(KEGG,paths,base_size =10,color = c('#7B68EE','#CD3333','#20B2AA','#FF8C00','#FF6666'),rel_heights = c(1.5, 0.3, 0.5),title=rownames(hub.exp)[i])
  fn1 = paste0("0", i, ".", rownames(hub.exp)[i],'(KEGG)', ".png")
  fn2 = paste0("0", i, ".", rownames(hub.exp)[i],'(KEGG)',".pdf")
  ggsave(fn1, p, width = 8, height = 5, units = "in", limitsize = 300)
  ggsave(fn2, p, width = 8, height = 5, units = "in", limitsize = 300)
  i<-i+1
}

