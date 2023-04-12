rm(list = ls())
# 设置当前工作环境 ---------------------------------------------------------------
setwd("/data/nas1/luchunlin/project/BJTC-334")
if (! dir.exists("./00_rawdata")){
  dir.create("./00_rawdata")
}
setwd("./00_rawdata")
library(GEOquery)    ##加载包
gset <- getGEO('GSE23561',   ##GSE73094
               AnnotGPL = F,
               getGPL = F)
#save(gset,file = '00_Rawdata/GSE23561.gset.Rdata')   
eset<-exprs(gset[[1]])   ##表达矩阵  
metadata<-pData(gset[[1]]) ##临床信息

# ID转换 --------------------------------------------------------------------
gpl <- getGEO('GPL10775')    ##数据集对应的平台文件，若无法下载可直接去GEO官网下载
anno1<-Table(gpl)[,c(1,5)]    #提取id和gene symbol 
colnames(anno1)<-c('ID','Gene Symbol')  
head(anno1)
eset<-as.data.frame(eset)  
eset$ID <- rownames(eset)   
head(eset)
merg<-merge(eset,anno1,by="ID")  ##平台ID和表达量矩阵ID合并
y<-merg$`Gene Symbol`
gene<-unlist(lapply(y,function(y) strsplit(as.character(y)," /// ")[[1]][1]))
merg$gene <- gene  
aggr<-aggregate(merg[,2:36],by=list(merg$gene),mean)   ##371为样本数+1
aggr<-aggr[!duplicated(aggr$Group.1),]
rownames(aggr)<-aggr[,1]    ###以基因名为列名 
aggr<-aggr[,-1]       
expr<-aggr
condition<-data.frame(sample=rownames(metadata),group=metadata$characteristics_ch1)
condition<-subset(condition,group=='disease state: Control'|group=='disease state: Coronary Artery Disease')
condition$group<-ifelse(condition$group=='disease state: Control','Control','CAD')
Group<-data.frame(row.names=condition$sample,group=condition$group)
expr<-subset(expr,select=rownames(Group))
#write.table(expr,'00_Rawdata/expr_traing.txt',sep='\t')
#write.table(Group,'00_Rawdata/Group_traing.txt',sep='\t')
mRNA<-read.table('mRNA.txt')
c<-intersect(rownames(expr),mRNA$gene_name)
exprSet = expr[c,]   ###提取mrna  ##13717个基因  15个样本（正常9，患病6）

# 差异分析 --------------------------------------------------------------------
library(limma)
group<- factor(Group$group,levels = c("CAD","Control"))
condition<-data.frame(group)
rownames(condition)<-colnames(exprSet )
design<-model.matrix(~-1+group)   
contrast.matrix<-makeContrasts(contrasts = "groupCAD-groupControl", levels = design) 
fit <- lmFit(log2(exprSet+1),design)     
fit1 <- contrasts.fit(fit, contrast.matrix)    
fit2 <- eBayes(fit1)  
tempOutput<- topTable(fit2, coef=1, n=nrow(fit2),lfc=log2(1),adjust="fdr")  ###所有基因检验结果
# write.table(tempOutput,'01_Degs/tempOutput.xls',sep='\t',quote=FALSE)   ##导出所有基因检验结果
dif<-tempOutput[which(tempOutput$adj.P.Val<0.05&abs(tempOutput$logFC)>0),]    ##差异基因检验结果
dim(dif)   ##9485个差异基因 
dif_up<-tempOutput[which(tempOutput$adj.P.Val<0.05&tempOutput$logFC>0),]    ##差异基因检验结果
dim(dif_up)  ##9159个上调，326个下调
###mRNA
#write.table(dif,'01_Degs/dif.xls',sep = '\t', quote = FALSE)

# 火山图 ---------------------------------------------------------------------
library(ggplot2)
library(ggrepel)
Res<-na.omit(tempOutput)
Res$sig[Res$adj.P.Val>=0.05 | abs(Res$logFC) <0] <- "NO"    ###GREEN
Res$sig[Res$adj.P.Val<0.05 & Res$logFC >0] <- "Up"     ###BLACK
Res$sig[Res$adj.P.Val<0.05 & Res$logFC <= 0] <- "Down" ###RED  
p<-ggplot(
  #设置数据
  Res, 
  aes(x = logFC, 
      y = -log10(adj.P.Val), 
      color=sig)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("blue", "grey","red"))+
  # 辅助线
  # geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  
  # 坐标轴
  labs(x="logFC",
       y="-log10(P-Value)")+
  theme_bw()+
  # 图例
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank()
  )
p
# df<-df[order(df$logFC, decreasing = T),]  ##差异基因根据logFC由大到小进行排序
# gene<-df[c(1:10,9476:9485),]  ##标出前10个基因和后10个基因
# gene$label<-rownames(gene)
# p+geom_text_repel(data = gene, aes(x=logFC,y=-log10(adj.P.Val), 
#                                    label = label),
#                   size = 3,box.padding = unit(0.5, "lines"),
#                   point.padding = unit(0.8, "lines"), 
#                   segment.color = "black", 
#                   show.legend = FALSE)
ggsave('01_DEGs/volcano.png',width=10,height=8)
ggsave('01_DEGs/volcano.pdf',width=10,height=8)

# 热图 ----------------------------------------------------------------------
library(edgeR)
library(pheatmap)
a<-condition
colnames(a)<-'group'  ##a为样本分组情况
df<-df[order(df$logFC , decreasing = T),]
gene<-df[c(1:10,9476:9485),]   ##前10个上调基因和后10个上调基因
Direction<-c(rep('Up',10),rep('Down',10))
type<-data.frame(row.names=rownames(gene),Direction)
select_df<-exprSet [rownames(gene),]
color.key <- c("#3300CC", "#3399FF", "white", "#FF3333", "#CC0000")
pheatmap(select_df,scale='row',annotation_col=a ,annotation_row=type,annotation_names_col = F,annotation_names_row = F,cluster_rows = FALSE,cluster_cols = F,show_colnames = FALSE,width=9,height=7,file='01_Degs/heat.png')
pheatmap(select_df ,scale='row',annotation_col=a,annotation_row=type,annotation_names_col = F,annotation_names_row = F,cluster_rows = FALSE,cluster_cols = F,show_colnames = FALSE,width=9,height=7,file='01_Degs/heat.pdf')

# GO富集分析 ------------------------------------------------------------------
GO_database <- 'org.Hs.eg.db'
KEGG_database <- 'hsa' 
library(clusterProfiler)
library(org.Hs.eg.db)
###GO富集分析
gene <- bitr(rownames(df),fromType = 'SYMBOL',toType = 'ENTREZID',OrgDb = GO_database)
GO<-enrichGO( gene$ENTREZID,
              OrgDb = GO_database,
              keyType = "ENTREZID",             
              ont = "ALL",  ##注释
              pAdjustMethod = "BH",##p值调整方法
              readable = T)
library(ggplot2)
barplot(GO, showCategory = 10,split="ONTOLOGY")+labs(title = "Gene Ontology ")+facet_grid(ONTOLOGY~., scale="free")  ##GO富集柱状图
ggsave(width=12,height=10,'01_Degs/GO.png')
ggsave(width=12,height=10,'01_Degs/GO.pdf')
GO<-data.frame(GO)  ##富集结果
write.table(GO,file="01_Degs/GO.xls",sep="\t",row.names=F,quote=F)

# KEGG富集分析 ----------------------------------------------------------------
KEGG<-enrichKEGG(gene$ENTREZID,
                 organism = KEGG_database)
dotplot(KEGG,showCategory = 10,title = 'KEGG Pathway')  ##KEGG气泡图
ggsave(width=8,height=6,'01_Degs/KEGG.pdf')
ggsave(width=8,height=6,'01_Degs/KEGG.png')
KEGG<-data.frame(KEGG)  ##KEGG富集结果
write.table(KEGG,file="01_Degs/KEGG.xls",sep="\t",row.names=F,quote=F)
