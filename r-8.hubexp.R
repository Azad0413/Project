## 09 hub.exp-----
rm(list = ls())
setwd("/data/nas1/luchunlin/project/LLZK-505")
if (! dir.exists("./08_hubexp")){
  dir.create("./08_hubexp")
}
setwd("./08_hubexp")
## 表达分析  KM  GSEA
dat<-read.delim2('/data/nas1/luchunlin/project/LLZK-505/00_rawdata/dat.xls',row.names = 1)%>%lc.tableToNum()
group<-read.delim2('/data/nas1/luchunlin/project/LLZK-505/00_rawdata/group.xls')
control.sample<-group$sample[which(group$group=='control')]
ov.sample<-group$sample[which(group$group=='OV')]
## 09-1 GRIM-19  NDUFS3------
### 09-1-1 表达分析------
hubgene<-c('EYA1','SOX21','DLX2','POU3F3','EMX1')
hub_exp<-dat[hubgene,]
hub_exp2<-hub_exp
hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))

## 样本分组
hub_exp2$Group<-ifelse(hub_exp2$sample%in%control.sample,'control','OV')

##分面图形
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
exp_plot <- ggplot(hub_exp2,aes(x = Symbol, y = expr, fill = Group)) +
  #geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#4682B4","#CD3700"), name = "Group")+
  labs(title="", x="", y = "",size=20) +
  stat_compare_means(data = hub_exp2,
                     mapping = aes(group = Group),
                     label ="p.signif",
                     method = 't.test') +
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=15),
        axis.text.x=element_text(angle=0,hjust=0.5,colour="black",face="bold",size=12), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+facet_wrap(~Symbol,scales = "free",nrow = 2) 
exp_plot



## 09-2 NDUFA4,LRPPRC-----
### 09-2-1 表达分析------
hubgene<-c('EYA1','SOX21','DLX2','GFAP','SOX3','LIN28A')
hub_exp<-dat[hubgene,]
hub_exp2<-hub_exp
hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))

## 样本分组
hub_exp2$Group<-ifelse(hub_exp2$sample%in%control.sample,'control','OV')

##分面图形
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
exp_plot <- ggplot(hub_exp2,aes(x = Symbol, y = expr, fill = Group)) +
  #geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#4682B4","#CD3700"), name = "Group")+
  labs(title="", x="", y = "",size=20) +
  stat_compare_means(data = hub_exp2,
                     mapping = aes(group = Group),
                     label ="p.signif",
                     method = 't.test') +
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=15),
        axis.text.x=element_text(angle=0,hjust=0.5,colour="black",face="bold",size=12), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+facet_wrap(~Symbol,scales = "free",nrow = 2) 
exp_plot


### 09-2-2 KM-----
cluster2$cluster<-as.vector(cluster2$cluster)
km.dat<-t(dat)%>%as.data.frame()
km.dat$group<-as.vector(cluster2$cluster,)
km.dat$sample<-rownames(km.dat)
km.dat<-merge(survival,km.dat,by='sample')
rownames(km.dat)<-km.dat$sample
km.dat<-km.dat[,-c(1,3)]
table(km.dat$group)
kmfit<-survfit(Surv(OS.time, OS) ~ group, data =  km.dat)
cluster_survival_median <- ggsurvplot(kmfit,
                                      pval = TRUE, 
                                      conf.int = F,
                                      legend.labs=c("Cluster1","Cluster2" ),
                                      legend.title="group",
                                      title="Train KM",
                                      font.main = c(15,"bold"),
                                      risk.table = TRUE, 
                                      risk.table.col = "strata", 
                                      linetype = "strata", 
                                      surv.median.line = "hv", 
                                      ggtheme = theme_bw(), 
                                      palette = c("#A73030FF", "#0073C2FF"))
cluster_survival_median
### 09-2-3 GSEA-----
colData<-cluster2
colnames(colData)<-c('group','sample')
table(colData$group)
colData$group<-factor(colData$group,levels = c('Cluster1','Cluster2'))
dds<-DESeqDataSetFromMatrix(countData = dat,colData=colData,design = ~group)
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
logFCcutoff <- 1
allGeneSets$change = as.factor(
  ifelse(allGeneSets$padj < 0.05 & abs(allGeneSets$log2FoldChange) > logFCcutoff,
         ifelse(allGeneSets$log2FoldChange > logFCcutoff,'UP','DOWN'),'NOT')
)

DEGeneSets <- subset(allGeneSets,
                     allGeneSets$padj < 0.05 & abs(allGeneSets$log2FoldChange) > logFCcutoff)
genelist <- allGeneSets$log2FoldChange
names(genelist) <- rownames(allGeneSets)
geneList <- sort(genelist, decreasing = T)
DEGeneSets <- DEGeneSets[order(DEGeneSets$padj),]
dim(DEGeneSets)
## GSEA KEGG
library(clusterProfiler)
library(enrichplot)
kegg_set<- read.gmt("c2.cp.kegg.v7.5.1.symbols.gmt")
set.seed(1)
kegg_gsea <- GSEA(geneList, TERM2GENE = kegg_set, pvalueCutoff = 0.05)
kegg_result <- kegg_gsea@result
dim(kegg_result)
write.table(kegg_result,file = 'GSEA(GRIM19&NDUFS3).xls',sep = '\t',quote = F,row.names = F)
gseaplot2(kegg_gsea,c(1:5),color = c('#7B68EE','#CD3333','#20B2AA','#FF8C00','#FF6666'),
          title = 'KEGG GSEA(GRIM19&NDUFS3)',
          base_size = 11,
          rel_heights = c(1.5, 0.3, 0.5))


