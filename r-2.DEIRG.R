## 02 差异基因鉴定-----
rm(list = ls())
setwd("/data/nas1/luchunlin/project/YCZK-127/")
if (! dir.exists("./02_DEIRG")){
  dir.create("./02_DEIRG")
}
setwd("./02_DEIRG")
sig.diff<-read.delim2('/data/nas1/luchunlin/project/YCZK-127/01_DEGs/DEG_sig.xls',row.names = 1)
invasion<-read.csv(file = 'cell search result.csv')
#invasion<-read_xlsx('invasion.xlsx')
DEIRG<-invasion[invasion$Symbol%in%rownames(sig.diff),]
DEIRG<-DEIRG[,c(1,2)]
##75
library(ggvenn)
mydata<-list('DEGs'=rownames(sig.diff),'Invasion'=invasion$Symbol)
ggvenn(mydata,c('DEGs','Invasion'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
write.table(DEIRG,file = 'DEIRG.xls',sep = '\t',quote = F,row.names = F)

## GO/KEGG
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(GOplot)
gene_transform <- bitr(DEIRG$Symbol,
                       fromType = "SYMBOL",
                       toType = c("ENTREZID"),
                       OrgDb = "org.Hs.eg.db")
ego <- enrichGO(gene = gene_transform$ENTREZID,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
CC <- enrichGO(gene = gene_transform$ENTREZID,
               OrgDb = org.Hs.eg.db,
               keyType = "ENTREZID",
               ont = "CC",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               readable = TRUE)
BP <- enrichGO(gene = gene_transform$ENTREZID,
               OrgDb = org.Hs.eg.db,
               keyType = "ENTREZID",
               ont = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               readable = TRUE)
MF <- enrichGO(gene = gene_transform$ENTREZID,
               OrgDb = org.Hs.eg.db,
               keyType = "ENTREZID",
               ont = "MF",
               pAdjustMethod = "BH",
               pvalueCutoff = 0.05,
               qvalueCutoff = 0.05,
               readable = TRUE)
write.table(ego,file = "GO.xls",sep = "\t",quote = F,row.names = F)

go2<-data.frame(Category = "ALL",ID = ego$ID,Term = ego$Description, Genes = gsub("/", ", ", ego$geneID), adj_pval = ego$p.adjust)
## 读取logFC文件
sigind<-sig.diff[rownames(sig.diff)%in%DEIRG$Symbol,]
genelist<-data.frame(ID=rownames(sigind),
                     logFC=sigind$log2FoldChange)
genelist$logFC<-as.numeric(genelist$logFC)
rownames(genelist)<-genelist$ID
genelist<-genelist[order(genelist$logFC,decreasing = T),]
library(enrichplot)
colnames(gene_transform)<-c('ID','ENTREZID')
geneList<-merge(genelist,gene_transform,by='ID')
geneList<-geneList[,c(2:3)]
m<-geneList$logFC
m<-as.numeric(m)
names(m)<-geneList$ENTREZID
m
p1<-heatplot(CC,foldChange = m,showCategory = 5)
p1
p2<-heatplot(BP,foldChange = m,showCategory = 5)
p2
p3<-heatplot(MF,foldChange = m,showCategory = 5)
p3
go_bp<-mutate(BP,qscore= -log(p.adjust,base = 10))%>%
  barplot(x='qscore',showCategory=5)
go_bp
go_cc<-mutate(CC,qscore= -log(p.adjust,base = 10))%>%
  barplot(x='qscore',showCategory=5)
go_cc

go_mf<-mutate(MF,qscore= -log(p.adjust,base = 10))%>%
  barplot(x='qscore',showCategory=5)
go_mf
# 展示富集最显著的 GO term
go_dot<-dotplot(ego, showCategory=5, split="ONTOLOGY") +
  facet_grid(ONTOLOGY ~ ., scales = "free")
go_dot

## KEGG富集分析（气泡图）
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "hsa",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk@result,file = "KEGG.xls",sep = "\t",quote = F,row.names = F)
kk2<-data.frame(Category = "ALL",ID = kk$ID,Term = kk$Description, Genes = gsub("/", ", ", kk$geneID), adj_pval = kk$p.adjust)
## 读取logFC文件
circ<-circle_dat(kk2,genelist)
# GOBubble(circ, labels = 3,table.legend =F)
# GOCircle(circ,rad1=2.5,rad2=3.5,label.size=4,nsub=10) 
#termNum<-20
#geneNum<-nrow(genelist)
chord<-chord_dat(circ,genelist)
kegg_chord<-GOChord(chord,
                    gene.order = 'logFC',
                    gene.space = 0.25,
                    gene.size = 9,
                    space = 0.01,
                    lfc.col = c('red','white','blue'),
                    process.label = 10)
kegg_chord

p4<-heatplot(kk,foldChange = m,showCategory = 10)
p4
keggbar<-mutate(kk,qscore= -log(p.adjust,base = 10))%>%
  barplot(x='qscore',showCategory=10)
keggbar
