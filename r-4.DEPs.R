rm(list = ls())
# 01 获取数据集--------------
setwd("/data/nas1/luchunlin/project/BJTC-356/")
if (! dir.exists("./04_DEPs")){
  dir.create("./04_DEPs")
}
setwd("./04_DEPs")
library(tidyverse)
library(lance)
library(readxl)
result <- read_xlsx('result_1_2.xlsx')
colnames(result)
result$logFC <- log2(result$`1_2`)
sig_DEP <- subset(result,result$p_value<0.05)%>%
  column_to_rownames(var = 'Accession')

# sig_DEP <- subset(result,result$p_value<0.05 & result$`1_2`<0.909 |result$p_value<0.05 &result$`1_2`>1.1)%>%
#   column_to_rownames(var = 'Accession')
colnames(sig_DEP)

### logFC t p.value adj.Pvalue B
sig_DEP <- sig_DEP[,c(13,11,12)]
colnames(sig_DEP)
colnames(sig_DEP) <- c('logFC','p.value','1_2')
sig_DEP$change <- ifelse(sig_DEP$logFC>0,'UP','DOWN')
table(sig_DEP$change)
result$change <- ifelse(result$logFC>0 & result$p_value<0.05,'UP',
                        ifelse(result$logFC<0 & result$p_value <0.05 ,'DOWN','NOT'))
write.table(sig_DEP,file = 'sig_DEP.xls',sep = '\t',row.names = T,quote = F)
# 
colnames(result)
colnames(sig_DEP)
volcano_plot<- ggplot(data = result, 
                      aes(x = logFC,
                          y = -log10(p_value), 
                          color =change)) +
  scale_color_manual(values = c("#20B2AA", "darkgray","red")) +
  scale_x_continuous(breaks = c(-10,-6,-3,-1,0,1,3,6,10)) +
  scale_y_continuous(trans = "log1p",
                     breaks = c(0,1,5,10,20,50, 100,200)) +
  geom_point(size = 1.2, alpha = 0.4, na.rm=T) +
  theme_bw(base_size = 12, base_family = "Times") +
  # geom_vline(xintercept = c(-0.5,0.5),
  #            lty = 4,
  #            col = "darkgray",
  #            lwd = 0.6)+
  geom_hline(yintercept = -log10(0.05),
             lty = 4,
             col = "darkgray",
             lwd = 0.6)+
  theme(legend.position = "right",
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(face="bold",
                                   color="black",
                                   family = "Times",
                                   size=13),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face = "bold",
                                   color = "black",
                                   size = 15),
        axis.text.y = element_text(face = "bold",
                                   color = "black",
                                   size = 15),
        axis.title.x = element_text(face = "bold",
                                    color = "black",
                                    size = 15),
        axis.title.y = element_text(face = "bold",
                                    color = "black",
                                    size = 15)) +
  # geom_label_repel(
  #   data = dat_rep,
  #   aes(label = rownames(dat_rep)),
  #   max.overlaps = 20,
  #   size = 4,
  #   box.padding = unit(0.5, "lines"),
  #   min.segment.length = 0,
  #   point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )+
  labs(x = "log (Fold Change)",
       y = "-log10 (P.Val)")
volcano_plot
ggsave('01.volcano.png', volcano_plot,width = 8, height = 6)
ggsave('01.volcano.pdf', volcano_plot,width = 8, height = 6)

# ### 热图-----
# dat <- read.delim2('../00_rawdata/protein.dat.xls',row.names = 1)%>%lc.tableToNum()
# #dat <- log2(dat+1)
# table(sig_DEP$change)
# 
# dat_rep<-rbind(head(sig_DEP2[order(sig_DEP2$logFC,decreasing = T),],20),
#                head(sig_DEP2[order(sig_DEP2$logFC,decreasing = F),],20))
# 
# df.group <- data.frame(sample=colnames(dat),group=c(rep('Persistent',3),rep('Paroxysmal',3)))
# group_rt<-df.group%>%as.data.frame()
# group_rt<-group_rt[order(group_rt$group),]
# rt<-dat[,group_rt$sample]
# group_rt<-data.frame(group_rt$group)
# colnames(group_rt)<-'group'
# rownames(group_rt)<-colnames(rt)
# heat<-rt[rownames(dat_rep),]%>%lc.tableToNum()
# 
# x<-heat
# #x<-t(scale(t(heat)))
# ann_colors<-list(
#   group = c(Paroxysmal="#00CED1",Persistent="#F08080"))
# pdf('02.pos.heatmap.pdf',w=6,h=6)
# pheatmap(mat=x,
#          annotation_col = group_rt,
#          color=bluered(100),
#          scale = "row",
#          annotation_colors = ann_colors,
#          fontsize = 10,
#          show_colnames = F,
#          cluster_cols = F,
#          cluster_rows = T,
#          show_rownames = F,
#          annotation_names_row = F)
# dev.off()
# png('02.pos.heatmap.png',w=500,h=500)
# pheatmap(mat=x,
#          annotation_col = group_rt,
#          color=bluered(100),
#          scale = "row",
#          annotation_colors = ann_colors,
#          fontsize = 10,
#          show_colnames = F,
#          cluster_cols = F,
#          cluster_rows = T,
#          show_rownames = F,
#          annotation_names_row = F)
# dev.off()
# 


# ##转化ID
# anno <- read_xlsx('/data/nas1/luchunlin/project/BJTC-321(drug)/01_Target/uniprot-filtered-organism__Homo+sapiens+(Human)+[9606]_.xlsx')
# protein2symbol <- anno[,c(1,4,5)]
# protein2symbol <- protein2symbol[protein2symbol$Entry%in%sig_DEP$Accession,]
# 
# 
library(clusterProfiler)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(GOplot)
gene_transform <- bitr(rownames(sig_DEP),
                       fromType = "UNIPROT",
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
write.table(ego,file = "GO.xls",sep = "\t",quote = F,row.names = F)
# 展示富集最显著的 GO term
go_bar<-barplot(ego, showCategory=5, split="ONTOLOGY",label_format=50) +
  facet_grid(ONTOLOGY ~ ., scales = "free")
go_bar
# ggsave(filename = '01.GO_bar.pdf',go_bar,w=7.5,h=6)
# ggsave(filename = '01.GO_bar.png',go_bar,w=7.5,h=6)

## GO 圈图
go_result<-read.table('GO.xls',header = T,sep = '\t',check.names = F)
go2=data.frame(Category='ALL',ID=go_result$ID,Term=go_result$Description,Genes=gsub("/", ", ", go_result$geneID), adj_pval = go_result$p.adjust)
df<-sig_DEP%>%rownames_to_column(var = 'UNIPROT')
gene_transform2 <- bitr(rownames(sig_DEP),
                       fromType = "UNIPROT",
                       toType = c("SYMBOL"),
                       OrgDb = "org.Hs.eg.db")

df <- merge(df,gene_transform2,by='UNIPROT')

idfc2<-na.omit(df)
idfc2 <- idfc2[!duplicated(idfc2$SYMBOL),]
rownames(idfc2) <- idfc2$SYMBOL

genelist<-data.frame(ID=rownames(idfc2),logFC=idfc2$logFC)
rownames(genelist)=genelist[,1]
circ<-circle_dat(go2,genelist)
circ<-na.omit(circ)
circ$logFC<-as.numeric(circ$logFC)
go_cir<-GOCircle(circ,rad1 = 2.5,rad2 = 3.5,label.size = 4,nsub = 10)
ggsave('02.GO_cir.png',go_cir,width =12,height = 7)
ggsave('02.GO_cir.pdf',go_cir,width =11,height = 6)

## KEGG富集分析（气泡图）
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "hsa",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 1)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk@result,file = "KEGG.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk, showCategory=15,label_format=50)
kk_dot
#ggsave(filename = '03.KEGG_dot.pdf',kk_dot,w=7,h=4)
#ggsave(filename = '03.KEGG_dot.png',kk_dot,w=7,h=4)
## KEGG 弦图
kegg_cir<-cnetplot(kk,circular = TRUE, colorEdge = TRUE,showCategory = 5,
                   color_category = "#E5C494",color_gene = "#B3B3B3",
)
kegg_cir

ggsave(filename = '03.KEGG_cir.pdf',kegg_cir,w=10,h=6)
ggsave(filename = '03.KEGG_cir.png',kegg_cir,w=10,h=6)

