rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-345/")
if (! dir.exists("./03_DELRG")){
  dir.create("./03_DELRG")
}
setwd("./03_DELRG")
library(readxl)
library(lance)
## 脂代谢
lipid<-read_xlsx('lipid.xlsx')
lipid<-lipid[!duplicated(lipid$Geneset),]
## 861
modgene<-read.delim2('/data/nas1/luchunlin/project/BJTC-345/02_WGCNA/DEmod.xls')
DELRG<-lipid[lipid$Geneset%in%modgene$.,]%>%as.data.frame()
df<-read.delim2('/data/nas1/luchunlin/project/BJTC-345/01_DEGs/DEG_sig.xls',row.names = 1)%>%lc.tableToNum()
df<-df[DELRG$Geneset,]
##64
write.table(DELRG,file = 'DELRG.xls',sep = '\t',row.names = F,quote = F)
library(ggvenn)
mydata<-list('DEmod'=modgene$.,'Lipid metabolism'=lipid$Geneset)
pdf("DELRG.venn.pdf",w=5,h=5)
ggvenn(mydata,c('DEmod','Lipid metabolism'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()
png("DELRG.venn.png",w=500,h=500)
ggvenn(mydata,c('DEmod','Lipid metabolism'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(GOplot)

gene_transform <- bitr(DELRG$Geneset,
                       fromType = "SYMBOL",
                       toType = c("ENTREZID",'ENSEMBL'),
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
go_bar<-barplot(ego, showCategory=5, split="ONTOLOGY",label_format = 50) +
  facet_grid(ONTOLOGY ~ ., scales = "free")
go_bar
## GO 圈图
go_result<-read.table('GO.xls',header = T,sep = '\t',check.names = F)
go2=data.frame(Category='ALL',ID=go_result$ID,Term=go_result$Description,Genes=gsub("/", ", ", go_result$geneID), adj_pval = go_result$p.adjust)
idfc2<-df
genelist<-data.frame(ID=rownames(idfc2),logFC=idfc2$logFC)
rownames(genelist)=genelist[,1]
circ<-circle_dat(go2,genelist)
go_cir<-GOCircle(circ,rad1 = 2.5,rad2 = 3.5,label.size = 4,nsub = 10)
ggsave('GO_cir.png',go_cir,width =9,height = 6)
ggsave('GO_cir.pdf',go_cir,width =9,height = 6)
## KEGG富集分析（气泡图）
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "hsa",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk@result,file = "KEGG.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk, showCategory=15,label_format = 50)
kk_dot
kk2<-data.frame(Category = "ALL",ID = kk$ID,Term = kk$Description, Genes = gsub("/", ", ", kk$geneID), adj_pval = kk$p.adjust)
## 读取logFC文件
circ<-circle_dat(kk2,genelist)
#termNum<-20
#geneNum<-nrow(genelist)
chord<-chord_dat(circ,genelist)
trace(GOChord, edit = TRUE)
# legend.box = "vertical",
# legend.direction = "horizontal"
# legend.key.size / height / width
kegg_chord<-GOChord(chord,
                    gene.order = 'logFC',
                    gene.space = 0.25,
                    gene.size = 5,
                    space = 0.01,
                    lfc.col = c('red','white','blue'),
                    process.label = 7)
kegg_chord


