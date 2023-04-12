rm(list = ls())
setwd("/data/nas1/luchunlin/project/LZZK-519-10/")
if (! dir.exists("./05_intersect")){
  dir.create("./05_intersect")
}
setwd("./05_intersect")

sig_diff1 <- read.delim2('../01_DEGs(GSE122063)/DEG_sig(GSE122063).xls')%>%lc.tableToNum()
up1 <- sig_diff1[sig_diff1$change=='UP',]
down1 <- sig_diff1[sig_diff1$change=='DOWN',]
sig_diff2 <- read.delim2('../03_DEGs(GSE25724)/DEG_sig(GSE25724).xls')%>%lc.tableToNum()
up2 <- sig_diff2[sig_diff2$change=='UP',]
down2 <- sig_diff2[sig_diff2$change=='DOWN',]

up <- data.frame(symbol=intersect(rownames(up1),rownames(up2)))
down <- data.frame(symbol=intersect(rownames(down1),rownames(down2)))
intersect <- rbind(up,down)
#intersect <- intersect(rownames(sig_diff1),rownames(sig_diff2))%>%as.data.frame()
#colnames(intersect) <- 'symbol'

write.table(intersect,file = 'intersect.xls',sep = '\t',row.names = F,quote = F)
library(ggvenn)
mydata1<-list('UP(GSE122063)'=rownames(up1),'UP(GSE25724)'=rownames(up2))
library(ggvenn)

pdf(file = '01.up.pdf',w=6,h=6)
ggvenn(mydata1,c('UP(GSE122063)','UP(GSE25724)'),
       fill_color = c("#ffb2b2","#b2e7cb","#F0E442"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()
png(file = '01.up.png',w=400,h=400)
ggvenn(mydata1,c('UP(GSE122063)','UP(GSE25724)'),
       fill_color = c("#ffb2b2","#b2e7cb","#F0E442"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()

mydata2<-list('DOWN(GSE122063)'=rownames(down1),'DOWN(GSE25724)'=rownames(down2))
pdf(file = '02.down.pdf',w=6,h=6)
ggvenn(mydata2,c('DOWN(GSE122063)','DOWN(GSE25724)'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()
png(file = '02.down.png',w=400,h=400)
ggvenn(mydata2,c('DOWN(GSE122063)','DOWN(GSE25724)'),
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

gene_transform <- bitr(intersect$symbol,
                       fromType = "SYMBOL",
                       toType = c("ENTREZID",'ENSEMBL'),
                       OrgDb = "org.Hs.eg.db")
ego <- enrichGO(gene = gene_transform$ENTREZID,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "ALL",
                pAdjustMethod = "none",
                pvalueCutoff = 0.05,
                qvalueCutoff = 1,
                readable = TRUE)
write.table(ego,file = "GO.xls",sep = "\t",quote = F,row.names = F)
# 展示富集最显著的 GO term
go_dot<-dotplot(ego, showCategory=5, split="ONTOLOGY",label_format = 70) +
  facet_grid(ONTOLOGY ~ ., scales = "free",space = 'free')
go_dot

ggsave('01.go_dot.png',go_dot,width =8,height = 6)
ggsave('01.go_dot.pdf',go_dot,width =8,height = 6)


## KEGG富集分析（气泡图）
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "hsa",
                 pAdjustMethod = "none",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 1)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk@result,file = "KEGG.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- barplot(kk, showCategory=10,label_format = 50)
kk_dot
# ggsave('04.KEGG_bar.png',kk_dot,width =8,height = 6)
# ggsave('04.KEGG_bar.pdf',kk_dot,width =8,height = 6)


##DO富集------

library(DOSE)
enrich.do <- enrichDO(gene = gene_transform$ENTREZID,
               ont = 'DO',
               pvalueCutoff = 0.05,
               qvalueCutoff = 1,
               readable = T,
               pAdjustMethod = 'none')
library(ggnewscale)

write.table(enrich.do@result,file = 'DO.xls',sep = '\t',quote = F,row.names = F)
pdf(file = '03.DO.pdf',w=8,h=4)
cnetplot(enrich.do,showCategory = 7)
dev.off()
png(file = '03.DO.png',w=600,h=300)
cnetplot(enrich.do,showCategory = 7)
dev.off()
