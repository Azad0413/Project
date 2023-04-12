rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-406-12/")
if (! dir.exists("./04_DEMPRG")){
  dir.create("./04_DEMPRG")
}
setwd("./04_DEMPRG")

sig_diff1 <- read.delim2('../01_DEG/DEG_sig.xls',row.names = 1)%>%lc.tableToNum()
# sig_diff1 <- sig_diff1[which(abs(sig_diff1$log2FoldChange)>0.6),]

sig_diff2 <- read.delim2('../03_DEG(subtype)/DEG_sig.xls',row.names = 1)%>%lc.tableToNum()
# sig_diff2 <- sig_diff2[which(abs(sig_diff2$log2FoldChange)>0.6),]

intersect <- intersect(rownames(sig_diff1),rownames(sig_diff2))%>%as.data.frame()
##2794
scgene <- read.delim2('../16_scRNA/03_SingleR/mac.gene.xls',row.names = 1)
scgene <- scgene[which(scgene$avg_log2FC>0.7),]
intersect <- intersect(intersect$.,scgene$gene)%>%as.data.frame()
##28
colnames(intersect) <- 'symbol'

write.table(intersect,file = 'DEMPRG.xls',sep = '\t',row.names = F,quote = F)
library(ggvenn)
mydata<-list('DEGs(Tumor)'=rownames(sig_diff1),'DEGs(Subtype)'=rownames(sig_diff2),'scRNA'=scgene$gene)
pdf('01.venn.pdf',w=5,h=5)
ggvenn(mydata,c('DEGs(Tumor)','DEGs(Subtype)','scRNA'),
       fill_color = c("#ffb2b2","#b2e7cb","#F0E442"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       text_size = 4,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb","#F0E442"),
       text_color = 'black')
dev.off()
png('01.venn.png',w=400,h=400)
ggvenn(mydata,c('DEGs(Tumor)','DEGs(Subtype)','scRNA'),
       fill_color = c("#ffb2b2","#b2e7cb","#F0E442"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       text_size = 4,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb","#F0E442"),
       text_color = 'black')
dev.off()

