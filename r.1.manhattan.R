rm(list = ls())
setwd('/data/nas1/luchunlin/project/BJTC-307(modify)/')
library(biomaRt)
library(org.Hs.eg.db)
library(tidyverse)
library(clusterProfiler)
library(ggplot2)
diff = read.delim2('/data/nas1/yangly/Project/BJTC307/02.DEG/01.DEG_all.xls')
gs.all <- as.character(diff$GeneID)
gs.all = AnnotationDbi::select(org.Hs.eg.db, gs.all, "ENSEMBL", "ENTREZID")

gs = read.delim2("/data/nas1/yangly/Project/BJTC307/03.WGCNA/10.WGCNA_DEG.xls")$GeneID %>% as.character()
gs = AnnotationDbi::select(org.Hs.eg.db, gs, "SYMBOL", "ENTREZID")
gs = gs[!duplicated(gs$ENTREZID),]

gs.all = gs.all[!duplicated(gs.all$ENTREZID),]


mart = useMart(biomart = "ENSEMBL_MART_ENSEMBL", 
               dataset = "hsapiens_gene_ensembl",
               host = "grch37.ensembl.org")
attributes = c(
  "ensembl_gene_id",
  "chromosome_name",
  "start_position"
)
gene.pos = getBM(attributes = attributes, mart = mart)

gs.pos = subset(gene.pos, ensembl_gene_id %in% gs.all$ENSEMBL)
gs.pos$entrez = gs.all$ENTREZID[match(gs.pos$ensembl_gene_id, gs.all$ENSEMBL)]
gs.pos$Gene = AnnotationDbi::select(org.Hs.eg.db, gs.pos$entrez, "SYMBOL", "ENTREZID")$SYMBOL
diff.p <- diff[,c(2,5)]
gs.pos <- merge(gs.pos,diff.p,by='Gene')

# install.packages('qqman')
library(qqman)
colnames(gs.pos)
man.dat <- gs.pos[,c(1,3,4,6)]
man.dat$chromosome_name <- gsub('X','23',man.dat$chromosome_name)
man.dat$chromosome_name <- gsub('Y','24',man.dat$chromosome_name)
table(man.dat$chromosome_name)
class(man.dat$chromosome_name)
man.dat$chromosome_name <- as.numeric(man.dat$chromosome_name)
man.dat <- na.omit(man.dat)
colnames(man.dat) <- c('SNP','CHR','BP','P')
man.dat$P <- as.numeric(man.dat$P)
manhattan(man.dat)

snpsOfInterest <-gs$SYMBOL
pdf(file = 'manhattan.pdf',w=10,h=5)
manhattan(man.dat, col = c("#F39B7FB2","#91D1C2B2",'#FFB6C1'),genomewideline = FALSE,
          suggestiveline = -log10(0.05),highlight = snpsOfInterest,annotateTop = T,annotatePval = 0.05)

dev.off()
#准备数据,使用基础函数
# data <- man.dat
# #根据目的基因的位置，新加gene和gene_annotate列
# 
# data$gene_annotate[data$SNP%in%gs$SYMBOL] <- "yes"
# 
# colnames(data)
# library(ggrepel)
# # 绘图
# p2 <- ggplot(data, aes(x=SNP, y=-log10(P))) +
#   geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
#   scale_color_manual(values = rep(c("grey", "skyblue"), 23 )) +
#   scale_x_continuous( label = data$SNP, breaks= data$CHR ) +
#   scale_y_continuous(expand = c(0, 0) ) +  
#   geom_label_repel( data=subset(data, gene_annotate=="yes"), aes(label=SNP), size=4, col = "red") +
#   theme_bw() + 
#   theme(
#     legend.position="none",
#     panel.border = element_blank(),
#     panel.grid.major.x = element_blank(),
#     panel.grid.minor.x = element_blank()
#   )
# p2
##统计每条染色体的候选基因和全部基因数量--------
data <- man.dat
data$gene_annotate <- ifelse(data$SNP%in%gs$SYMBOL,'yes','no')

freq.dat <- data.frame(chr=rep(1:24),candidate=NA,allgenes=NA)
i <- 1
for (i in c(1:24)) {
  mydata <- data[which(data$CHR==i),]
  freq.dat$candidate[i] <-length(rownames(mydata[which(mydata$gene_annotate=='yes'),]))
  freq.dat$allgenes[i] <- length(rownames(mydata))
  i <- i+1
}
write.table(freq.dat,file = 'freq.dat.xls',sep = '\t',row.names = F,quote = F)

freq.dat$freq = freq.dat$candidate/freq.dat$allgenes

# freq.dat1 <- freq.dat[which(freq.dat$chr=='1'|freq.dat$chr=='2'),]
ka.dat <- xtabs(~freq.dat$freq+freq.dat$chr,data = freq.dat)


chisq.test(ka.dat)
fisher.test(ka.dat)

