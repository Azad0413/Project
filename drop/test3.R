## 02 差异基因鉴定-----
rm(list = ls())
setwd("/data/nas1/luchunlin/project/TY0307-11")
if (! dir.exists("./01_DEG(validation3)")){
  dir.create("./01_DEGs(validation3)")
}
setwd("./01_DEGs(validation3)")
library(lance)
library(tidyverse)
library(magrittr)
library(stringr)
library(limma)
df = read.delim2("../00_rawdata/dat(GSE20681).xls", row.names = 1) %>% lc.tableToNum
df.group = read.delim2("../00_rawdata/group(GSE20681).xls")
df = df[df.group$sample]
table(df.group$group)
df.group$group = factor(df.group$group, levels = c("control", "CAD"))
design.mat = cbind(control = ifelse(df.group$group == "control", 1, 0), 
                   CAD = ifelse(df.group$group == "control", 0, 1))
contrast.mat = makeContrasts(contrasts="CAD-control", levels=design.mat)

fit = lmFit(df, design.mat)
fit = contrasts.fit(fit, contrast.mat)
fit = eBayes(fit)
fit = topTable(fit, coef = 1, number = Inf, adjust.method = "fdr")
#fit = fit[c(1,4,5)]
DEG=na.omit(fit)
logFC_cutoff <- 0
DEG$change = as.factor(
  ifelse(DEG$P.Value & abs(DEG$logFC) > logFC_cutoff,
         ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
sig_diff <- subset(DEG,
                   DEG$P.Value < 0.05 & abs(DEG$logFC) > logFC_cutoff)

dim(DEG)
dim(sig_diff)
##842   6
summary(sig_diff$change)
# DOWN  NOT   UP 
# 297  545 
write.table(DEG,file = "DEG_all(GSE20681).xls",quote = F,sep = "\t",row.names = T)
write.table(sig_diff,file = "DEG_sig((GSE20681).xls",quote = F,sep = "\t",row.names = T)

diff1 <- read.delim2('../01_DEGs/DEG_sig.xls')
up1 <- diff1[which(diff1$change=='UP'),] 
down1 <- diff1[which(diff1$change=='DOWN'),]

up2 <- sig_diff[which(sig_diff$change=='UP'),]
down2 <- sig_diff[which(sig_diff$change=='DOWN'),]

up <- data.frame(symbol=intersect(rownames(up1),rownames(up2)))
down <- data.frame(symbol=intersect(rownames(down1),rownames(down2)))
all <- rbind(up,down)
#geneset <- read.csv('../02_DEERS/GeneCards-SearchResults.csv')
#ERS <- geneset[which(geneset$Relevance.score>=7),]

#intersce <- intersect(all$symbol,ERS$Gene.Symbol)

hub_gene <- all$symbol

train_expr <- read.delim2("../00_rawdata/dat(GSE113079).xls",row.names = 1)%>%lc.tableToNum()
train_group <- read.delim2('../00_rawdata/group(GSE113079).xls')
train_group <- train_group[order(train_group$group),]
expr <- train_expr[hub_gene, train_group$sample]

# ROC
library(geomROC)
library(ROCit)
dat5 <- expr %>% t %>% as.data.frame() %>% tibble::rownames_to_column(var = "sample")
dat5 <- merge(dat5, train_group, by.x="sample")
roc_all_formula <-  ""
n <- 1
for (i in hub_gene){
  prefix <- sprintf("%02d", n)
  roc_gene <- paste("roc_",i,sep="")
  assign(roc_gene, rocit(score=dat5[,i],class=dat5$group, method = "binormal", negref = "control"))
  roc_gene_plot <- paste("roc_", i, "_plot", sep="")
  # assign(roc_gene_plot, ggplot()+geom_roc(get(roc_gene),color='red')+
  #          roc_diagonal()+
  #          theme_bw() +
  #          annotate('text',x=0.75,y=0.25,label=paste('AUC = ', round(get(roc_gene)[["AUC"]],3),sep=""),alpha=1,size=6,family = "Times")+
  #          labs(title = i) +
  #          theme(panel.grid = element_blank(),
  #                plot.title = element_text(size = 20, face = "bold", family = "Times", hjust = 0.5),
  #                panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
  #                axis.title = element_text(size = 18, family = "Times", face = "bold"),
  #                axis.text = element_text(face = "bold", size = 12),
  #                legend.position = "none",
  #                strip.text = element_text(colour = 'black', family = "Times", face = 'bold', size = rel(1)),
  #                strip.background = element_rect(fill = "white"))
  # )
  # ggsave(filename = paste0(prefix,".",i,"_roc.png"), width = 5, height = 5, get(roc_gene_plot))
  # ggsave(filename = paste0(prefix,".",i,"_roc.pdf"), width = 5, height = 5, get(roc_gene_plot))
  # roc_all_formula <- paste(roc_all_formula, roc_gene_plot, sep="+")
  n <- n + 1
}
# 
# library(patchwork)
# roc_all <- roc_TP53_plot+roc_STAT3_plot+roc_CAV1_plot+roc_TGFBR1_plot+roc_HIF1A_plot+roc_PRKCA_plot+roc_PML_plot+roc_RB1_plot+roc_HELLS_plot+roc_RRM2_plot+plot_layout(ncol =5) 
# roc_all
# ggsave(filename = "11.all_gene_roc.pdf", height = 10, width = 25, roc_all)
# ggsave(filename = "11.all_gene_roc.png", height = 10, width = 25, roc_all)

# ref: https://www.sohu.com/a/235533810_788208
# 检查下是否为病例组（positive） vs 对照组（negative）
# 提取roc >= 0.7 的基因
pass_hub <- c()
for (i in hub_gene){
  result <- summary(get(paste("roc", i, sep="_")))
  print(result)
  auc_num <- as.numeric(strsplit(as.character(result[4,1]), " ")[[1]][4])
  if (auc_num >= 0.7){
    pass_hub <- c(pass_hub, i)
  }
}
length(pass_hub)
pass_hub
# [1] "CLN6"   "PPARG"  "CAPN3"  "ERLIN1" "VCPIP1" "MGAT2"  "USE1"   "STX5"   "ATR"    "YIPF5"  "ATL1"  
# write.table(pass_hub, file = "12.ROCpass_gene.txt", quote = F, row.names = F, col.names = F)


test_expr <- read.delim2("../00_rawdata/dat(GSE20681).xls",row.names = 1)%>%lc.tableToNum()
test_group <- read.delim2('../00_rawdata/group(GSE20681).xls')
table(test_group$group)
test_group <- test_group[order(test_group$group),]
expr <- test_expr[hub_gene, test_group$sample]

# ROC
library(geomROC)
library(ROCit)
dat5 <- expr %>% t %>% as.data.frame() %>% tibble::rownames_to_column(var = "sample")
dat5 <- merge(dat5, test_group, by.x="sample")
roc_all_formula <-  ""
n <- 1
for (i in hub_gene){
  prefix <- sprintf("%02d", n)
  roc_gene <- paste("roc_",i,sep="")
  assign(roc_gene, rocit(score=dat5[,i],class=dat5$group, method = "binormal", negref = "control"))
  roc_gene_plot <- paste("roc_", i, "_plot", sep="")
  # assign(roc_gene_plot, ggplot()+geom_roc(get(roc_gene),color='red')+
  #          roc_diagonal()+
  #          theme_bw() +
  #          annotate('text',x=0.75,y=0.25,label=paste('AUC = ', round(get(roc_gene)[["AUC"]],3),sep=""),alpha=1,size=6,family = "Times")+
  #          labs(title = i) +
  #          theme(panel.grid = element_blank(),
  #                plot.title = element_text(size = 20, face = "bold", family = "Times", hjust = 0.5),
  #                panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"),
  #                axis.title = element_text(size = 18, family = "Times", face = "bold"),
  #                axis.text = element_text(face = "bold", size = 12),
  #                legend.position = "none",
  #                strip.text = element_text(colour = 'black', family = "Times", face = 'bold', size = rel(1)),
  #                strip.background = element_rect(fill = "white"))
  # )
  # ggsave(filename = paste0(prefix,".",i,"_roc.png"), width = 5, height = 5, get(roc_gene_plot))
  # ggsave(filename = paste0(prefix,".",i,"_roc.pdf"), width = 5, height = 5, get(roc_gene_plot))
  # roc_all_formula <- paste(roc_all_formula, roc_gene_plot, sep="+")
  n <- n + 1
}

# library(patchwork)
# roc_all <- roc_TP53_plot+roc_STAT3_plot+roc_CAV1_plot+roc_TGFBR1_plot+roc_HIF1A_plot+roc_PRKCA_plot+roc_PML_plot+roc_RB1_plot+roc_HELLS_plot+roc_RRM2_plot+plot_layout(ncol =5) 
# roc_all
# ggsave(filename = "11.all_gene_roc.pdf", height = 10, width = 25, roc_all)
# ggsave(filename = "11.all_gene_roc.png", height = 10, width = 25, roc_all)

# ref: https://www.sohu.com/a/235533810_788208
# 检查下是否为病例组（positive） vs 对照组（negative）
# 提取roc >= 0.7 的基因
pass_hub2 <- c()
for (i in hub_gene){
  result <- summary(get(paste("roc", i, sep="_")))
  print(result)
  auc_num <- as.numeric(strsplit(as.character(result[4,1]), " ")[[1]][4])
  if (auc_num >= 0.7){
    pass_hub2 <- c(pass_hub2, i)
  }
}
length(pass_hub2)
# write.table(pass_hub, file = "12.ROCpass_gene.txt", quote = F, row.names = F, col.names = F)
