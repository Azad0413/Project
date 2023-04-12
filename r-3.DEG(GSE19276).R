rm(list = ls())
# 02 差异分析---------
setwd("/data/nas1/luchunlin/project/BJTC-300-8/")
if (! dir.exists("./03_DEGs(GSE19276)")){
  dir.create("./03_DEGs(GSE19276)")
}
setwd("./03_DEGs(GSE19276)")
## 读取数据
library(magrittr)
library(stringr)
library(lance)
library(limma)
df = read.delim2("../00_rawdata/dat(GSE19276).xls", row.names = 1) %>% lc.tableToNum
df.group = read.delim2("../00_rawdata/group(GSE19276).xls")
df = df[df.group$sample]
df.group$group = factor(df.group$group, levels = c("control", "OS"))
design.mat = cbind(control = ifelse(df.group$group == "control", 1, 0), 
                   OS = ifelse(df.group$group == "control", 0, 1))

contrast.mat = makeContrasts(contrasts="OS-control", levels=design.mat)
fit = lmFit(df, design.mat)
fit = contrasts.fit(fit, contrast.mat)
fit = eBayes(fit)
fit = topTable(fit, coef = 1, number = Inf, adjust.method = "fdr")
#fit = fit[c(1,4,5)]
DEG=na.omit(fit)
logFC_cutoff <- 1
DEG$change = as.factor(
  ifelse(DEG$adj.P.Val <0.05 & abs(DEG$logFC) > logFC_cutoff,
         ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
sig_diff <- subset(DEG,
                   DEG$adj.P.Val < 0.05 & abs(DEG$logFC) > logFC_cutoff)

dim(DEG)
dim(sig_diff)
## Gene Symbol
summary(sig_diff$change)
# DOWN  NOT   UP 
#  747    0  501 
write.table(DEG,file = "DEG_all(GSE19276).xls",quote = F,sep = "\t",row.names = T)
write.table(sig_diff,file = "DEG_sig(GSE19276).xls",quote = F,sep = "\t",row.names = T)



## 曼哈顿图------
#BiocManager::install('qqman')
library(qqman)
data(gwasResults)
head(gwasResults)
manhattan(gwasResults,col = c("royalblue4", "darksalmon"),suggestiveline = -log10(1e-05), genomewideline = -log10(5e-08), annotatePval = 5e-08, annotateTop = FALSE)

##ggplot2 绘制 GWAS 曼哈顿图
library(ggplot2)

#计算染色体刻度坐标
gwasResults$SNP1 <- seq(1, nrow(gwasResults), 1)
gwasResults$CHR <- factor(gwasResults$CHR, levels = unique(gwasResults$CHR))
chr <- aggregate(gwasResults$SNP1, by = list(gwasResults$CHR), FUN = median)

#ggplot2 作图
#定义 p < 1e-05 为临界显著性，p < 5e-08 为高可信显著性
p <- ggplot(gwasResults, aes(SNP1, -log(P, 10))) +
  annotate('rect', xmin = 0, xmax = max(gwasResults$SNP1), ymin = -log10(1e-05), ymax = -log10(5e-08), fill = 'gray98') +
  geom_hline(yintercept = c(-log10(1e-05), -log10(5e-08)), color = c("#F39B7FB2","#91D1C2B2"), size = 0.35) +
  geom_point(aes(color = CHR), show.legend = FALSE) +
  scale_color_manual(values = rep(c("grey", "skyblue"), 11)) +
  scale_x_continuous(breaks = chr$x, labels = chr$Group.1, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(1, 9, 2), labels = as.character(seq(1, 9, 2)), expand = c(0, 0), limits = c(0, 9)) +
  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'), panel.background = element_rect(fill = 'transparent')) +
  labs(x = 'Chromosome', y = expression(''~-log[10]~'(P)'))

p


diff1 <- read.delim2('../01_DEGs(GSE16088)/DEG_sig(GSE16088).xls')%>%select('logFC','adj.P.Val','change')

up1 <- diff1[which(diff1$change=='UP'),]
up1$anno <- 'UP(GSE16088)'
down1 <- diff1[which(diff1$change=='DOWN'),]
down1$anno <- 'DOWN(GSE16088)'
diff2 <- read.delim2('../02_DEGs(GSE99671)/DEG_sig(GSE99671).xls',row.names = 1)%>%select('log2FoldChange','padj','change')
colnames(diff2) <- colnames(diff1)
up2 <- diff2[which(diff2$change=='UP'),]
up2$anno <- 'UP(GSE99671)'
down2 <- diff2[which(diff2$change=='DOWN'),]
down2$anno <- 'DOWN(GSE99671)'
diff3 <- read.delim2('../03_DEGs(GSE19276)/DEG_sig(GSE19276).xls')%>%select('logFC','adj.P.Val','change')
up3 <- diff3[which(diff3$change=='UP'),]
up3$anno <- 'UP(GSE19276)'
down3 <- diff3[which(diff3$change=='DOWN'),]
down3$anno <- 'DOWN(GSE19276)'
man.dat <- rbind(up1,down1,up2,down2,up3,down3)

table(man.dat$anno)
man.dat <- rownames_to_column(man.dat,var = 'symbol')
man.dat$x <- seq(1,nrow(man.dat),1)
anno <- aggregate(man.dat$x, by = list(man.dat$anno), FUN = median)
man.dat$anno <- factor(man.dat$anno,levels = unique(man.dat$anno))
man.dat$adj.P.Val <- as.numeric(man.dat$adj.P.Val)
dat_rep<-head(man.dat[order(man.dat$adj.P.Val,decreasing = F),],20)
#dat_rep <- data.frame(symbol=dat_rep$symbol)
manhattan <- ggplot(man.dat, aes(x, -log10(adj.P.Val))) +
  annotate('rect', xmin = 0, xmax = max(man.dat$x), ymin = -log10(1e-05), ymax = -log10(5e-08), fill = 'gray98') +
  geom_hline(yintercept = c(-log10(0.05), -log10(5e-30)), color = c("#F39B7FB2","#91D1C2B2"), size = 0.35) +
  geom_point(aes(color = man.dat$anno), show.legend = FALSE) +
  scale_color_manual(values = rainbow(7)) +
  scale_x_continuous(breaks = anno$x, labels = anno$Group.1, expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(1, 20, 2), labels = as.character(seq(1, 20, 2)), expand = c(0, 0), limits = c(0, 20)) +
  theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'), panel.background = element_rect(fill = 'transparent'),
        axis.text.x = element_text(angle = 45,hjust = 1,size = 13)) +
  labs(x = '', y = expression(''~-log[10]~'(padj)'))+
  geom_label_repel(
    data = dat_rep,
    aes(label = dat_rep$symbol),
    max.overlaps = 20,
    size = 4,
    box.padding = unit(0.5, "lines"),
    min.segment.length = 0,
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )

manhattan
ggsave(manhattan,filename = '01.manhattan.pdf',w=9,h=7)
ggsave(manhattan,filename = '01.manhattan.png',w=9,h=7)

# test <- data.frame(SNP=man.dat$symbol,CHR=man.dat$anno,BP=man.dat$x,P=man.dat$adj.P.Val)
# table(test$CHR)
# test$CHR <- ifelse(test$CHR=='UP(GSE16088)',1,
#                    ifelse(test$CHR=='DOWN(GSE16088)',2,
#                           ifelse(test$CHR=='UP(GSE99671)',3,
#                                  ifelse(test$CHR=='DOWN(GSE99671)',4,
#                                         ifelse(test$CHR=='UP(GSE19276)',5,6)))))
# test$CHR <- as.numeric(test$CHR)
# manhattan(test, col = c("royalblue4", "darksalmon"), suggestiveline = -log10(0.05), genomewideline = -log10(5e-08), annotatePval = 5e-08, annotateTop = FALSE)
# 
# man.dat <- man.dat[order(man.dat$anno),]
# test$SNP1 <- seq(1,nrow(test),1)
# test$CHR <- factor(test$CHR, levels = unique(test$CHR))
# test <- test[order(test$SNP),]
# chr <- aggregate(test$SNP1, by = list(test$CHR), FUN = median)
# p <- ggplot(test, aes(SNP1, -log(P, 10))) +
#   annotate('rect', xmin = 0, xmax = max(test$SNP1), ymin = -log10(1e-05), ymax = -log10(5e-08), fill = 'gray98') +
#   geom_hline(yintercept = c(-log10(0.05), -log10(5e-30)), color = c("#F39B7FB2","#91D1C2B2"), size = 0.35) +
#   geom_point(aes(color = CHR), show.legend = FALSE) +
#   scale_color_manual(values = rep(c("grey", "skyblue"), 11)) +
#   scale_x_continuous(breaks = chr$x, labels = chr$Group.1, expand = c(0, 0)) +
#   scale_y_continuous(breaks = seq(1, 20, 2), labels = as.character(seq(1, 20, 2)), expand = c(0, 0), limits = c(0, 20)) +
#   theme(panel.grid = element_blank(), axis.line = element_line(colour = 'black'), panel.background = element_rect(fill = 'transparent')) +
#   labs(x = 'Chromosome', y = expression(''~-log[10]~'(P)'))
# 
# p
