rm(list = ls())
#08 风险评分与临床指标相关性分析-------------
setwd("/data/nas1/luchunlin/project/BJTC-258/")
if (! dir.exists("./11_GSVA")){
  dir.create("./11_GSVA")
}
setwd("./11_GSVA")
dat.tcga<-read.delim2("/data/nas1/luchunlin/project/BJTC-258/00_rawdata/dat.fpkm.xls", row.names = 1)%>% lc.tableToNum
colname<-data.frame(sample=colnames(dat.tcga))
colname$sample<-gsub('.','-',colname$sample,fixed = T)
colnames(dat.tcga)<-colname$sample
library(GSVA)
library(GSEABase)
library(limma)

group<-read.delim2('/data/nas1/luchunlin/project/BJTC-258/05_survival/clinical.xls')%>%dplyr::select(c('sample','group'))
high.sample<-group$sample[which(group$group=='High')]
group2 <- group
colnames(group)
colnames(group2)<-c('id','label')
gsva_exp<-dat.tcga[,group2$id]
all(colnames(gsva_exp) == group2$id)
dim(gsva_exp)
## 19712   378
# 分组
group_score <- group2$label %>% as.factor()
design_score <- model.matrix(~0 + group_score)
rownames(design_score) <- colnames(gsva_exp)
colnames(design_score) <- levels(group_score)
compare_score <- makeContrasts("High-Low", levels = design_score)
# 09-1 GO-BP---------
GOBP_ref <- getGmt("/data/nas1/luchunlin/pipeline/GSVA/c5.go.bp.v7.4.symbols.gmt")
es_GOBP <- gsva(as.matrix(gsva_exp), GOBP_ref,
                min.sz=10, max.sz=500, verbose=TRUE)
es_GOBP <- as.data.frame(es_GOBP)

fit <- lmFit(es_GOBP, design_score)
fit2 <- contrasts.fit(fit ,compare_score)
fit3 <- eBayes(fit2)
allGeneSets <- topTable(fit3, coef = 1, number = Inf)
logFCcutoff <- 0
allGeneSets$change = as.factor(
  ifelse(allGeneSets$adj.P.Val < 0.05 & abs(allGeneSets$logFC) > logFCcutoff,
         ifelse(allGeneSets$logFC > logFCcutoff,'UP','DOWN'),'NOT')
)

DEGeneSets <- subset(allGeneSets,
                     allGeneSets$adj.P.Val < 0.05 & abs(allGeneSets$logFC) > 0)
write.table(allGeneSets,
            file = "GOBP.xls",
            quote = F,
            sep = "\t",
            row.names = T)

DEGeneSets <- DEGeneSets[order(DEGeneSets$adj.P.Val),]
dim(DEGeneSets)
# [1] 3145    6
up_DEG <- DEGeneSets[which(DEGeneSets$change == "UP"),]
down_DEG <- DEGeneSets[which(DEGeneSets$change == "DOWN"),]
up_DEG$gobp <- rownames(up_DEG)
down_DEG$gobp <- rownames(down_DEG)
# Volcano
library(ggplot2)
library(ggthemes)
library(Ipaper)
library(ggrepel)
allGeneSets$change <- factor(allGeneSets$change, levels = c("UP", "NOT", "DOWN"))
volcano_gobp <- ggplot(data = allGeneSets,
                       aes(y = logFC,
                           x = -log10(adj.P.Val), 
                           colour = change)) +
  scale_color_manual(values = c("#A73030FF", "darkgray","#0073C2FF")) +
  geom_point(size = 1.5, alpha = 0.5, na.rm=T) +
  xlim(0, 65) +
  geom_label_repel(
    data =up_DEG[1:5,],
    aes(label = gobp),
    fontface = "italic",
    size = 3.5,
    color = "black",
    segment.color = "black", 
    show.legend = FALSE,
    direction = "y",
    hjust = 0,
    force = 0.5,
    force_pull = 0,
    # nudge_x = 5,
    box.padding = 0.2,
    max.overlaps = 20,
    xlim = c(20,55)) +
  geom_label_repel(
    data = down_DEG[1:5,],
    aes(label = gobp),
    fontface = "italic",
    size = 3.5,
    color = "black",
    segment.color = "black", show.legend = FALSE,
    direction = "y",
    # nudge_x = 5,
    hjust = 0,
    force = 1,
    force_pull = 0,
    xlim = c(20,40),
    max.overlaps = 20) +
  theme_bw(base_size = 12, base_family = "Times") +
  geom_vline(xintercept = -log10(0.05),
             lty = 4,
             col = "darkgray",
             lwd = 0.6)+
  theme(legend.position = "right",
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(face="bold",
                                   color="black",
                                   family = "Times",
                                   size=8),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face = "bold",
                                   color = "black",
                                   size = 12),
        axis.text.y = element_text(face = "bold",
                                   color = "black",
                                   size = 12),
        axis.title.x = element_text(face = "bold",
                                    color = "black",
                                    size = 12),
        axis.title.y = element_text(face = "bold",
                                    color = "black",
                                    size = 12)) +
  labs(y = "log2 (Fold Change)",
       x = "-log10 (p-adjust)",
       title = "GSVA for Biological Process in GO Database")

volcano_gobp
# 09-2 GO-CC--------
GOCC_ref <- getGmt("/data/nas1/luchunlin/pipeline/GSVA/c5.go.cc.v7.4.symbols.gmt")
es_GOCC <- gsva(as.matrix(gsva_exp), GOCC_ref,
                min.sz=10, max.sz=500, verbose=TRUE)
es_GOCC <- as.data.frame(es_GOCC)
fit <- lmFit(es_GOCC, design_score)
fit2 <- contrasts.fit(fit ,compare_score)
fit3 <- eBayes(fit2)
allGeneSets <- topTable(fit3, coef = 1, number = Inf)
logFCcutoff <- 0
allGeneSets$change = as.factor(
  ifelse(allGeneSets$adj.P.Val < 0.05 & abs(allGeneSets$logFC) > logFCcutoff,
         ifelse(allGeneSets$logFC > logFCcutoff,'UP','DOWN'),'NOT')
)

DEGeneSets <- subset(allGeneSets,
                     allGeneSets$adj.P.Val < 0.05 & abs(allGeneSets$logFC) > 0)

write.table(allGeneSets,
            file = "GOCC.xls",
            quote = F,
            sep = "\t",
            row.names = T)

DEGeneSets <- DEGeneSets[order(DEGeneSets$adj.P.Val),]
dim(DEGeneSets)
# [1] 353   6
up_DEG <- DEGeneSets[which(DEGeneSets$change == "UP"),]
down_DEG <- DEGeneSets[which(DEGeneSets$change == "DOWN"),]
up_DEG$gobp <- rownames(up_DEG)
down_DEG$gobp <- rownames(down_DEG)

# Volcano
library(ggplot2)
library(ggthemes)
library(Ipaper)
library(ggrepel)
allGeneSets$change <- factor(allGeneSets$change, levels = c("UP", "NOT", "DOWN"))
volcano_gocc <- ggplot(data = allGeneSets,
                       aes(y = logFC,
                           x = -log10(adj.P.Val), 
                           colour = change)) +
  scale_color_manual(values = c("#A73030FF", "darkgray","#0073C2FF")) +
  geom_point(size = 1.5, alpha = 0.5, na.rm=T) +
  xlim(0, 60) +
  geom_label_repel(
    data = up_DEG[1:5,],
    aes(label = gobp),
    fontface = "italic",
    size = 3.5,
    color = "black",
    segment.color = "black", 
    show.legend = FALSE,
    direction = "y",
    hjust = 0,
    force = 0.5,
    force_pull = 0,
    # nudge_x = 5,
    box.padding = 0.2,
    max.overlaps = 20,
    xlim = c(20,60)) +
  geom_label_repel(
    data = down_DEG[1:5,],
    aes(label = gobp),
    fontface = "italic",
    size = 3.5,
    color = "black",
    segment.color = "black", show.legend = FALSE,
    direction = "y",
    # nudge_x = 5,
    hjust = 0,
    force = 1,
    force_pull = 0,
    xlim = c(20,35),
    max.overlaps = 20) +
  theme_bw(base_size = 12, base_family = "Times") +
  geom_vline(xintercept = -log10(0.05),
             lty = 4,
             col = "darkgray",
             lwd = 0.6)+
  theme(legend.position = "right",
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(face="bold",
                                   color="black",
                                   family = "Times",
                                   size=8),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face = "bold",
                                   color = "black",
                                   size = 12),
        axis.text.y = element_text(face = "bold",
                                   color = "black",
                                   size = 12),
        axis.title.x = element_text(face = "bold",
                                    color = "black",
                                    size = 12),
        axis.title.y = element_text(face = "bold",
                                    color = "black",
                                    size = 12)) +
  labs(y = "log2 (Fold Change)",
       x = "-log10 (p-adjust)",
       title = "GSVA for Cellular Component in GO Database")

volcano_gocc
# 09-3 GO-MF---------
GOMF_ref <- getGmt("/data/nas1/luchunlin/pipeline/GSVA/c5.go.mf.v7.4.symbols.gmt")
es_GOMF <- gsva(as.matrix(gsva_exp), GOMF_ref,
                min.sz=10, max.sz=500, verbose=TRUE)
es_GOMF <- as.data.frame(es_GOMF)
fit <- lmFit(es_GOMF, design_score)
fit2 <- contrasts.fit(fit ,compare_score)
fit3 <- eBayes(fit2)
allGeneSets <- topTable(fit3, coef = 1, number = Inf)

logFCcutoff <- 0
allGeneSets$change = as.factor(
  ifelse(allGeneSets$adj.P.Val < 0.05 & abs(allGeneSets$logFC) > logFCcutoff,
         ifelse(allGeneSets$logFC > logFCcutoff,'UP','DOWN'),'NOT')
)

DEGeneSets <- subset(allGeneSets,
                     allGeneSets$adj.P.Val < 0.05 & abs(allGeneSets$logFC) > 0)

write.table(allGeneSets,
            file = "GOMF.xls",
            quote = F,
            sep = "\t",
            row.names = T)
DEGeneSets <- DEGeneSets[order(DEGeneSets$adj.P.Val),]
dim(DEGeneSets)
# [1] 533   6
up_DEG <- DEGeneSets[which(DEGeneSets$change == "UP"),]
down_DEG <- DEGeneSets[which(DEGeneSets$change == "DOWN"),]
up_DEG$gobp <- rownames(up_DEG)
down_DEG$gobp <- rownames(down_DEG)
# data_repel <- rbind(DEGeneSets[1:5,], DEGeneSets[1:5,])
# data_repel$gobp <- rownames(data_repel)

# Volcano
library(ggplot2)
library(ggthemes)
library(Ipaper)
library(ggrepel)
allGeneSets$change <- factor(allGeneSets$change, levels = c("UP", "NOT", "DOWN"))
volcano_gomf <- ggplot(data = allGeneSets,
                       aes(y = logFC,
                           x = -log10(adj.P.Val), 
                           colour = change)) +
  scale_color_manual(values = c("#A73030FF", "darkgray","#0073C2FF")) +
  geom_point(size = 1.5, alpha = 0.5, na.rm=T) +
  xlim(0, 65) +
  geom_label_repel(
    data = up_DEG[1:5,],
    aes(label = gobp),
    fontface = "italic",
    size = 3.5,
    color = "black",
    segment.color = "black", 
    show.legend = FALSE,
    direction = "y",
    hjust = 0,
    force = 0.5,
    force_pull = 0,
    # nudge_x = 5,
    box.padding = 0.2,
    max.overlaps = 20,
    xlim = c(20,70)) +
  geom_label_repel(
    data = down_DEG[1:5,],
    aes(label = gobp),
    fontface = "italic",
    size = 3.5,
    color = "black",
    segment.color = "black", show.legend = FALSE,
    direction = "y",
    # nudge_x = 5,
    hjust = 0,
    force = 1,
    force_pull = 0,
    xlim = c(20,30),
    max.overlaps = 20) +
  theme_bw(base_size = 12, base_family = "Times") +
  geom_vline(xintercept = -log10(0.05),
             lty = 4,
             col = "darkgray",
             lwd = 0.6)+
  theme(legend.position = "right",
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(face="bold",
                                   color="black",
                                   family = "Times",
                                   size=8),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face = "bold",
                                   color = "black",
                                   size = 12),
        axis.text.y = element_text(face = "bold",
                                   color = "black",
                                   size = 12),
        axis.title.x = element_text(face = "bold",
                                    color = "black",
                                    size = 12),
        axis.title.y = element_text(face = "bold",
                                    color = "black",
                                    size = 12)) +
  labs(y = "log2 (Fold Change)",
       x = "-log10 (p-adjust)",
       title = "GSVA for Molecular Function in GO Database")

volcano_gomf
# 09-4 KEGG------
KEGG_ref <- getGmt("/data/nas1/luchunlin/pipeline/GSVA/c2.cp.kegg.v7.4.symbols.gmt")
es_KEGG <- gsva(as.matrix(gsva_exp), KEGG_ref,
                min.sz=10, max.sz=500, verbose=TRUE)
es_KEGG <- as.data.frame(es_KEGG)
fit <- lmFit(es_KEGG, design_score)
fit2 <- contrasts.fit(fit ,compare_score)
fit3 <- eBayes(fit2)
allGeneSets <- topTable(fit3, coef = 1, number = Inf)

logFCcutoff <- 0
allGeneSets$change = as.factor(
  ifelse(allGeneSets$adj.P.Val < 0.05 & abs(allGeneSets$logFC) > logFCcutoff,
         ifelse(allGeneSets$logFC > logFCcutoff,'UP','DOWN'),'NOT')
)

DEGeneSets <- subset(allGeneSets,
                     allGeneSets$adj.P.Val < 0.05 & abs(allGeneSets$logFC) > 0)

write.table(allGeneSets,
            file = "KEGG.xls",
            quote = F,
            sep = "\t",
            row.names = T)
DEGeneSets <- DEGeneSets[order(DEGeneSets$adj.P.Val),]
dim(DEGeneSets)
# [1]103  7
up_DEG <- DEGeneSets[which(DEGeneSets$change == "UP"),]
down_DEG <- DEGeneSets[which(DEGeneSets$change == "DOWN"),]
up_DEG$gobp <- rownames(up_DEG)
down_DEG$gobp <- rownames(down_DEG)
# data_repel <- rbind(DEGeneSets[1:5,], DEGeneSets[1:5,])
# data_repel$gobp <- rownames(data_repel)

# Volcano
library(ggplot2)
library(ggthemes)
library(Ipaper)
library(ggrepel)
allGeneSets$change <- factor(allGeneSets$change, levels = c("UP", "NOT", "DOWN"))
volcano_kegg <- ggplot(data = allGeneSets,
                       aes(y = logFC,
                           x = -log10(adj.P.Val), 
                           colour = change)) +
  scale_color_manual(values = c("#A73030FF", "darkgray","#0073C2FF")) +
  geom_point(size = 1.5, alpha = 0.5, na.rm=T) +
  xlim(0, 55) +
  geom_label_repel(
    data = up_DEG[1:5,],
    aes(label = gobp),
    fontface = "italic",
    size = 3.5,
    color = "black",
    segment.color = "black", 
    show.legend = FALSE,
    direction = "y",
    hjust = 0,
    force = 0.5,
    force_pull = 0,
    # nudge_x = 5,
    box.padding = 0.2,
    max.overlaps = 20,
    xlim = c(20,55)) +
  geom_label_repel(
    data = down_DEG[1:5,],
    aes(label = gobp),
    fontface = "italic",
    size = 3.5,
    color = "black",
    segment.color = "black", show.legend = FALSE,
    direction = "y",
    # nudge_x = 5,
    hjust = 0,
    force = 1,
    force_pull = 0,
    xlim = c(20,35),
    max.overlaps = 20) +
  theme_bw(base_size = 12, base_family = "Times") +
  geom_vline(xintercept = -log10(0.05),
             lty = 4,
             col = "darkgray",
             lwd = 0.6)+
  theme(legend.position = "right",
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(face="bold",
                                   color="black",
                                   family = "Times",
                                   size=8),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face = "bold",
                                   color = "black",
                                   size = 12),
        axis.text.y = element_text(face = "bold",
                                   color = "black",
                                   size = 12),
        axis.title.x = element_text(face = "bold",
                                    color = "black",
                                    size = 12),
        axis.title.y = element_text(face = "bold",
                                    color = "black",
                                    size = 12)) +
  labs(y = "log2 (Fold Change)",
       x = "-log10 (p-adjust)",
       title = "GSVA for Pathways in KEGG Database")
volcano_kegg
