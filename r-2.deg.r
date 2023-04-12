rm(list = ls())
library(magrittr)
library(stringr)
library(lance)
library(limma)
library(org.Hs.eg.db)

df = read.delim2("clean.data/01.Expression.xls", row.names = 1) %>% lc.tableToNum
df.pheno = read.delim2("clean.data/02.MetaInfo.xls")
df = df[df.pheno$geo_accession]
df.pheno$group = ifelse(str_detect(df.pheno$title, "Normal"), "Control", "Tumor")
design.mat = cbind(Control = ifelse(df.pheno$group == "Control", 1, 0), 
                   Patient = ifelse(df.pheno$group == "Control", 0, 1))
contrast.mat = makeContrasts(contrasts="Patient-Control", levels=design.mat)

fit = lmFit(df, design.mat)
fit = contrasts.fit(fit, contrast.mat)
fit = eBayes(fit)
fit = topTable(fit, coef = 1, number = Inf, adjust.method = "fdr")
fit = fit[c(1,3,4)]
fit.annot = AnnotationDbi::select(org.Hs.eg.db, keys = rownames(fit), keytype = "ENTREZID", columns = "SYMBOL")
fit.annot = subset(fit.annot, !duplicated(fit.annot$ENTREZID))
fit = cbind(GeneID = fit.annot$ENTREZID, Gene = fit.annot$SYMBOL, fit)
fit = na.omit(fit)
fit$Direction = ifelse(fit$logFC > 0, "Up", "No")
fit$Direction = ifelse(fit$logFC < 0, "Down", fit$Direction)
fit$Direction = ifelse(fit$adj.P.Val < 0.05, fit$Direction, "No")
fit$Direction = factor(fit$Direction, levels = c("Down", "No", "Up"))
fit = fit[order(fit$adj.P.Val, decreasing = F), ]
fit.out = subset(fit, Direction != "No")
write.table(fit.out, "01.DEG/01.DEG.xls", sep = "\t", quote = F, col.names = T, row.names = F)

library(ggplot2)
ggplot(fit, aes(x = logFC, y = -log10(adj.P.Val), color = Direction)) + 
  geom_point(size = 1, alpha = 0.7, shape = 20) + 
  scale_color_manual(values=c("blue", "grey30", "red"), labels = c("Down", "No Change", "Up")) +
  geom_hline(yintercept=-log10(0.05), col="grey60", linetype = 2) +
  #geom_vline(xintercept = c(-1,1), col="grey60", linetype = 2) +
  ggtitle("Tumor vs Control") +
  ylab("-log10(pval)") +
  #xlim(c(-2,2)) +
  theme_bw(base_size = 16) + 
  theme(legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
        legend.text = element_text(size = 12), 
        axis.text = element_text(face = "bold"), legend.position = "bottom")
ggsave("01.DEG/02.Volcano.png", width = 10, height = 5, dpi = 300, bg = "white")
ggsave("01.DEG/02.Volcano.pdf", width = 10, height = 5, dpi = 300, bg = "white")

library(pheatmap)
fit = fit[order(fit$Direction, fit$P.Value, decreasing = F),]
fit = na.omit(fit)
fit.sub = rbind(subset(fit, Direction == "Up")[1:15,],
                subset(fit, Direction == "Down")[1:15,])
df.sub = df[rownames(fit.sub),]
rownames(df.sub) = fit.sub$Gene

col.mat = data.frame(Sample = df.pheno$geo_accession,
                     Group = ifelse(df.pheno$group =="Control", "Control", "Tumor"))
col.mat$Group = factor(col.mat$Group, levels = c("Control","Tumor"))
col.mat = col.mat[order(col.mat$Group),]
rownames(col.mat) = col.mat$Sample

row.mat = data.frame(Direction = factor(fit.sub$Direction, levels = c("Up", "Down")))
rownames(row.mat) = rownames(df.sub)
df.sub = df.sub[rownames(col.mat)]

p = pheatmap(df.sub, border_color = "black", scale = "row",
             color = colorRampPalette(c("green","green","black","red","red"))(50),
             cluster_rows = T, cluster_cols = F, show_colnames = F, 
             annotation_colors = list(Group = c(Control = "darkgreen", Tumor = "darkorange"),
                                      Direction = c(Up = "red", Down = "blue")),
             annotation_names_col = F, annotation_names_row = F,
             annotation_col = col.mat[-1], annotation_row = row.mat)

png("01.DEG/03.Heatmap.png", width = 8, height = 4, units = "in", res = 300)
print(p)
dev.off()
pdf("01.DEG/03.Heatmap.pdf", width = 8, height = 4)
print(p)
dev.off()

rm(list = ls())
gs = read.delim2("invasion.txt", header = F)[[1]]
gs = unique(gs)
deg = read.delim2("01.DEG/01.DEG.xls")

p.venn = lc.vennFromList(list(DEG = deg$Gene, Invasion = gs), label.dist = -20)
png("01.DEG/04.Venn.png", width = 5, height = 5, res = 300, units = "in", bg = "white")
grid.newpage();grid.draw(p.venn);dev.off()
pdf("01.DEG/04.Venn.pdf", width = 5, height = 5)
grid.newpage();grid.draw(p.venn);dev.off()

deg.sub = subset(deg, Gene %in% gs)
write.table(deg.sub, "01.DEG/05.Invasion.DEG.xls", quote = F, sep = "\t", row.names = F, col.names = T)

rm(list = ls())
library(ggplot2)
gene.inter = read.delim2("01.DEG/05.Invasion.DEG.xls")
direct = gene.inter$Direction %>% table
lbs = paste0(names(direct), "(", direct, ")")
png("01.DEG/06.Pie.png", width = 4, height = 4, units = "in", res = 300)
pie(direct, labels = lbs, col = c("blue","red"))
dev.off()
pdf("01.DEG/06.Pie.pdf", width = 4, height = 4)
pie(direct, labels = lbs, col = c("blue","red"))
dev.off()


rm(list = ls())
library(clusterProfiler)
library(ggplot2)
gene.inter = read.delim2("01.DEG/05.Invasion.DEG.xls")
res.kegg = enrichKEGG(as.character(gene.inter$GeneID) %>% as.character, organism = "hsa", keyType = "ncbi-geneid")
res.go = enrichGO(as.character(gene.inter$GeneID), org.Hs.eg.db, keyType = "ENTREZID", ont = "ALL")
res.kegg = setReadable(res.kegg, org.Hs.eg.db, "ENTREZID")
res.go = setReadable(res.go, org.Hs.eg.db, "ENTREZID")
df.kegg = res.kegg@result
df.go = res.go@result

dotplot(res.go, showCategory = 10, label_format = 50) 
ggsave("01.DEG/07.GO.png", width = 8, height = 5, units = "in", bg = "white", dpi = 300)
ggsave("01.DEG/07.GO.pdf", width = 8, height = 5, units = "in", bg = "white", dpi = 300)
write.table(df.go, "01.DEG/08.GO.xls", sep = "\t", quote = F, row.names = F)

dotplot(res.kegg, showCategory = 10, label_format = 50)
ggsave("01.DEG/08.KEGG.png", width = 8, height = 4, units = "in", bg = "white", dpi = 300)
ggsave("01.DEG/08.KEGG.pdf", width = 8, height = 4, units = "in", bg = "white", dpi = 300)
write.table(df.kegg, "01.DEG/09.KEGG.xls", sep = "\t", quote = F, row.names = F)








