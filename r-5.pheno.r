rm(list = ls())
library(lance)
library(magrittr)
library(stringr)
library(ggplot2)
library(org.Hs.eg.db)
library(ggpubr)

df.risk = read.delim2("full.table.tsv")
df.pheno = read.delim2("clean.data/pheno.tsv")
df.full = merge(df.risk, df.pheno, by = "sample") %>% lc.tableToNum()

df.full$t = str_match(df.full$pathologic_T, "T[1-4]") %>% as.vector
dict.t = c("T1"="T1-T2","T2"="T1-T2","T3"="T3-T4","T4"="T3-T4")
df.full$tt = dict.t[df.full$t]
ggplot(data = subset(df.full, !is.na(t)), aes(x = t, y = risk)) + 
  geom_violin(aes(fill = t)) + 
  geom_boxplot(width = 0.3) + 
  xlab(NULL) + ylab("Risk score") +
  stat_compare_means(size = 5, label.x = 2.5, label.y = 1.2, hjust = 0.5) +
  scale_fill_manual(values = grDevices::hcl.colors(4, palette = "pastel 1", rev = T)) +
  theme_classic() +
  theme(panel.grid.major.x = element_line(color = "grey", linetype = 3), 
        panel.border = element_rect(color = "black", fill = "transparent"), 
        legend.position = "none", axis.text = element_text(size = 12), axis.title = element_text(size = 14))
ggsave("04.Phenotype/01.Tstage.png", width = 5, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("04.Phenotype/01.Tstage.pdf", width = 5, height = 5, units = "in", dpi = 300, bg = "white")

df.full$n = str_match(df.full$pathologic_N, "N[0-3]")
dict.n = c("N0"="N0","N1"="N1-N3","N2"="N1-N3","N3"="N1-N3")
df.full$nn = dict.n[df.full$n]
ggplot(data = subset(df.full, !is.na(nn)), aes(x = nn, y = risk)) + 
  geom_violin(aes(fill = nn)) + 
  geom_boxplot(width = 0.3) + 
  xlab(NULL) + ylab("Risk score") +
  stat_compare_means(size = 5, label.x = 1.5, label.y = 1.2, hjust = 0.5) +
  scale_fill_manual(values = grDevices::hcl.colors(4, palette = "pastel 1", rev = T)) +
  theme_classic() +
  theme(panel.grid.major.x = element_line(color = "grey", linetype = 3), 
        panel.border = element_rect(color = "black", fill = "transparent"), 
        legend.position = "none", axis.text = element_text(size = 12), axis.title = element_text(size = 14))
ggsave("04.Phenotype/02.Nstage.png", width = 5, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("04.Phenotype/02.Nstage.pdf", width = 5, height = 5, units = "in", dpi = 300, bg = "white")

df.full$m = str_match(df.full$pathologic_M, "M[0-1]")
ggplot(data = subset(df.full, !is.na(m)), aes(x = m, y = risk)) + 
  geom_violin(aes(fill = m)) + 
  geom_boxplot(width = 0.3) + 
  xlab(NULL) + ylab("Risk score") +
  stat_compare_means(size = 5, label.x = 1.5, label.y = 0.8, hjust = 0.5) +
  scale_fill_manual(values = grDevices::hcl.colors(4, palette = "pastel 1", rev = T)) +
  theme_classic() +
  theme(panel.grid.major.x = element_line(color = "grey", linetype = 3), 
        panel.border = element_rect(color = "black", fill = "transparent"), 
        legend.position = "none", axis.text = element_text(size = 12), axis.title = element_text(size = 14))
ggsave("04.Phenotype/03.Mstage.png", width = 5, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("04.Phenotype/03.Mstage.pdf", width = 5, height = 5, units = "in", dpi = 300, bg = "white")

ggplot(data = subset(df.full, !is.na(gender.demographic)), aes(x = str_to_sentence(gender.demographic), y = risk)) + 
  geom_violin(aes(fill = str_to_sentence(gender.demographic))) + 
  geom_boxplot(width = 0.3) + 
  xlab(NULL) + ylab("Risk score") +
  stat_compare_means(size = 5, label.x = 1.5, label.y = 0.8, hjust = 0.5) +
  scale_fill_manual(values = grDevices::hcl.colors(4, palette = "pastel 1", rev = T)) +
  theme_classic() +
  theme(panel.grid.major.x = element_line(color = "grey", linetype = 3), 
        panel.border = element_rect(color = "black", fill = "transparent"), 
        legend.position = "none", axis.text = element_text(size = 12), axis.title = element_text(size = 14))
ggsave("04.Phenotype/04.Gender.png", width = 5, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("04.Phenotype/04.Gender.pdf", width = 5, height = 5, units = "in", dpi = 300, bg = "white")

rm(list = ls())
library(estimate)
df.exp = read.delim2("clean.data/fpkm.tsv")
write.table(df.exp, file = "expression4estimate.tsv", sep = '\t', quote = F, row.names = F)
filterCommonGenes(input.f = "expression4estimate.tsv", output.f = "filtered4estimate.tsv", id="EntrezID")
estimateScore(input.ds = "filtered4estimate.tsv", output.ds = "estimate.tsv", platform="illumina")
df.estimate = read.delim2("estimate.tsv", skip = 2, row.names = 1)
df.estimate = df.estimate[-1] %>% t %>% as.data.frame %>% lc.tableToNum
df.estimate$sample = rownames(df.estimate)
df.full = read.delim2("full.table.tsv")
df.full = merge(df.full, df.estimate, by = "sample", all = F)
df.full = lc.tableToNum(df.full)

library(ggstatsplot)
cor.immune = cor.test(df.full$ImmuneScore, df.full$risk, method = "spearman")
cor.immune.output = paste0("Spearman Correlation Coefficient: ", round(cor.immune$estimate,3), 
                           ", p < 0.0001")#= ", signif(cor.immune$p.value,3))
png("04.Phenotype/05.ImmuneScore.png", width = 7, height = 7, units = "in", res = 300, bg = "white")
ggscatterstats(df.full, x = "risk", y = "ImmuneScore", type = "nonparametric", 
               results.subtitle = F, 
               ggplot.component = geom_label(aes(x = -0.8, y = 4100, label = cor.immune.output), size = 4), 
               xlab = "Risk Score")
dev.off()
pdf("04.Phenotype/05.ImmuneScore.pdf", width = 7, height = 7)
ggscatterstats(df.full, x = "risk", y = "ImmuneScore", type = "nonparametric", 
               results.subtitle = F, 
               ggplot.component = geom_label(aes(x = -0.8, y = 4100, label = cor.immune.output), size = 4), 
               xlab = "Risk Score")
dev.off()

cor.str = cor.test(df.full$StromalScore, df.full$risk, method = "spearman")
cor.str.output = paste0("Spearman Correlation Coefficient: ", round(cor.str$estimate,3), 
                        ", p < 0.0001")#= ", signif(cor.str$p.value,3))
png("04.Phenotype/06.StromalScore.png", width = 7, height = 7, units = "in", res = 300, bg = "white")
ggscatterstats(df.full, x = "risk", y = "StromalScore", type = "nonparametric", 
               results.subtitle = F, 
               ggplot.component = geom_label(aes(x = -0.8, y = 2200, label = cor.str.output), size = 4), 
               xlab = "Risk Score")
dev.off()
pdf("04.Phenotype/06.StromalScore.pdf", width = 7, height = 7)
ggscatterstats(df.full, x = "risk", y = "StromalScore", type = "nonparametric", 
               results.subtitle = F, 
               ggplot.component = geom_label(aes(x = -0.8, y = 2200, label = cor.str.output), size = 4), 
               xlab = "Risk Score")
dev.off()

rm(list = ls())
library(IOBR)
df.exp = read.delim2("clean.data/fpkm.tsv")
df.pheno = read.delim2("full.table.tsv")

gene.annot = AnnotationDbi::select(org.Hs.eg.db, as.character(df.exp$gene), "SYMBOL", "ENTREZID")
gene.annot = na.omit(gene.annot)
gene.annot = subset(gene.annot, !duplicated(gene.annot$SYMBOL))
df.exp$gene = gene.annot$SYMBOL[match(df.exp$gene, gene.annot$ENTREZID)]
df.exp = na.omit(df.exp)
rownames(df.exp) = df.exp[[1]]; df.exp = df.exp[-1]

df.im = IOBR::deconvo_cibersort(df.exp, absolute = F, perm = 1, arrays = F)
df.im = as.data.frame(df.im)
rownames(df.im) = df.im[[1]]
df.im = df.im[2:23]
colnames(df.im) = str_remove(colnames(df.im), "_CIBERSORT")
df.im = t(df.im) %>% as.data.frame()
df.im = cbind(im = rownames(df.im), df.im)
write.table(df.im, "cibersort.tsv", row.names = T, col.names = T, sep = "\t", quote = F)

df.pheno = read.delim2("full.table.tsv")
library(reshape2)
df.im.plot = melt(df.im, id.vars = "im", variable.name = "sample", value.name = "value")
df.im.plot$group = df.pheno$risk.group[match(df.im.plot$sample, df.pheno$sample)]
df.im.plot = na.omit(df.im.plot)
df.im.plot$group = factor(df.im.plot$group, levels = c("Low", "High"), labels = c("Low risk", 'High risk'))
df.pheno = df.pheno[order(df.pheno$risk, decreasing = F),]
df.im.plot$sample = factor(df.im.plot$sample, levels = df.pheno$sample)

library(ggpubr)
ggplot(df.im.plot, aes(x = sample, y = value, fill = im)) + 
  geom_bar(stat = "identity", position = position_fill()) +
  facet_wrap(.~group, scales = "free_x") +
  xlab("Sample") + ylab("Relative Abundance") +
  scale_fill_manual(values = IOBR::palettes(category = "random", show_col = F)) +
  scale_y_continuous(expand = c(0,0)) +
  guides(fill = guide_legend(title = "Immune Cells", ncol = 7)) +
  theme_pubclean(base_size = 12) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave("04.Phenotype/07.CIBERSORT.Stack.png", width = 14, height = 8, units = "in", dpi = 300)
ggsave("04.Phenotype/07.CIBERSORT.Stack.pdf", width = 14, height = 8, units = "in", dpi = 300)

df.im.test = df.im[-1]
df.pheno = df.pheno[order(df.pheno$risk.group, decreasing = T),]
rownames(df.pheno) = df.pheno$sample
df.im.test = df.im.test[df.pheno$sample]

res.im = lc.wilcox.test(df.im.test, df.pheno$risk.group, p.adj.intra = "bonferroni")$Result
stat.im = lc.groupRowMedian(df.im.test, df.pheno$risk.group)
res.im = cbind(stat.im, res.im[4:5])
colnames(res.im) = c("High-risk","Low-risk","pval","padj")
write.table(res.im, "04.Phenotype/08.CIBERSORT.Diff.xls", quote = F, sep = "\t", row.names = T, col.names = T)

res.im = subset(res.im, padj < 0.05)
res.im$im = rownames(res.im)
res.im$p = res.im$padj %>% signif(3) %>% paste0("p = ", .)
df.im.box = df.im[rownames(res.im),]
df.im.box = melt(df.im.box, id.vars = "im", variable.name = "sample", value.name = "val")
df.im.box$group = df.pheno$risk.group[match(df.im.box$sample, df.pheno$sample)]
df.im.box = na.omit(df.im.box)
df.im.box$p = res.im$p[match(df.im.box$im, res.im$im)]
df.im.box$label = paste(df.im.box$im, df.im.box$p, sep = "\n")

ggplot() + 
  geom_boxplot(data = df.im.box, mapping = aes(x = im, y = val, color = group), width = 0.5) + 
  xlab(NULL) + scale_x_discrete(expand = c(0.15,0.15)) +
  ylab("Score") + 
  scale_color_manual(values = c("darkred","darkblue"), labels = c("High risk", "Low risk"), 
                     guide = guide_legend(title = NULL)) +
  facet_wrap(.~label, scales = "free") +
  theme_pubclean(base_size = 12) + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
ggsave("04.Phenotype/09.CIBERSORT.Box.png", width = 8, height = 7, units = "in", dpi = 300)
ggsave("04.Phenotype/09.CIBERSORT.Box.pdf", width = 8, height = 7, units = "in", dpi = 300)

rm(list = ls())
library(lance)
library(magrittr)
library(stringr)
library(ggplot2)
library(org.Hs.eg.db)
library(GSVA)
library(msigdbr)

df.exp = read.delim2("clean.data/fpkm.tsv", row.names = 1) %>% lc.tableToNum()
msigdb = msigdbr()
msigdb.kegg = subset(msigdb, gs_subcat == "CP:KEGG")
msigdb.bp = subset(msigdb, gs_subcat =="GO:BP")
msigdb.mf = subset(msigdb, gs_subcat == "GO:MF")
msigdb.cc = subset(msigdb, gs_subcat == "GO:CC")

gs.kegg = split(msigdb.kegg$entrez_gene, msigdb.kegg$gs_name)
gs.bp = split(msigdb.bp$entrez_gene, msigdb.bp$gs_name)
gs.cc = split(msigdb.cc$entrez_gene, msigdb.cc$gs_name)
gs.mf = split(msigdb.mf$entrez_gene, msigdb.mf$gs_name)

df.es.kegg = gsva(as.matrix(df.exp), gs.kegg, kcdf = "Poisson")
df.es.bp = gsva(as.matrix(df.exp), gs.bp, kcdf = "Poisson")
df.es.cc = gsva(as.matrix(df.exp), gs.cc, kcdf = "Poisson")
df.es.mf = gsva(as.matrix(df.exp), gs.mf, kcdf = "Poisson")

df.es.kegg = as.data.frame(df.es.kegg)
df.es.bp = as.data.frame(df.es.bp)
df.es.cc = as.data.frame(df.es.cc)
df.es.mf = as.data.frame(df.es.mf)

df.es.kegg = cbind(Pathway = rownames(df.es.kegg), df.es.kegg)
df.es.bp = cbind(BP = rownames(df.es.bp), df.es.bp)
df.es.cc = cbind(CC = rownames(df.es.cc), df.es.cc)
df.es.mf = cbind(MF = rownames(df.es.mf), df.es.mf)

df.grp = read.delim2("full.table.tsv")
df.es.kegg = df.es.kegg[df.grp$sample] %>% lc.tableToNum
df.es.bp = df.es.bp[df.grp$sample] %>% lc.tableToNum
df.es.cc = df.es.cc[df.grp$sample] %>% lc.tableToNum
df.es.mf = df.es.mf[df.grp$sample] %>% lc.tableToNum

design.mat = cbind(low = ifelse(df.grp$risk.group == "Low", 1, 0), 
                   high = ifelse(df.grp$risk.group == "Low", 0, 1))
contrast.mat = makeContrasts(contrasts="high-low", levels=design.mat)

fit.kegg = lmFit(df.es.kegg, design.mat)
fit.bp = lmFit(df.es.bp, design.mat)
fit.cc = lmFit(df.es.cc, design.mat)
fit.mf = lmFit(df.es.mf, design.mat)

fit.kegg = contrasts.fit(fit.kegg, contrast.mat)
fit.bp = contrasts.fit(fit.bp, contrast.mat)
fit.cc = contrasts.fit(fit.cc, contrast.mat)
fit.mf = contrasts.fit(fit.mf, contrast.mat)

fit.kegg = eBayes(fit.kegg)
fit.bp = eBayes(fit.bp)
fit.cc = eBayes(fit.cc)
fit.mf = eBayes(fit.mf)

fit.kegg = topTable(fit.kegg, coef = 1, number = Inf, adjust.method = "fdr")
fit.bp = topTable(fit.bp, coef = 1, number = Inf, adjust.method = "fdr")
fit.cc = topTable(fit.cc, coef = 1, number = Inf, adjust.method = "fdr")
fit.mf = topTable(fit.mf, coef = 1, number = Inf, adjust.method = "fdr")

fit.kegg = cbind(Pathway = rownames(fit.kegg), fit.kegg)
fit.bp = cbind(BP = rownames(fit.bp), fit.bp)
fit.cc = cbind(CC = rownames(fit.cc), fit.cc)
fit.mf = cbind(MF = rownames(fit.mf), fit.mf)

fit.kegg$Direction = ifelse(fit.kegg$logFC > 0, "Up", "Down")
fit.kegg$Direction = ifelse(fit.kegg$adj.P.Val < 0.05, fit.kegg$Direction, "No")
fit.bp$Direction = ifelse(fit.bp$logFC > 0, "Up", "Down")
fit.bp$Direction = ifelse(fit.bp$adj.P.Val < 0.05, fit.bp$Direction, "No")
fit.cc$Direction = ifelse(fit.cc$logFC > 0, "Up", "Down")
fit.cc$Direction = ifelse(fit.cc$adj.P.Val < 0.05, fit.cc$Direction, "No")
fit.mf$Direction = ifelse(fit.mf$logFC > 0, "Up", "Down")
fit.mf$Direction = ifelse(fit.mf$adj.P.Val < 0.05, fit.mf$Direction, "No")

fit.kegg$Pathway = str_remove(fit.kegg$Pathway, "KEGG_")
fit.bp$BP = str_remove(fit.bp$BP, "GOBP_")
fit.cc$CC = str_remove(fit.cc$CC, "GOCC_")
fit.mf$MF = str_remove(fit.mf$MF, "GOMF_")

fit.kegg = fit.kegg[order(fit.kegg$adj.P.Val, decreasing = F),];fit.kegg$label = NA
fit.bp = fit.bp[order(fit.bp$adj.P.Val, decreasing = F),];fit.bp$label = NA
fit.cc = fit.cc[order(fit.cc$adj.P.Val, decreasing = F),];fit.cc$label = NA
fit.mf = fit.mf[order(fit.mf$adj.P.Val, decreasing = F),];fit.mf$label = NA

fit.kegg = fit.kegg[order(fit.kegg$Direction, decreasing = F),]; fit.kegg$label[1:5] = fit.kegg$Pathway[1:5]
fit.kegg = fit.kegg[order(fit.kegg$Direction, decreasing = T),]; fit.kegg$label[1:5] = fit.kegg$Pathway[1:5]
fit.bp = fit.bp[order(fit.bp$Direction, decreasing = F),]; fit.bp$label[1:5] = fit.bp$BP[1:5]
fit.bp = fit.bp[order(fit.bp$Direction, decreasing = T),]; fit.bp$label[1:5] = fit.bp$BP[1:5]
fit.cc = fit.cc[order(fit.cc$Direction, decreasing = F),]; fit.cc$label[1:5] = fit.cc$CC[1:5]
fit.cc = fit.cc[order(fit.cc$Direction, decreasing = T),]; fit.cc$label[1:5] = fit.cc$CC[1:5]
fit.mf = fit.mf[order(fit.mf$Direction, decreasing = F),]; fit.mf$label[1:5] = fit.mf$MF[1:5]
fit.mf = fit.mf[order(fit.mf$Direction, decreasing = T),]; fit.mf$label[1:5] = fit.mf$MF[1:5]

library(ggrepel)
ggplot(fit.kegg, aes(x = -log10(adj.P.Val), y = logFC)) + 
  geom_point(aes(color = Direction)) + 
  geom_vline(xintercept = -log10(0.05), linetype = 2, color = "grey60") +
  geom_label_repel(aes(label = label), size = 2.5, xlim = c(40,65), direction = "y") +
  scale_color_manual(values=c("blue", "grey30", "red"), labels = c("Down", "No Change", "Up")) +
  xlim(c(0,65)) +
  theme_classic(base_size = 14) + 
  theme(axis.text = element_text(face = "bold"), legend.position = "none")
ggsave("04.Phenotype/10.GSVA.KEGG.png", width = 8, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("04.Phenotype/10.GSVA.KEGG.pdf", width = 8, height = 5, units = "in", dpi = 300, bg = "white")

fit.bp$label[4083] = "REGULATION_OF_CELLULAR_RESPONSE_TO\nMACROPHAGE_COLONY_STIMULATING_FACTOR_STIMULUS"
ggplot(fit.bp, aes(x = -log10(adj.P.Val), y = logFC)) + 
  geom_point(aes(color = Direction)) + 
  geom_vline(xintercept = -log10(0.05), linetype = 2, color = "grey60") +
  geom_label_repel(aes(label = label), size = 2.5, xlim = c(55, 110), direction = "y", hjust = 0) +
  scale_color_manual(values=c("blue", "grey30", "red"), labels = c("Down", "No Change", "Up")) +
  xlim(c(0,115)) +
  theme_classic(base_size = 14) + 
  theme(axis.text = element_text(face = "bold"), legend.position = "none")
ggsave("04.Phenotype/11.GSVA.GOBP.png", width = 8, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("04.Phenotype/11.GSVA.GOBP.pdf", width = 8, height = 5, units = "in", dpi = 300, bg = "white")

ggplot(fit.cc, aes(x = -log10(adj.P.Val), y = logFC)) + 
  geom_point(aes(color = Direction)) + 
  geom_vline(xintercept = -log10(0.05), linetype = 2, color = "grey60") +
  geom_label_repel(aes(label = label), size = 2.5, xlim = c(35, 70), direction = "y", hjust = 0) +
  scale_color_manual(values=c("blue", "grey30", "red"), labels = c("Down", "No Change", "Up")) +
  xlim(c(0,70)) +
  theme_classic(base_size = 14) + 
  theme(axis.text = element_text(face = "bold"), legend.position = "none")
ggsave("04.Phenotype/12.GSVA.GOCC.png", width = 8, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("04.Phenotype/12.GSVA.GOCC.pdf", width = 8, height = 5, units = "in", dpi = 300, bg = "white")

ggplot(fit.mf, aes(x = -log10(adj.P.Val), y = logFC)) + 
  geom_point(aes(color = Direction)) + 
  geom_vline(xintercept = -log10(0.05), linetype = 2, color = "grey60") +
  geom_label_repel(aes(label = label), size = 2.5, xlim = c(60, 100), direction = "y") +
  scale_color_manual(values=c("blue", "grey30", "red"), labels = c("Down", "No Change", "Up")) +
  xlim(c(0,100)) + ylim(-0.75,0.75) +
  theme_classic(base_size = 14) + 
  theme(axis.text = element_text(face = "bold"), legend.position = "none")
ggsave("04.Phenotype/13.GSVA.GOMF.png", width = 8, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("04.Phenotype/13.GSVA.GOMF.pdf", width = 8, height = 5, units = "in", dpi = 300, bg = "white")

rm(list=ls())
imchk = read.delim2("imchk", header = F)[[1]]
annot = AnnotationDbi::select(org.Hs.eg.db, imchk, "ENTREZID", "ALIAS")
annot = subset(annot, !duplicated(ALIAS))

df.imck = read.delim2("clean.data/fpkm.tsv", row.names = 1) %>% lc.tableToNum
df.imck = df.imck[annot$ENTREZID,] %>% add(1) %>% log2
rownames(df.imck) = annot$ALIAS

df.risk = read.delim2("full.table.tsv")
df.imck = df.imck[df.risk$sample] %>% lc.tableToNum
df.risk$risk.group = factor(df.risk$risk.group, levels = c("High", "Low"))

df.diff = lc.wilcox.test(df.imck, df.risk$risk.group, p.adj.intra = "bonferroni")
df.diff = df.diff$Result
df.diff = cbind(Gene = rownames(df.diff), df.diff)
df.diff = df.diff[c(1,5,6)]
colnames(df.diff) = c("Gene","p.val","p.adj")
write.table(df.diff, "04.Phenotype/14.ImmuneCheckpoint.xls", quote = F, col.names = T, row.names = F, sep = "\t")

df.diff = subset(df.diff, p.adj < 0.05)
df.diff = df.diff[order(df.diff$p.adj, decreasing = F),]
df.diff = df.diff[1:10,]

library(reshape2)
df4plot = cbind(Gene = rownames(df.imck), df.imck)
df4plot = df4plot[df.diff$Gene,]
df4plot = melt(df4plot, id.vars = "Gene", variable.name = "Sample", value.name = "Expr")
df4plot$Group = df.risk$risk.group[match(df4plot$Sample, df.risk$sample)]
df4plot$Group = factor(df4plot$Group, levels = c("Low","High"))
df4plot$Gene = factor(df4plot$Gene, levels = df.diff$Gene)

df.diff$p = signif(df.diff$p.adj, 2)
df.diff$p = paste0("p = ", df.diff$p)

ggplot() + 
  geom_boxplot(data = df4plot, mapping = aes(x = Gene, y = Expr, color = Group)) +
  geom_text(inherit.aes = F, mapping = aes(x = df.diff$Gene, y = 10.5, label = df.diff$p), hjust = 0.5, size = 4, angle = 60, color = "red") +
  scale_y_continuous(expand = c(0,0), limits = c(-1,12)) + 
  theme_minimal() + 
  ylab("Expression log2(FPKM+1)") +
  xlab(NULL) +
  scale_color_manual(values = c("darkblue","darkred")) +
  theme_bw() +
  theme(axis.text = element_text(size = 11),
        axis.title = element_text(size = 12), legend.text = element_text(size = 14),
        legend.title = element_blank(), legend.position = "top") 
ggsave("04.Phenotype/15.ImmuneCheckpoint.png", width = 10, height = 5, units = "in", dpi = 300, bg = "white")
ggsave("04.Phenotype/15.ImmuneCheckpoint.pdf", width = 10, height = 5, units = "in", dpi = 300, bg = "white")

