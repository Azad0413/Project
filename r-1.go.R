# 1.补充外部验证。
# 2.列线图增加c-index。
# 3.GO图改为BP/CC/MF分开做。
# 4.风险模型临床分析，以M分期为例，做成M0,M1分期的高低风险组分布

rm(list = ls())
library(lance)
library(magrittr)
library(stringr)
library(ggplot2)
library(org.Hs.eg.db)
library(DESeq2)
library(pheatmap)
setwd('/data/nas1/luchunlin/project/JNZK-204(modify)/')
df.count = read.delim2("/data/nas1/yangly/Project/JNZK204/clean.data/count.tsv", row.names = 1) %>% lc.tableToNum
grps = ifelse(colnames(df.count) %>% str_ends(".1.."),"normal","tumor")
table(grps)
#normal  tumor 
#50    374 

df.meta = data.frame(sample = colnames(df.count), group = grps)
df.meta$group = factor(df.meta$group, levels = c("normal","tumor"))
dds = DESeqDataSetFromMatrix(df.count, df.meta, design = ~group)
dds = DESeq(dds)
res.dds = results(dds) %>% as.data.frame
res.dds = na.omit(res.dds)
res.dds = res.dds[c(2,5,6)]
res.dds = cbind(EntrezID = rownames(res.dds), res.dds)
res.dds$Direction = ifelse(res.dds$log2FoldChange > 1, "Up", "No")
res.dds$Direction = ifelse(res.dds$log2FoldChange < -1, "Down", res.dds$Direction)
res.dds$Direction = ifelse(res.dds$padj < 0.05, res.dds$Direction, "No")
annot.dds = AnnotationDbi::select(org.Hs.eg.db, res.dds$EntrezID, "SYMBOL", "ENTREZID")
annot.dds = na.omit(annot.dds)
annot.dds = subset(annot.dds, !duplicated(annot.dds$ENTREZID))
res.dds$Gene = annot.dds$SYMBOL[match(res.dds$EntrezID, annot.dds$ENTREZID)]
res.dds = na.omit(res.dds)
res.dds = subset(res.dds, padj > 0)

ggplot(res.dds, aes(x = log2FoldChange, y = -log10(padj), color = Direction)) + 
  geom_point(size = 2, alpha = 0.7) + 
  scale_color_manual(values=c("blue", "grey30", "red"), labels = c("Down", "No Change", "Up")) +
  geom_hline(yintercept=-log10(0.05), col="grey60", linetype = 2) +
  geom_vline(xintercept = c(-1,1), col="grey60", linetype = 2) +
  ggtitle("Tumor vs Normal") +
  ylab("-log10(pvalue)") +
  #xlim(c(-5,10)) + ylim(c(0,150)) +
  theme_minimal(base_size = 16) + 
  theme(legend.title = element_blank(), 
        plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5),
        legend.text = element_text(size = 12), 
        axis.text = element_text(face = "bold"), legend.position = "top")
ggsave("01.DEG/01.Volcano.png", width = 10, height = 7, units = "in", dpi = 300, bg = "white")
ggsave("01.DEG/01.Volcano.pdf", width = 10, height = 7, units = "in", dpi = 300, bg = "white")

res.dds = res.dds[order(res.dds$padj, decreasing = F),]
top.genes = rbind(subset(res.dds, Direction == "Up")[1:15,],
                  subset(res.dds, Direction == "Down")[1:15,])
df.fpkm = read.delim2("/data/nas1/yangly/Project/JNZK204/clean.data/fpkm.tsv", row.names = 1)
df.meta = df.meta[order(df.meta$group, decreasing = F),]
df.fpkm = df.fpkm[df.meta$sample] %>% lc.tableToNum()
df.fpkm = df.fpkm[top.genes$EntrezID,]
rownames(df.fpkm) = top.genes$Gene
rowmat = top.genes["Direction"]
rownames(rowmat) = top.genes$Gene
colnames(df.meta)[2] = "Group"
df.meta$Group = str_to_sentence(df.meta$Group)
rownames(df.meta) = df.meta$sample

p = pheatmap(log2(df.fpkm+1), annotation_col = df.meta[-1], cluster_rows = T, cluster_cols = F,scale = "row", 
             color = colorRampPalette(c("green","green","green","black","red","red","red"))(100),
             annotation_row = rowmat, annotation_names_col = F, annotation_names_row = F,
             annotation_colors = list(Group = c(Normal = "darkgreen", Tumor = "darkorange"),
                                      Direction = c(Up = "red", Down = "blue")),
             show_colnames = F)
ggsave("01.DEG/02.Heatmap.png", p, width = 12, height = 6, units = "in", dpi = 300, bg = "white")
ggsave("01.DEG/02.Heatmap.pdf", p, width = 12, height = 6, units = "in", dpi = 300, bg = "white")

res.dds.out = res.dds[c(1,6,2:5)]
res.dds.out = subset(res.dds.out, Direction != "No")
write.table(res.dds.out, "01.DEG/03.DEG.xls", sep = "\t", quote = F, col.names = T, row.names = F)

rm(list = ls())
library(VennDiagram)
library(msigdbr)
df.deg = read.delim2("01.DEG/03.DEG.xls") %>% lc.tableToNum()
df.m = msigdbr()
#df.m = subset(df.m, gs_name %in% c("GOBP_LACTATE_METABOLIC_PROCESS",
#                                   "HP_INCREASED_SERUM_LACTATE",
#                                   "HP_LACTIC_ACIDOSIS",
#                                   "HP_LACTIC_ACIDURIA",
#                                   "HP_SEVERE_LACTIC_ACIDOSIS"))


#df.m = subset(df.m, gs_name %in% c(
#  "GOBP_LACTATE_METABOLIC_PROCESS",
#  "HP_INCREASED_SERUM_LACTATE",
#  "HP_LACTIC_ACIDOSIS",
#  "HP_LACTICACIDURIA",
#  "HP_SEVERE_LACTIC_ACIDOSIS",
#  "GOBP_LACTATE_TRANSMEMBRANE_TRANSPORT",
#  "GOMF_LACTATE_DEHYDROGENASE_ACTIVITY",
#  "GOMF_LACTATE_TRANSMEMBRANE_TRANSPORTER_ATER_ACTIVITY",
#  "HP_ABNORMAL_LACTATE_DEHYDROGENASE_LEVEL",
#  "HP_ELEVATED_LACTATE_PYRUVATE_RATIO",
#  "HP_INCREASED_CIRCULATING_LACTATE_DEHYDRO",
#  "HYDROGENASE_CONCENTRATION",
#  "HP_INCREASED_SERUM_LACTATE")) 
gs = read.delim2("/data/nas1/yangly/Project/JNZK204/gs.txt")[[1]]
gs = AnnotationDbi::select(org.Hs.eg.db, gs, "ENTREZID", "SYMBOL")

sel.genes = gs$ENTREZID %>% unique %>% na.omit() %>% as.character()
sel.genes = AnnotationDbi::select(org.Hs.eg.db, sel.genes, "SYMBOL", "ENTREZID")
df.sel = sel.genes
p.venn = lc.vennFromList(list(DEG = df.deg$EntrezID, Selected = df.sel$ENTREZID), label.dist = -20)
png("01.DEG/04.Venn.png", width = 5, height = 5, res = 300, units = "in", bg = "white")
grid.newpage();grid.draw(p.venn);dev.off()
pdf("01.DEG/04.Venn.pdf", width = 5, height = 5)
grid.newpage();grid.draw(p.venn);dev.off()
inter.gene = subset(df.deg, EntrezID %in% df.sel$ENTREZID)
write.table(inter.gene, "01.DEG/05.Selected.DEG.xls", sep = "\t", quote = F, col.names = T, row.names = F)

rm(list = ls())
library(clusterProfiler)
library(ggplot2)
gene.inter = read.delim2("01.DEG/05.Selected.DEG.xls")
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
library(org.Hs.eg.db)
gene.inter = read.delim2("01.DEG/05.Selected.DEG.xls")

ego <- enrichGO(gene = gene.inter$EntrezID,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
write.table(ego,file = "GO.xls",sep = "\t",quote = F,row.names = F)
# 展示富集最显著的 GO term
go_bar<-dotplot(ego, showCategory=5, split="ONTOLOGY",label_format=50) +
  facet_grid(ONTOLOGY ~ ., scales = "free")
go_bar
ggsave(filename = '01.GO_dot.pdf',go_bar,w=9,h=6)
ggsave(filename = '01.GO_dot.png',go_bar,w=9,h=6)


dotplot(res.kegg, showCategory = 10, label_format = 50)
ggsave("01.DEG/08.KEGG.png", width = 8, height = 5, units = "in", bg = "white", dpi = 300)
ggsave("01.DEG/08.KEGG.pdf", width = 8, height = 5, units = "in", bg = "white", dpi = 300)
write.table(df.kegg, "01.DEG/09.KEGG.xls", sep = "\t", quote = F, row.names = F)

