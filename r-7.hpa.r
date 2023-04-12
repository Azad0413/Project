rm(list = ls())
library(tidyverse)
library(ggpubr)

df = read.delim2("/nas1/yangly/Database/HPA/rna_tissue_consensus.tsv")
df.coef = read.delim2("02.RiskScore/05.Coefficients.xls")
df.annot = AnnotationDbi::select(org.Hs.eg.db, df.coef$gene, "ENSEMBL", "SYMBOL")
df.organ = read.delim2("/nas1/yangly/Database/HPA/organs.txt", header = F)
df = subset(df, Gene %in% df.annot$ENSEMBL)
df = lc.tableToNum(df)
df$Tissue = str_to_sentence(df$Tissue)
df$Organ = df.organ$V1[match(df$Tissue, df.organ$V2)]
df$Organ = factor(df$Organ, levels = c("Brain", "Eye", "Endocrine tissues", 
                                       "Respiratory system", "Proximal digestive tract",
                                       "Gastrointestinal tract", "Liver and gallbladder",
                                       "Pancreas", "Kideny and urinary urinary bladder",
                                       "Male reproductive system","Breast and female reproductive system",
                                       "Muscle tissues","Connective and soft tissues","Skin",
                                       "Bone marrow and lymphoid tissues"))
df = df[order(df$Organ, decreasing = F),]
df$Tissue = factor(df$Tissue, levels = df.organ$V2)

for(g in unique(df$Gene.name)){
  d = subset(df, Gene.name == g)
  ggplot(d, aes(x = Tissue, y = nTPM, fill = Organ)) + 
    geom_bar(stat = "identity") + 
    ggtitle(g) +
    scale_fill_manual(values = c("gold","greenyellow","darkmagenta","darkcyan","pink","steelblue","thistle",
                                 "mediumaquamarine","orange","lightblue","pink","palevioletred",
                                 "turquoise","peachpuff","darkgrey")) +
    theme_pubclean(base_size = 14) +
    theme(axis.text.x.bottom = element_text(angle = 60, hjust = 1, vjust = 1), legend.position = "none")
  fn1 = paste0("02.RiskScore/10.HLA/", g, ".png")
  fn2 = paste0("02.RiskScore/10.HLA/", g, ".pdf")
  ggsave(fn1, width = 14, height = 7, units = "in", dpi = 300, bg = "white")
  ggsave(fn2, width = 14, height = 7, units = "in", dpi = 300, bg = "white")
}






