rm(list = ls())
library(lance)
library(magrittr)
library(stringr)
library(org.Hs.eg.db)

## Read all tables
df.count = read.delim2("raw.data/TCGA-SKCM.htseq_counts.tsv.gz", row.names = 1)
df.fpkm = read.delim2("raw.data/TCGA-SKCM.htseq_fpkm.tsv.gz", row.names = 1)
df.pheno = read.delim2("raw.data/TCGA-SKCM.GDC_phenotype.tsv.gz")
df.survival = read.delim2("raw.data/TCGA-SKCM.survival.tsv")

df.pheno$sample = make.names(df.pheno$submitter_id.samples)
df.survival$sample = make.names(df.survival$sample)
df.pheno = subset(df.pheno, sample %in% colnames(df.count))
df.survival= subset(df.survival, sample %in% colnames(df.count))

df.count = df.count[!str_starts(rownames(df.count), "__"),]

rownames(df.count) = sapply(rownames(df.count), function(x){
  str_split(x, fixed(".")) %>% unlist %>% .[1]
})
rownames(df.fpkm) = sapply(rownames(df.fpkm), function(x){
  str_split(x, fixed(".")) %>% unlist %>% .[1]
})

df.annot = AnnotationDbi::select(org.Hs.eg.db, rownames(df.count), "ENTREZID", "ENSEMBL")
df.annot = na.omit(df.annot)
df.annot = subset(df.annot, !duplicated(ENTREZID))

df.count$gene = df.annot$ENTREZID[match(rownames(df.count), df.annot$ENSEMBL)]
df.fpkm$gene = df.annot$ENTREZID[match(rownames(df.fpkm), df.annot$ENSEMBL)]
df.count = na.omit(df.count)
df.fpkm = na.omit(df.fpkm)
df.count = subset(df.count, !duplicated(gene))
df.fpkm = subset(df.fpkm, !duplicated(gene))
rownames(df.count) = df.count$gene
rownames(df.fpkm) = df.fpkm$gene
df.count = df.count[-ncol(df.count)]
df.fpkm = df.fpkm[-ncol(df.fpkm)]

## Change back
df.count = lc.tableToNum(df.count)
df.fpkm = lc.tableToNum(df.fpkm)

df.count = apply(df.count, c(1,2), function(x){
  x %>% as.numeric %>% raise_to_power(2,.) %>% as.integer %>% subtract(1)
}) %>% as.data.frame

df.fpkm = apply(df.fpkm, c(1,2), function(x){
  x %>% as.numeric %>% raise_to_power(2,.) %>% subtract(1)
}) %>% as.data.frame

df.count = cbind(gene = rownames(df.count), df.count)
df.fpkm = cbind(gene = rownames(df.fpkm), df.fpkm)

write.table(df.count, "clean.data/count.tsv", sep = "\t", row.names = F, col.names = T, quote = F)
write.table(df.fpkm, "clean.data/fpkm.tsv", sep = "\t", row.names = F, col.names = T, quote = F)
write.table(df.survival, "clean.data/survival.tsv", sep = "\t", row.names = F, col.names = T, quote = F)
write.table(df.pheno, "clean.data/pheno.tsv", sep = "\t", row.names = F, col.names = T, quote = F)



rm(list = ls())
library(magrittr)
library(stringr)
library(lance)
library(GEOquery)
library(dplyr)

options(timeout = 100000)
Sys.setenv("VROOM_CONNECTION_SIZE"=131072 * 50)
gs = getGEO("GSE15605", destdir = ".", getGPL = T, AnnotGPL = T)
df.exp = exprs(gs[[1]]) %>% as.data.frame
df.pheno = pData(gs[[1]]) %>% as.data.frame

df.pheno = subset(df.pheno, !str_detect(title, "metastasis"))
df.exp = df.exp[df.pheno$geo_accession]

df.annot = getGEO("GPL570", destdir = ".") %>% Table
df.annot = subset(df.annot, !str_detect(ENTREZ_GENE_ID, "///"))
df.exp$gene = df.annot$ENTREZ_GENE_ID[match(rownames(df.exp), df.annot$ID)] %>% as.character
df.exp = aggregate(.~gene, FUN = max, data = df.exp)
df.exp = subset(df.exp, gene != "")

write.table(df.exp, "clean.data/01.Expression.xls", sep = "\t", col.names = T, row.names = F)
write.table(df.pheno, "clean.data/02.MetaInfo.xls", sep = "\t", col.names = T, row.names = F)



rm(list = ls())
library(lance)
library(magrittr)
library(stringr)
library(org.Hs.eg.db)

## Read all tables
df.count = read.delim2("uvm/TCGA-UVM.htseq_counts.tsv.gz", row.names = 1)
df.fpkm = read.delim2("uvm/TCGA-UVM.htseq_fpkm.tsv.gz", row.names = 1)
df.pheno = read.delim2("uvm/TCGA-UVM.GDC_phenotype.tsv.gz")
df.survival = read.delim2("uvm/TCGA-UVM.survival.tsv")

df.pheno$sample = make.names(df.pheno$submitter_id.samples)
df.survival$sample = make.names(df.survival$sample)
df.pheno = subset(df.pheno, sample %in% colnames(df.count))
df.survival= subset(df.survival, sample %in% colnames(df.count))

df.count = df.count[!str_starts(rownames(df.count), "__"),]

rownames(df.count) = sapply(rownames(df.count), function(x){
  str_split(x, fixed(".")) %>% unlist %>% .[1]
})
rownames(df.fpkm) = sapply(rownames(df.fpkm), function(x){
  str_split(x, fixed(".")) %>% unlist %>% .[1]
})

df.annot = AnnotationDbi::select(org.Hs.eg.db, rownames(df.count), "ENTREZID", "ENSEMBL")
df.annot = na.omit(df.annot)
df.annot = subset(df.annot, !duplicated(ENTREZID))

df.count$gene = df.annot$ENTREZID[match(rownames(df.count), df.annot$ENSEMBL)]
df.fpkm$gene = df.annot$ENTREZID[match(rownames(df.fpkm), df.annot$ENSEMBL)]
df.count = na.omit(df.count)
df.fpkm = na.omit(df.fpkm)
df.count = subset(df.count, !duplicated(gene))
df.fpkm = subset(df.fpkm, !duplicated(gene))
rownames(df.count) = df.count$gene
rownames(df.fpkm) = df.fpkm$gene
df.count = df.count[-ncol(df.count)]
df.fpkm = df.fpkm[-ncol(df.fpkm)]

## Change back
df.count = lc.tableToNum(df.count)
df.fpkm = lc.tableToNum(df.fpkm)

df.count = apply(df.count, c(1,2), function(x){
  x %>% as.numeric %>% raise_to_power(2,.) %>% as.integer %>% subtract(1)
}) %>% as.data.frame

df.fpkm = apply(df.fpkm, c(1,2), function(x){
  x %>% as.numeric %>% raise_to_power(2,.) %>% subtract(1)
}) %>% as.data.frame

df.count = cbind(gene = rownames(df.count), df.count)
df.fpkm = cbind(gene = rownames(df.fpkm), df.fpkm)

write.table(df.count, "uvm/count.tsv", sep = "\t", row.names = F, col.names = T, quote = F)
write.table(df.fpkm, "uvm/fpkm.tsv", sep = "\t", row.names = F, col.names = T, quote = F)
write.table(df.survival, "uvm/survival.tsv", sep = "\t", row.names = F, col.names = T, quote = F)
write.table(df.pheno, "uvm/pheno.tsv", sep = "\t", row.names = F, col.names = T, quote = F)
