rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-149(modify)/")
if (! dir.exists("./01_validation")){
  dir.create("./01_validation")
}
setwd("./01_validation")

library(magrittr)
library(stringr)
library(lance)
library(DESeq2)
library(org.Hs.eg.db)
library(GSVA)
library(msigdbr)
library(ggplot2)

##炎症评分------

df = read.delim2("../00_rawdata/dat(GSE130970).xls", row.names = 1) %>% lc.tableToNum
df <- log2(df)
df.m = msigdbr()
df.m = subset(df.m, gs_name == "HALLMARK_INFLAMMATORY_RESPONSE")
gene.inflm = df.m$gene_symbol %>% as.character()
score.inflm = gsva(as.matrix(df), list(gene.inflm), kcdf = "Poisson")
score.inflm = t(score.inflm) %>% as.data.frame()
score.inflm$sample = rownames(score.inflm)
score.inflm$group = ifelse(score.inflm$V1 > median(score.inflm$V1), "high", "low")
score.inflm$group = factor(score.inflm$group, levels = c("low","high"))


### 表达水平-------
library(magrittr)
library(stringr)
library(lance)
library(reshape2)

#df.deg = read.delim2("/data/nas1/yangly/Project/BJTC149/02.DEG/01.DEG.xls", row.names = 1)
hub.gene = read.delim2("/data/nas1/yangly/Project/BJTC149/hub.gene", header = F)[[1]]
#hub.gene = AnnotationDbi::select(org.Hs.eg.db, hub.gene, "ENTREZID", "SYMBOL")

df.exp = read.delim2("../00_rawdata/dat(GSE130970).xls", row.names = 1) %>% lc.tableToNum()
df.exp <- log2(df.exp)
df.exp = df.exp[hub.gene,]
df.exp$gene = hub.gene
df.exp <- na.omit(df.exp)
df.exp.plot = melt(df.exp, id.vars = "gene", variable.name = "sample", value.name = "value")

df.pheno = score.inflm[,-1]
write.table(df.pheno,file = 'group.xls',sep = '\t',row.names = F,quote = F)
df.exp.plot$group = df.pheno$group[match(df.exp.plot$sample, df.pheno$sample)]
df.exp.plot$group = factor(df.exp.plot$group, levels = c("low", "high"))
library(rstatix)
stat.test<-df.exp.plot%>%
  group_by(gene)%>%
  wilcox_test(value ~ group)%>%
  adjust_pvalue(method = 'fdr')
write.table(stat.test,file = 'wilcox.validation.xls',sep = '\t',row.names = F,quote = F)
stat.test$p <- round(stat.test$p,digits = 3)%>%paste0("p = ", .)
stat.test$label = paste0(stat.test$gene,"\n",stat.test$p)
# df.diff = read.delim2("02.DEG/01.DEG.xls", row.names = 1)
# df.diff = df.diff[as.character(hub.genes$GeneID),]
# df.diff = lc.tableToNum(df.diff)
# df.diff$p = df.diff$padj %>% signif(3) %>% paste0("p = ", .)
# df.diff$label = paste0(rownames(df.diff),"\n",df.diff$p)

df.exp.plot = lc.tableToNum(df.exp.plot)
df.exp.plot$label = stat.test$label[match(df.exp.plot$gene, stat.test$gene)]
df.exp.plot$group = factor(df.exp.plot$group, levels = c("low", "high"))

ggplot() +
  geom_boxplot(data = df.exp.plot, aes(x = group, y = value, color = group)) + 
  facet_wrap(label ~ ., scales = "free", ncol = 4) +
  theme_bw() + 
  xlab("Inflammation Score") + 
  ylab("Expression(log2TPM)") + 
  scale_color_manual(values = c("darkgreen", "darkorange"), guide = guide_legend(ncol = 1), labels = c("Low","High")) +
  theme(text = element_text(size = 16), legend.title = element_blank(), legend.position = "right", axis.text = element_blank())
ggsave("Validation.png", width = 10, height = 10, units = "in", dpi = 300)
ggsave("Validation.pdf", width = 10, height = 10, units = "in", dpi = 300)
