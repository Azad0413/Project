rm(list = ls())
# 02 差异分析---------
setwd("/data/nas1/luchunlin/project/BJTC-300-8/")
if (! dir.exists("./01_DEGs(GSE16088)")){
  dir.create("./01_DEGs(GSE16088)")
}
setwd("./01_DEGs(GSE16088)")
## 读取数据
library(magrittr)
library(stringr)
library(lance)
library(limma)
df = read.delim2("../00_rawdata/dat(GSE16088).xls", row.names = 1) %>% lc.tableToNum
df.group = read.delim2("../00_rawdata/group(GSE16088).xls")
table(df.group$group)
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
## 5345
summary(sig_diff$change)
# DOWN  NOT   UP 
# 2440    0 2905 
write.table(DEG,file = "DEG_all(GSE16088).xls",quote = F,sep = "\t",row.names = T)
write.table(sig_diff,file = "DEG_sig(GSE16088).xls",quote = F,sep = "\t",row.names = T)


