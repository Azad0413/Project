rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-300-8/")
if (! dir.exists("./20_drug")){
  dir.create("./20_drug")
}
setwd("./20_drug")
dat<-read.delim2("../00_rawdata/dat.fpkm.xls", row.names = 1)%>% lc.tableToNum
colnames(dat)<-gsub('.','-',colnames(dat),fixed = T)
risk<-read.delim2('../08_risk/risk.xls')
high.sample<-risk$id[which(risk$risk==0)]
low.sample<-risk$id[which(risk$risk==1)]
dat<-dat[,risk$id]
library(pRRophetic)
library(ggplot2)
set.seed(12345)
model_geneids<-read.delim2('../04_DEARG/DEARG.xls',header = T)
model_expr<-dat[model_geneids$symbol, risk$id]
#model_expr <- log2(model_expr+1)
riskscore<-data.frame(risk$id,risk$risk)
colnames(riskscore)<-c('sample','risk')
riskscore$risk[which(riskscore$risk==1)] <-'Low risk'
riskscore$risk[which(riskscore$risk==0)] <-'High risk'
head(riskscore)
drug<-read.table(file = '/data/nas1/luchunlin/pipeline/Medicinal_Sensity/drugs.txt',sep='\t',header=F)
ic50<-data.frame(riskscore$sample)
a<-data.frame(row.names=riskscore$sample,riskscore$risk)
colnames(a)<-'risk'
cnt<-1
while (cnt < 139) {
  predictedPtype <- pRRopheticPredict(as.matrix(model_expr), drug[cnt,],selection=1)
  Tipifarnib<-data.frame(predictedPtype)
  colnames(Tipifarnib)<-drug[cnt,]
  a<-cbind(a,Tipifarnib)
  cnt = cnt + 1
}
write.table(a,'IC50.xls',sep='\t',quote=F)
b<-a
# b[b<0]<-NA
# 先写成函数的形式，方便调用
removeRowsAllNa  <- function(x){x[apply(x, 1, function(y) any(!is.na(y))),]}
removeColsAllNa  <- function(x){x[, apply(x, 2, function(y) any(!is.na(y)))]}
c<-removeColsAllNa(b)
na_flag <- apply(is.na(c), 2, sum)
x <- c[, which(na_flag == 0)]
dim(x)
medicinal_result <- t(subset(x, select = -risk)) 
high_group <- risk$id[which(risk$risk==0)]
low_group <- risk$id[which(risk$risk==1)]
pvalue = padj = log2FoldChange <- matrix(0, nrow(medicinal_result), 1)
for (i in 1:nrow(medicinal_result)){
  pvalue[i, 1] = p.value = wilcox.test(medicinal_result[i, high_group],
                                       medicinal_result[i, low_group])$p.value
  log2FoldChange[i, 1] = mean(medicinal_result[i, high_group]) - 
    mean(medicinal_result[i, low_group])
}
padj <- p.adjust(as.vector(pvalue), "fdr", n = length(pvalue))
rTable <- data.frame(log2FoldChange, 
                     pvalue, 
                     padj,
                     row.names = rownames(medicinal_result))
high_group_res <- signif(apply(medicinal_result[rownames(rTable), high_group], 
                               1,
                               median), 4)
low_group_res <- signif(apply(medicinal_result[rownames(rTable), low_group], 
                              1, 
                              median), 4)
rTable <- data.frame(high_group_res, 
                     low_group_res,
                     rTable[, c("padj", "pvalue", "log2FoldChange")])
rTable$drugs <- rownames(rTable)
rTable$sig <- ifelse(rTable$pvalue < 0.05,
                     ifelse(rTable$pvalue < 0.01, 
                            ifelse(rTable$pvalue < 0.001,
                                   ifelse(rTable$pvalue < 0.0001,
                                          paste(rTable$drugs, "****",  sep = ""),
                                          paste(rTable$drugs, "***", sep = "")),
                                   paste(rTable$drugs, "**", sep = "")),
                            paste(rTable$drugs, "*",  sep = "")), 
                     rTable$drugs)

write.table(rTable,
            file = "drugs_wilcox_test.xls",
            quote = F,
            row.names = F)
DE.drug<-rTable[which(rTable$pvalue<0.05),]
DE.drug_up <- DE.drug[which(DE.drug$log2FoldChange>0),]
DE.drug_down <- DE.drug[which(DE.drug$log2FoldChange<0),]

library(tidyr)
library(ggplot2)
library(ggpubr)
library(Ipaper)
all(rownames(rTable) == rownames(medicinal_result))
drugs_res <- data.frame(drugs=rownames(medicinal_result), medicinal_result, pvalue=rTable$pvalue)
drugs_res <- drugs_res[rownames(drugs_res)%in%rownames(DE.drug_down),]
# drugs_res <- drugs_res[rownames(rTable2),]

drugs_res <- drugs_res[which(rownames(drugs_res) %in% rownames(DE.drug)),]
violin_dat <- gather(drugs_res, key=indivs, value=score, -c("drugs","pvalue"))
violin_dat$indivs <- ifelse(gsub("\\.","-",violin_dat$indivs) %in% high_group,
                            "High risk", "Low risk") 
violin_dat$indivs <- factor(violin_dat$indivs, levels = c("High risk", "Low risk"))
head(violin_dat)
drugs_hub_boxplot1 <- ggboxplot(violin_dat, x = "indivs", y = "score",
                                color = "indivs", palette = c("#A73030FF", "#0073C2FF"),
                                add = "jitter",
                                short.panel.labs = T,
                                ggtheme = theme_bw()) +
  stat_compare_means(label = "p.signif", label.x = 1.4, vjust = 0.5)
drugs_hub_boxplot <- facet(drugs_hub_boxplot1,
                           facet.by = "drugs",
                           short.panel.labs = T,
                           panel.labs.background = list(fill = "white"),
                           ncol = 4,
                           scales = "free_y") + xlab("") + ylab("IC(50)") +
  # geom_text(data=data_text,
  #           mapping=aes(x=x,y=y,label=label),nudge_x=0.1,nudge_y=0.1)+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        # strip.background = element_blank(),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 10, face = "bold",angle = 45,hjust = 1),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        text = element_text(size = 13, face = "bold"))
drugs_hub_boxplot
ggsave(filename = "01.drugs.plot(down).pdf", height = 9, width = 9)
ggsave(filename = "01.drugs.plot(down).png", height = 9, width = 9)


drugs_res <- data.frame(drugs=rownames(medicinal_result), medicinal_result, pvalue=rTable$pvalue)
drugs_res <- drugs_res[rownames(drugs_res)%in%rownames(DE.drug_up),]
# drugs_res <- drugs_res[rownames(rTable2),]

drugs_res <- drugs_res[which(rownames(drugs_res) %in% rownames(DE.drug)),]
violin_dat <- gather(drugs_res, key=indivs, value=score, -c("drugs","pvalue"))
violin_dat$indivs <- ifelse(gsub("\\.","-",violin_dat$indivs) %in% high_group,
                            "High risk", "Low risk") 
violin_dat$indivs <- factor(violin_dat$indivs, levels = c("High risk", "Low risk"))
head(violin_dat)
drugs_hub_boxplot1 <- ggboxplot(violin_dat, x = "indivs", y = "score",
                                color = "indivs", palette = c("#A73030FF", "#0073C2FF"),
                                add = "jitter",
                                short.panel.labs = T,
                                ggtheme = theme_bw()) +
  stat_compare_means(label = "p.signif", label.x = 1.4, vjust = 0.5)
drugs_hub_boxplot <- facet(drugs_hub_boxplot1,
                           facet.by = "drugs",
                           short.panel.labs = T,
                           panel.labs.background = list(fill = "white"),
                           ncol = 7,
                           scales = "free_y") + xlab("") + ylab("IC(50)") +
  # geom_text(data=data_text,
  #           mapping=aes(x=x,y=y,label=label),nudge_x=0.1,nudge_y=0.1)+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        # strip.background = element_blank(),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 10, face = "bold",angle = 45,hjust = 1),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        text = element_text(size = 13, face = "bold"))
drugs_hub_boxplot
ggsave(filename = "02.drugs.plot(up).pdf", height = 13, width = 14)
ggsave(filename = "02.drugs.plot(up).png", height = 13, width = 14)
