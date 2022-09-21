rm(list = ls())
#08 风险评分与临床指标相关性分析-------------
setwd("/data/nas1/luchunlin/project/BJTC-258/")
if (! dir.exists("./14_TMB")){
  dir.create("./14_TMB")
}
setwd("./14_TMB")
mut_ov<-read_tsv(file = 'TCGA-OV.varscan2_snv.tsv')
table(mut_ov$effect)
head(mut_ov)
maf <- TCGAmutations::tcga_load(study = "OV")
library(maftools)

synonymous<-maf@maf.silent

table(synonymous$Variant_Classification)

##同义突变保留silent
synonymous<-synonymous[synonymous$Variant_Type=='SNP'&synonymous$Variant_Classification=='Silent',]
colnames(synonymous)
synonymous<-select(synonymous,c('Tumor_Sample_Barcode','Variant_Type','Variant_Classification'))
table(synonymous$Tumor_Sample_Barcode)
synonymous<-data.frame(table(synonymous$Tumor_Sample_Barcode))
colnames(synonymous)<-c('Tumor_Sample_Barcode','Synonymous')
##非同义突变 缺失突变和无义突变
non.synoymous<-maf@variant.classification.summary
colnames(non.synoymous)
View(maf@variant.classification.summary)
non.synoymous<-select(non.synoymous,c('Tumor_Sample_Barcode','Missense_Mutation','Nonsense_Mutation','total'))
non.synoymous$Non_synoymous<-rowSums(non.synoymous[,c(2,3)])
non.synoymous<-non.synoymous[,c(1,5,4)]
all.mut<-merge(synonymous,non.synoymous,by='Tumor_Sample_Barcode')
## 加上aging score
all.mut$Tumor_Sample_Barcode<-stringr::str_sub(all.mut$Tumor_Sample_Barcode,1,16)
score <- read.delim2('/data/nas1/luchunlin/project/BJTC-258/05_survival/clinical.xls')%>%
  dplyr::select(c('sample','scores','group'))
all.mut <- merge(all.mut, score, by.x="Tumor_Sample_Barcode", by.y="sample")
all.mut<-all.mut[-15,]
library(tidyr)
library(ggplot2)
library(ggpubr)
theme_nice <- theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
                    panel.grid = element_blank()) + 
  theme(axis.text.x=element_text(family="Times",size=12,face="bold"),
        axis.text.y=element_text(family="Times",size=12,face="bold"), 
        axis.title.y=element_text(family="Times",size = 15,face="bold"),
        axis.title.x=element_text(family="Times",size = 15,face="bold"),
        plot.title = element_text(hjust = 0.5,family="Times",size = 15,face="bold"),
        legend.text = element_text(face = "bold", size = 15),
        legend.title = element_text(face = "bold", size = 15),
        legend.position = "top",
        text = element_text(family = "Times"))
colnames(all.mut)
all.mut$scores<-as.numeric(all.mut$scores)
##all-------
mutall <- ggplot(all.mut, aes(x = scores,y = total)) +
  geom_point(aes(x = scores,
                 y = total,
                 color = group),
             all.mut) +
  stat_cor(method = "spearman") +
  scale_color_manual(values = c("#357EBDFF", "grey"), name = "Group") +
  geom_smooth(method = "lm",
              se = T,
              color = "red",
              formula = 'y ~ x') +
  theme_bw() +
  theme_nice
mutall
my_comparisons = list(c("High", "Low"))
colnames(all.mut)
all_dis <- ggboxplot(all.mut,
                     x = "group",
                     y = "total",
                     fill = "group",
                     palette =c("#357EBDFF", "grey")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  theme_minimal() +
  labs(title = "", x = "", y = "") +
  theme(legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")
all_dis

library(cowplot)
p1 <- insert_yaxis_grob(mutall,all_dis,grid::unit(0.8, "in"), position = "right")
ggdraw(p1)
ggsave(p1, filename = "01.mut_all.png", width = 5, height = 4)
ggsave(p1, filename = "01.mut_all.pdf", width = 5, height = 4)

## synonymous------
colnames(all.mut)
Synonymous <- ggplot(all.mut, aes(x = scores,y = Synonymous)) +
  geom_point(aes(x = scores,
                 y = Synonymous,
                 color = group),
             all.mut) +
  stat_cor(method = "spearman") +
  scale_color_manual(values = c("#357EBDFF", "grey"), name = "Group") +
  geom_smooth(method = "lm",
              se = T,
              color = "red",
              formula = 'y ~ x') +
  theme_bw() +
  theme_nice
Synonymous
my_comparisons = list(c("High", "Low"))
colnames(all.mut)
Synonymous_dis <- ggboxplot(all.mut,
                     x = "group",
                     y = "Synonymous",
                     fill = "group",
                     palette =c("#357EBDFF", "grey")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  theme_minimal() +
  labs(title = "", x = "", y = "") +
  theme(legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")
Synonymous_dis
p2 <- insert_yaxis_grob(Synonymous,Synonymous_dis,grid::unit(0.8, "in"), position = "right")
ggdraw(p2)
ggsave(p2, filename = "02.Synonymous.png", width = 5, height = 4)
ggsave(p2, filename = "02.Synonymous.pdf", width = 5, height = 4)

## Non-synonymous------
colnames(all.mut)
Non_Synonymous <- ggplot(all.mut, aes(x = scores,y = Non_synoymous)) +
  geom_point(aes(x = scores,
                 y = Non_synoymous,
                 color = group),
             all.mut) +
  stat_cor(method = "spearman") +
  scale_color_manual(values = c("#357EBDFF", "grey"), name = "Group") +
  geom_smooth(method = "lm",
              se = T,
              color = "red",
              formula = 'y ~ x') +
  theme_bw() +
  theme_nice
Non_Synonymous
my_comparisons = list(c("High", "Low"))
colnames(all.mut)
Non_synoymous_dis <- ggboxplot(all.mut,
                            x = "group",
                            y = "Non_synoymous",
                            fill = "group",
                            palette =c("#357EBDFF", "grey")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  theme_minimal() +
  labs(title = "", x = "", y = "") +
  theme(legend.title = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        legend.position = "none")
Non_synoymous_dis
p3 <- insert_yaxis_grob(Non_Synonymous,Non_synoymous_dis,grid::unit(0.8, "in"), position = "right")
ggdraw(p3)
ggsave(p3, filename = "03.Non_Synonymous.png", width = 5, height = 4)
ggsave(p3, filename = "03.Non_Synonymous.pdf", width = 5, height = 4)

gene <- read.table("/data/nas1/luchunlin/project/BJTC-258/00_rawdata/dat.fpkm.xls",row.names = 1)
gene <- rownames(gene)

maf_sample <- data.frame(barcode = maf@clinical.data$Tumor_Sample_Barcode)
maf_sample$sample <- stringr::str_sub(maf_sample$barcode, 1, 16)
score <- read.delim2('/data/nas1/luchunlin/project/BJTC-258/05_survival/clinical.xls')%>%
  dplyr::select(c('sample','group'))
colnames(score)<-c('sample','score')
maf_sample <- merge(maf_sample, score, by = "sample")
high <- subset(maf_sample, score == "High")$barcode
low <- subset(maf_sample, score == "Low")$barcode
maf_high <- subsetMaf(maf, tsb = high, gene = gene)
maf_low <- subsetMaf(maf, tsb = low, gene = gene)
pt.vs.rt <- mafCompare(m1 = maf_low, 
                       m2 = maf_high, 
                       m1Name = 'Low', 
                       m2Name = 'High', 
                       minMut = 1)
ptresult<-pt.vs.rt$results%>%subset(pval<0.05)

pdf("04.Mutation_Compare.pdf", width = 8, height = 7)
forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.05, color = c('royalblue', 'maroon'), geneFontSize = 0.8, titleSize = 1.2, lineWidth = 1.2)
dev.off()
png("04.Mutation_Compare.png", width = 8, height = 7, units = "in", res = 300)
forestPlot(mafCompareRes = pt.vs.rt, pVal = 0.05, color = c('royalblue', 'maroon'), geneFontSize = 0.8, titleSize = 1.2, lineWidth = 1.2)
dev.off()

## 互斥性和共现性-----
diff.gene<-ptresult$Hugo_Symbol
maf_high <- subsetMaf(maf, tsb = high, gene = gene)
maf.diff<- subsetMaf(maf,tsb = maf_sample$barcode,genes = diff.gene)
pdf(file = '05.somaticInteractions.pdf',w=8,h=8)
output<-somaticInteractions(maf = maf.diff,pvalue = c(0.05,0.01),showSum = F)
dev.off()

png(file = '05.somaticInteractions.png',w=600,h=600)
output<-somaticInteractions(maf = maf.diff,pvalue = c(0.05,0.01),showSum = F)
dev.off()
write.table(output$pairs, file="somaticInteractions.pairwise.tsv", quote=FALSE, row.names=FALSE, sep="\t")
write.table(output$gene_sets, file="somaticInteractions.comet.tsv", quote=FALSE, row.names=FALSE, sep="\t")
