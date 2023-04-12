rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-441-3/")
if (! dir.exists("./08_TMB")){
  dir.create("./08_TMB")
}
setwd("./08_TMB")

## 高低表达组患者TMB评分箱线图
group <- read.delim2('../05_survival/group(UGCG).xls')

library(TCGAmutations)
tmp<-as.data.frame(tcga_available())
dt<-TCGAmutations::tcga_load(study = 'LAML')
dt<-dt@data
dt1<-as.data.frame(table(dt$Tumor_Sample_Barcode))
names(dt1)<-c('Barcode','Freq')
dt1$tmb<-dt1$Freq/38
tmb<-dt1
tmb$Barcode <- substr(tmb$Barcode,1,16)
tmb <- tmb[tmb$Barcode%in%group$sample,]
colnames(tmb) <- c('sample','Freq','tmb')
tmb<-merge(tmb,group,by='sample')



## 38是人类基因外显子的长度
class(tmb$group)
tmb$group <- factor(tmb$group)
class(tmb$tmb)
my_comparisons <- list(c("High UGCG","Low UGCG"))
tmb.bar<-ggplot(tmb,aes(x = group, y = tmb, fill = group)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "TMB",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("TMB score") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              map_signif_level = F,
              y_position = c(1.3),
              test = 'wilcox.test')+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 20))+
  guides(fill='none')
tmb.bar
# ggsave(tmb.bar,filename = '01.tmb.pdf',w=4,h=5)
# ggsave(tmb.bar,filename = '01.tmb.png',w=4,h=5)

## PART A TMB------

maf <- TCGAmutations::tcga_load(study = "LAML")
library(maftools)
maf_sample <- data.frame(barcode = maf@clinical.data$Tumor_Sample_Barcode)
maf_sample$sample <- stringr::str_sub(maf_sample$barcode, 1, 16)
maf_sample <- merge(maf_sample, group, by.x = "sample",by.y = 'sample')

sample<-subset(maf_sample)$barcode
maf<-subsetMaf(maf,tsb = sample)
#maf<-subsetMaf(maf,tsb = sample,genes = model_gene$V1)
High.sample <- group$sample[which(group$group=='High UGCG')]
##HIGH LOW分开
high <- maf_sample[(maf_sample$sample)%in%High.sample,]
high <- subset(high)$barcode
maf.high <- subsetMaf(maf,tsb = high)
pdf(file = "01.oncoplot.high.pdf", height = 8, width = 10)
oncoplot(maf = maf.high, top = 20)
dev.off()

png(file = "01.oncoplot.high.png", height = 8, width = 10, units = "in", res = 600)
oncoplot(maf = maf.high, top = 20)
dev.off()

Low.sample <- group$sample[which(group$group=='Low UGCG')]

low <- maf_sample[(maf_sample$sample)%in%Low.sample,]
low <- subset(low)$barcode
maf.low <- subsetMaf(maf,tsb = low)
pdf(file = "02.oncoplot.low.pdf", height = 8, width = 10)
oncoplot(maf = maf.low, top = 20)
dev.off()

png(file = "02.oncoplot.low.png",height = 8, width = 10, units = "in", res = 600)
oncoplot(maf = maf.low, top = 20)
dev.off()

# ##预后基因突变率--------
# maf.mod<-subsetMaf(maf,tsb = sample,genes = 'UGCG')
# 
# pdf(file = "mut_modgene.pdf", family = "Times", height = 5, width = 6)
# plotmafSummary(maf = maf.mod, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
# dev.off()
# 
# png(file = "mut_modgene.png", family = "Times", height = 5, width = 6, units = "in", res = 600)
# plotmafSummary(maf = maf.mod, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
# dev.off()