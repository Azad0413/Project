rm(list = ls())
setwd("/data/nas1/luchunlin/project/JNZK-255-2/")
if (! dir.exists("./14_TMB")){
  dir.create("./14_TMB")
}
setwd("./14_TMB")
## 高低风险组患者TMB评分箱线图
risk <- read.delim2('../06_risk/risk.xls')
group <- select(risk,c('id','risk','riskScore'))
# group$riskScore <- as.numeric(group$riskScore)

group$group <- ifelse(group$risk==0,'High risk','Low risk')
colnames(group) <- c('sample','risk','riskScore','group')
library(TCGAmutations)
tmp<-as.data.frame(tcga_available())
dt<-TCGAmutations::tcga_load(study = 'BLCA')
dt<-dt@data
dt1<-as.data.frame(table(dt$Tumor_Sample_Barcode))
names(dt1)<-c('Barcode','Freq')
dt1$tmb<-dt1$Freq/38
tmb<-dt1
tmb$Barcode <- substr(tmb$Barcode,1,16)
tmb <- tmb[tmb$Barcode%in%risk$id,]
colnames(tmb) <- c('sample','Freq','tmb')
tmb<-merge(tmb,group,by='sample')
# tmb<-tmb[-c(80,97),]
## 38是人类基因外显子的长度
my_comparisons <- list(c("High risk","Low risk"))
tmb.bar<-ggplot(tmb,aes(x = group, y = tmb, fill = group)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "TMB",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("TMB score") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              map_signif_level = T,
              y_position = c(100))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 20),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 20))+
  guides(fill='none')
tmb.bar
ggsave(tmb.bar,filename = '01.tmb.pdf',w=4,h=5)
ggsave(tmb.bar,filename = '01.tmb.png',w=4,h=5)

##风险评分与TMB相关性------
cor_dat <- tmb

cor_dat$riskScore <- as.numeric(cor_dat$riskScore)

p1<-ggplot(cor_dat,aes(x=riskScore,y=tmb))+geom_point(aes(x=riskScore,y=tmb,color =group), size=2,alpha=0.9)+geom_smooth(method = 'lm', formula = y ~ x, se = T,color='black')+ stat_cor(data=cor_dat, method = "spearman")+
  scale_fill_manual(values =c("black")) +
  theme_bw()+
  scale_color_manual(values = c("#DB3587","#35B3B2"), name = "Risk") +
  theme(axis.title = element_text(size = 20, face = "bold", family = "Times",color='black'),
        axis.text.x = element_text(size = 16, face = "bold", family = "Times",color='black'),
        axis.text.y = element_text(size = 16, face = "bold", family = "Times",color='black'),
        legend.text = element_text(size = 14, face = "bold", family = "Times",color='black'),
        legend.title= element_text(size =16, face = "bold", family = "Times",color='black'),
        legend.position = "bottom"
  )+
  labs(x="Risk score",y="Tumor Burden Mutation")+
  theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
p1

pdf('02.correlation.pdf',width=7,height=6,family='Times')
print(p1)
dev.off()

png('02.correlation.png',width=7,height=6,family='Times',units='in',res=600)
print(p1)
dev.off()

### km曲线
km.dat <- cor_dat
km.dat$tmb.group <- ifelse(km.dat$tmb>median(km.dat$tmb),'High_TMB','Low_TMB')
table(km.dat$tmb.group)
km.dat$group.all <- paste0(km.dat$tmb.group,'+',km.dat$group)
table(km.dat$group.all)

survival <- read.delim2('../04_univariate_cox/survival.xls')
km.dat <- merge(survival,km.dat,by='sample')

kmfit <- survfit(Surv(OS.time,OS) ~ group.all, data = km.dat)
table(km.dat$group.all)
# 绘制KM曲线
KM <- ggsurvplot(kmfit,
                 pval = TRUE, 
                 conf.int = F,
                 legend.labs=c("H-TMB+High risk","H-TMB+Low risk","L-TMB+High risk","Low-TMB+Low risk"),
                 legend.title="group",
                 title="TMB+riskScore",
                 font.main = c(15,"bold"),
                 risk.table = TRUE, 
                 risk.table.col = "strata",  
                 linetype = "strata", 
                 surv.median.line = "hv", 
                 ggtheme = theme_bw(), 
                 palette = c("#A73030FF","#0073C2FF","darkgreen","purple"))
KM

## PART A TMB------
library(TCGAmutations)
maf <- TCGAmutations::tcga_load(study = "BLCA")
library(maftools)
maf_sample <- data.frame(barcode = maf@clinical.data$Tumor_Sample_Barcode)
maf_sample$sample <- stringr::str_sub(maf_sample$barcode, 1, 16)
sample<-subset(maf_sample)$barcode
maf<-subsetMaf(maf,tsb = sample)
model_gene <- read.delim2('../05_Lasso/lasso_genes.csv',header = F)

maf<-subsetMaf(maf,tsb = sample,genes = model_gene$V1)
sample<-subset(maf_sample)$barcode
maf<-subsetMaf(maf,tsb = sample)
#maf<-subsetMaf(maf,tsb = sample,genes = model_gene$V1)
High.sample <- group$sample[which(group$group=='High risk')]
Low.sample <- group$sample[which(group$group=='Low risk')]
##HIGH LOW分开
high <- maf_sample[(maf_sample$sample)%in%High.sample,]
high <- subset(high)$barcode
maf.high <- subsetMaf(maf,tsb = high)
# pdf(file = "01.oncoplot.high.pdf", height = 8, width = 10)
# oncoplot(maf = maf.high, top = 20)
# dev.off()
# 
# png(file = "01.oncoplot.high.png", family = "Times", height = 8, width = 10, units = "in", res = 600)
# oncoplot(maf = maf.high, top = 20)
# dev.off()
# 
low <- maf_sample[(maf_sample$sample)%in%Low.sample,]
low <- subset(low)$barcode
maf.low <- subsetMaf(maf,tsb = low)
# pdf(file = "02.oncoplot.low.pdf", height = 8, width = 10)
# oncoplot(maf = maf.low, top = 20)
# 
# dev.off()

# png(file = "02.oncoplot.low.png", family = "Times", height = 8, width = 10, units = "in", res = 600)
# oncoplot(maf = maf.low, top = 20)
# dev.off()

##预后基因突变率--------

pdf(file = "04.mut_high.pdf", family = "Times", height = 5, width = 6)
plotmafSummary(maf = maf.high, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

png(file = "04.mut_high.png", family = "Times", height = 5, width = 6, units = "in", res = 600)
plotmafSummary(maf = maf.high, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

pdf(file = "05.mut_low.pdf", family = "Times", height = 5, width = 6)
plotmafSummary(maf = maf.low, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

png(file = "05.mut_low.png", family = "Times", height = 5, width = 6, units = "in", res = 600)
plotmafSummary(maf = maf.low, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

### IDO1 IDO2 TDO1---------
intreast <- c('IDO1','IDO2','TDO2')

library(TCGAmutations)
maf <- TCGAmutations::tcga_load(study = "BLCA")
library(maftools)
maf_sample <- data.frame(barcode = maf@clinical.data$Tumor_Sample_Barcode)
maf_sample$sample <- stringr::str_sub(maf_sample$barcode, 1, 16)
sample<-subset(maf_sample)$barcode
maf<-subsetMaf(maf,tsb = sample)

maf<-subsetMaf(maf,tsb = sample,genes = intreast)
sample<-subset(maf_sample)$barcode
maf<-subsetMaf(maf,tsb = sample)
#maf<-subsetMaf(maf,tsb = sample,genes = model_gene$V1)
High.sample <- group$sample[which(group$group=='High risk')]
Low.sample <- group$sample[which(group$group=='Low risk')]
##HIGH LOW分开
high <- maf_sample[(maf_sample$sample)%in%High.sample,]
high <- subset(high)$barcode
maf.high <- subsetMaf(maf,tsb = high)
pdf(file = "06.oncoplot.high.pdf", height = 4, width = 6)
oncoplot(maf = maf.high, top = 20)
dev.off()
# 
png(file = "06.oncoplot.high.png", family = "Times", height = 4, width = 6, units = "in", res = 600)
oncoplot(maf = maf.high, top = 20)
dev.off()
# 
low <- maf_sample[(maf_sample$sample)%in%Low.sample,]
low <- subset(low)$barcode
maf.low <- subsetMaf(maf,tsb = low)
pdf(file = "07.oncoplot.low.pdf", height = 4, width = 6)
oncoplot(maf = maf.low, top = 20)
# 
dev.off()

png(file = "07.oncoplot.low.png", family = "Times", height = 4, width = 6, units = "in", res = 600)
oncoplot(maf = maf.low, top = 20)
dev.off()

