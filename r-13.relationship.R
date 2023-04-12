rm(list = ls())
#10 诊断基因功能相似性分析（GOSemSim包）----------
setwd("/data/nas1/luchunlin/project/BJTC-302")
if (! dir.exists("./13_relationship")){
  dir.create("./13_relationship")
}
setwd("./13_relationship")
library(GOSemSim)
library(reshape2)
library(org.Hs.eg.db)
bp<-godata('org.Hs.eg.db',ont = "BP",computeIC = F)
cc<-godata('org.Hs.eg.db',ont = "CC",computeIC = F)
mf<-godata('org.Hs.eg.db',ont = "MF",computeIC = F)
hubgene <- read.delim2('/data/nas1/luchunlin/project/BJTC-302/10_validation/hubgene.final.xls')
geneid2<-hubgene$hubgene
gene_transform2<-bitr(geneid2,fromType = "SYMBOL",
                      toType = "ENTREZID",
                      OrgDb = org.Hs.eg.db)
genelist<-gene_transform2
simbp<-mgeneSim(genelist$ENTREZID,
                semData = bp,
                measure = "Wang",
                drop = NULL,
                combine = "BMA")
simcc<-mgeneSim(genelist$ENTREZID,
                semData = cc,
                measure = "Wang",
                drop = NULL,
                combine = "BMA")
simmf<-mgeneSim(genelist$ENTREZID,
                semData = mf,
                measure = "Wang",
                drop = NULL,
                combine = "BMA")

## 基于相似度结果，进一步计算基因在BP、CC、MF层面的集合平均值，得到最终评分。
fsim<-(simmf*simcc*simbp)^(1/3)
View(fsim)
colnames(fsim)<-genelist$SYMBOL
rownames(fsim)<-genelist$SYMBOL

## 进一步去除基因与本身基因之间的相关性，使用melt函数将宽格式数据转化为长格式数据
for (i in 1:ncol(fsim)) {
  fsim[i,i] <- NA
}
dat_sim<-melt(fsim)
View(dat_sim)
dat_sim<-dat_sim[!is.na(dat_sim$value),]
dat_sim<-dat_sim[,c(1,3)]
head(dat_sim)
## 绘制boxplot图
dat.median<-aggregate(value~Var1,dat_sim,median)
dat.mean<-aggregate(value~Var1,dat_sim,mean)
View(dat.mean)
## 根据相似性评分的平均值或中位值（根据自身需要），对其进行排序，并根据评分的高低，将基因名设置为因子（factor）格式
#m<-dat.mean$value
m<-dat.median$value
#names(m)<-dat.mean$Var1
names(m)<-dat.median$Var1
dat_sim$Var1<-factor(dat_sim$Var1,
                     levels = names(sort(m)))
str(dat_sim)
write.table(fsim,file = 'relationship.xls',
            sep = '\t',
            row.names = T)
## 使用ggplot 对分析结果进行可视化。
ggplot(dat_sim,
       aes(x=Var1,y=value,fill=factor(Var1)))+
  scale_fill_brewer(palette = "Set3")+
  geom_boxplot()+
  coord_flip()+       # 坐标轴互换
  xlab("")+ylab("")+
  theme_bw()+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 18),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13))+
  guides(fill='none')
write.table(fsim,file = 'relationship.xls',
            sep = '\t',
            quote = F)
ggsave('relationship.pdf',width = 7,height = 5)
ggsave('relationship.png',width = 7,height = 5)
