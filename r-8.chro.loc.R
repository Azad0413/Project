rm(list = ls())
# 染色体定位----------
setwd("/data/nas1/luchunlin/project/BJTC-317")
if (! dir.exists("./08_chromosome_loc")){
  dir.create("./08_chromosome_loc")
}
setwd("./08_chromosome_loc")
DEMID<-read.delim2('/data/nas1/luchunlin/project/BJTC-317/03_DEMTG/DEMID.xls')
dat<-read.delim2("/data/nas1/luchunlin/project/BJTC-317/00_rawdata/dat.final.xls", row.names = 1)  %>% lc.tableToNum
group = read.delim2("/data/nas1/luchunlin/project/BJTC-317/00_rawdata/group.xls")
## 'OmicCircos'包 
# BiocManager::install('OmicCircos')
library(OmicCircos)
## 想要的结果：从外到内依次是基因名(SNP位点)、染色体定位、热图。
## 需要的数据：（1）一个gene list（chr po gene）(2)染色体图用自带的 （3）热图（chr po gene exp）
## 先将11个基因的表达矩阵提取出来
exp_chromosome<-dat[rownames(dat)%in%DEMID$.,]%>%log2()
## 查找并整理染色体位置信息
## 节段数据（segment data）chrom start ened 
## 映射数据（mapping data）前面两列固定为节段名和位置 后面是表达量和甲基化
location<-read_xlsx('/data/nas1/luchunlin/project/BJTC-317/08_chromosome_loc/location.xlsx')

## 构建热图数据(分成两组画一下)
chr_heatmap<-cbind(location,exp_chromosome)
## 实验组样本
chr_heatmap1<-chr_heatmap[,c(1:3,12:23)]
## 对照组样本
chr_heatmap2<-chr_heatmap[,c(1:11)]
## gene fusion列表
gene_fus<-cbind(location,location)
colnames(gene_fus)<-c('chr1','po1','gene1','chr2','po2','gene2')
gene_fus<-gene_fus[,c(1:3,6)]
par(mar=c(2,2,2,2))
plot(c(1,800),c(1,800),type='n',axes=F,xlab='',ylab='',main='')
circos(R=300,cir = 'hg18',W=15,type = 'chr',print.chr.lab = T,scale = T)
circos(R=330,cir = 'hg18',W=10,mapping=gene_fus,type = 'label',side = 'out',col=c('black','blue','red'),cex = 0.5)
circos(R=190,cir = 'hg18',W=100,mapping = chr_heatmap1,col.v = 4,type = 'heatmap2',cluster = F,col.bar = F,lwd = 0.5,col = 'blue')
circos(R=100,cir = 'hg18',W=100,mapping = chr_heatmap2,col.v = 4,type = 'heatmap2',cluster = F,col.bar = F,lwd = 0,col = 'blue')
