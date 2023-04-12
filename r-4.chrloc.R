## 02 差异基因鉴定-----
rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-321/")
if (! dir.exists("./04_chrloc")){
  dir.create("./04_chrloc")
}
setwd("./04_chrloc")
library(OmicCircos)
library(rtracklayer)
library(RCircos)
dat<-read.delim2('/data/nas1/luchunlin/project/BJTC-321/00_rawdata/dat.final.xls',row.names = 1)%>%lc.tableToNum()
group<-read.delim2("/data/nas1/luchunlin/project/BJTC-321/00_rawdata/group.xls")
table(group$group)
control.sample<-group$sample[which(group$group=='control')]
asthma.sample<-group$sample[which(group$group=='asthma')]
DEAEG<-read.delim2('/data/nas1/luchunlin/project/BJTC-321/03_DEAEG/DEAEG.xls')
gencode.gtf<-import.gff(con = '/data/nas1/luchunlin/pipeline/gtf/gencode.v19.annotation.gtf.gz')
chr<-data.frame(gencode.gtf)
chr<-subset(chr,source=='HAVANA'&type=='gene')  ##HAVANA
chr<-chr[chr$gene_name%in%DEAEG$.,]
unmap<-DEAEG[!DEAEG$.%in%chr$gene_name,]
data(RCircos.Gene.Label.Data)
Circos.Gene.Label<-data.frame(chr$seqnames,chr$start,chr$end,chr$gene_name)  
colnames(Circos.Gene.Label)<-colnames(RCircos.Gene.Label.Data)
Circos.Gene.Label$pos<-NA

cnt<-1
while (cnt < 62) {
  Circos.Gene.Label$pos[cnt]<-mean(c(Circos.Gene.Label$chr.start[cnt],Circos.Gene.Label$chromEnd[cnt]))
  cnt = cnt + 1
}

location<-Circos.Gene.Label[,c('Chromosome','pos','Gene')]
unmap<-data.frame(Chromosome=c('chrX'),pos=c('135551585'),Gene=c('INTS6L'))
location<-rbind(location,unmap)
exp_chromosome<-dat[rownames(dat)%in%location$Gene,]
## 构建热图数据(分成两组画一下)
exp.asthma<-exp_chromosome[,asthma.sample]
exp.control<-exp_chromosome[,control.sample]
## 实验组样本
chr_heatmap1<-cbind(location,exp.asthma)
## 对照组样本
chr_heatmap2<-cbind(location,exp.control)
## gene fusion列表
gene_fus<-cbind(location,location)
colnames(gene_fus)<-c('chr1','po1','gene1','chr2','po2','gene2')
gene_fus<-gene_fus[,c(1:3,6)]

par(mar=c(2,2,2,2))
plot(c(1,800),c(1,800),type='n',axes=F,xlab='',ylab='',main='')
circos(R=300,cir = 'hg18',W=15,type = 'chr',print.chr.lab = T,scale = T)
circos(R=330,cir = 'hg18',W=10,mapping=gene_fus,type = 'label',side = 'out',col=c('black','blue','red'),cex = 0.4)
circos(R=190,cir = 'hg18',W=100,mapping = chr_heatmap1,col.v = 4,type = 'heatmap2',cluster = F,col.bar = F,lwd = 0.5,col = 'blue')
circos(R=100,cir = 'hg18',W=100,mapping = chr_heatmap2,col.v = 4,type = 'heatmap2',cluster = F,col.bar = F,lwd = 0,col = 'blue')
