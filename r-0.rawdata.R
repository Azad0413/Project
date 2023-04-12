rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/BJTC-399-11/")
if (! dir.exists("./00_rawdata")){
  dir.create("./00_rawdata")
}
setwd("./00_rawdata")
##01-1 TCGA -------
library(TCGAbiolinks)
library(readr)
library(readxl)
library(tidyverse)
## 读取从xena下载的数据
tcga.expr<-read_tsv(file = '/data/nas1/luchunlin/TCGA.matrix/TCGA-GBM.htseq_counts.tsv.gz')%>%as.data.frame()%>%
  column_to_rownames(var = 'Ensembl_ID')
## xena下载的数据经过了log2+1转化，需要将其还原
tcga.expr<-2^tcga.expr-1
## 对数据进行id转化
genecode<-read.table(file = '/data/nas1/luchunlin/pipeline/GENEANNO/gencode.v22.annotation.gene.probeMap')
probe2symbol<-genecode[,(1:2)]
colnames(probe2symbol)<-c('ID','symbol')
probe2symbol<-probe2symbol[-1,]
dat.tcga<-tcga.expr
dat.tcga$ID <- rownames(dat.tcga)
dat.tcga$ID<-as.character(dat.tcga$ID)
probe2symbol$ID<-as.character(probe2symbol$ID)
dat.tcga<-dat.tcga %>%
  inner_join(probe2symbol,by='ID')%>% 
  select(-ID)%>%     ## 去除多余信息
  select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('TCGA',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
dim(dat.tcga)

## 筛选癌症组织，去掉癌旁组织。01-09为肿瘤，10-19为正常对照
mete=data.frame(colnames(dat.tcga))  # 取第一行样本id
for (i in 1:length(mete[,1])) {
  num=as.numeric(as.character(substring(mete[i,1],14,15)))
  if(num %in% seq(1,9)){mete[i,2]="T"}
  if(num %in% seq(10,29)){mete[i,2]="N"}
}
names(mete)=c("id","group")
table(mete$group)
mete$group=as.factor(mete$group)
mete=subset(mete,mete$group=="T")
exp_tumor<-dat.tcga[,which(colnames(dat.tcga)%in%mete$id)]
exp_tumor<-as.data.frame(exp_tumor)
#   168
## 保留有生存数据的
survival<-read.delim2('/data/nas1/luchunlin/TCGA_survival/TCGA-GBM.survival.tsv')
exp_tumor<-exp_tumor[,colnames(exp_tumor)%in%survival$sample]
##167
exp_control<-dat.tcga[,which(!colnames(dat.tcga)%in%mete$id)]
exp_control<-as.data.frame(exp_control)
# 5
dat.final<-cbind(exp_control,exp_tumor)
##172
write.table(dat.final,file = 'dat.tcga.xls',sep = '\t',quote = F,row.names = T)
##fpkm
expr_fpkm<-read_tsv(file = '/data/nas1/luchunlin/TCGA.matrix/TCGA-GBM.htseq_fpkm.tsv.gz')%>%as.data.frame()%>%
  column_to_rownames(var = 'Ensembl_ID')
## xena下载的数据经过了log2+1转化，需要将其还原
expr_fpkm<-2^expr_fpkm-1
## 对数据进行id转化
dat_fpkm<-expr_fpkm
dat_fpkm$ID <- rownames(dat_fpkm)
dat_fpkm$ID<-as.character(dat_fpkm$ID)
probe2symbol$ID<-as.character(probe2symbol$ID)
dat_fpkm<-dat_fpkm %>%
  inner_join(probe2symbol,by='ID')%>% 
  select(-ID)%>%     ## 去除多余信息
  select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('TCGA',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
dim(dat_fpkm)
# dat_fpkm<-dat_fpkm[mRNA$gene_name,]
dat_fpkm<-dat_fpkm[,colnames(dat.final)]
dat_fpkm<-dat_fpkm[rownames(dat.final),]
#dat_fpkm<-na.omit(dat_fpkm)

write.table(dat_fpkm,file = 'dat.fpkm.xls',sep = '\t',row.names = T,quote = F)

# # fpkm转TPM
# FPKM2TPM <- function(fpkm){
#   exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
# }
# 
# dat_tpm <- apply(dat_fpkm,2,FPKM2TPM)
# write.table(dat_tpm,file = 'dat.tpm.xls',sep = '\t',row.names = T,quote = F)


# 02-GTEx数据的获取 -------------------------------------------------------------

GTEx <- read.table('GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz',
                   header = T,
                   sep = '\t',
                   skip = 2)

# 上传读取注释文件
GTEx_anno <- read.table('GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt',
                        header = T,sep = '\t',quote = '')
table(GTEx_anno$SMTS)
### 从注释文件GTEx_anno中提取的样本ID
GTEx_brain <- GTEx_anno[which(GTEx_anno$SMTS == "Brain"),] 
### 将GTEx_muscle文件中SAMPID列的'-'替换成'.'
GTEx_brain$SAMPID <- gsub('-','.',GTEx_brain$SAMPID,fixed = T)
head(GTEx)
### 提取矩阵
#GTEx_dat <- GTEx[colnames(GTEx)%in%GTEx_brain$SAMPID,]
cols <- intersect(GTEx_brain$SAMPID, colnames(GTEx))
GTEx_dat <- subset(GTEx, select = cols)
Name <- GTEx$Name
GTEx_dat <- cbind(Name, GTEx_dat)
### 调整格式:去除GTEx_dat文件Name列Ensembl号后面的点号
#GTEx_dat <- separate(GTEx_dat,
#                     Name,
#                     into=c('Name'),
#                     sep='\\.')
ucsc <- read.delim2('gtex_gene_expected_count.gz')
ucsc <- column_to_rownames(ucsc,var = 'sample')
GTEx_dat <- column_to_rownames(GTEx_dat,var = 'Name')
GTEx_dat <- GTEx_dat[,colnames(GTEx_dat)%in%colnames(ucsc)]

### 表达量不为0的矩阵(去除基因表达量为0的行),即:56200-2162=54038
#GTEx_dat <- GTEx_dat[which(rowSums(GTEx_dat) > 0),] 

genecode<-read.table(file = 'probeMap_gencode.v23.annotation.gene.probemap')
probe2symbol<-genecode[,(1:2)]
colnames(probe2symbol)<-c('ID','symbol')
probe2symbol<-probe2symbol[-1,]
dat.GTEx<-GTEx_dat
dat.GTEx$ID <- rownames(dat.GTEx)
dat.GTEx$ID<-as.character(dat.GTEx$ID)
probe2symbol$ID<-as.character(probe2symbol$ID)
dat.GTEx<-dat.GTEx %>%
  inner_join(probe2symbol,by='ID')%>% 
  select(-ID)%>%     ## 去除多余信息
  select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GTEX',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
dim(dat.GTEx)
# 51211  1034


# 03-3去除批次效应 --------------------------------------------------------------
library(sva)
# BiocManager::install("bladderbatch")
library(bladderbatch)
## 1.样本分组
group_GTEx <- data.frame(sample = colnames(dat.GTEx),
                         group = 'normal')
group_TCGA <- data.frame(sample = colnames(dat.final),
                         group = c(rep('normal',5),rep('tumor',167)))
group_list<- rbind(group_GTEx, group_TCGA)
write.table(group_list,file = 'group.merge.xls',sep = '\t',row.names = F,quote = F)
## 2.设置批次

batch <- group_list
batch$batch <- (c(rep('1',1034),rep('2',172)))
class(batch$batch)
batch$group <- factor(batch$group,levels = c('normal','tumor'))

library(lance)
dat.GTEx <- lc.tableToNum(dat.GTEx)
dat.GTEx1 <- dat.GTEx%>%rownames_to_column(var = 'symbol')
dat.final1 <- dat.final%>%rownames_to_column(var = 'symbol')
merge <- merge(dat.GTEx1,dat.final1,by='symbol')%>%column_to_rownames(var = 'symbol')
### 注:batch含样本名、分组和批次,共3列信息

## 3.去除批次效应
library(sva)
Combat_merge <- ComBat_seq(as.matrix(merge),
                           batch = batch$batch,group = batch$group)
Combat_merge <- as.data.frame(Combat_merge)
dat.all <- Combat_merge
write.table(dat.all,file = 'dat.merge.xls',sep = '\t',row.names = T,quote = F)
dat.all <- read.delim2('dat.merge.xls',row.names = 1)%>%lc.tableToNum()
pcg <- read.delim2('/data/nas1/luchunlin/pipeline/PCG/PCG.xls(v22)')
dat.pcg <- dat.all[pcg$gene_name,]
dat.pcg <- na.omit(dat.pcg)
write.table(dat.pcg,file = 'dat.pcg.xls',sep = '\t',row.names = T,quote = F)

keep<-rowSums(dat.all>0)>=floor(0.75*ncol(dat.all))

dat.fliter.all<-dat.all[keep,]
#row <- data.frame(rownames(dat.all))
write.table(dat.fliter.all,file = 'dat.merge.fliter.xls',sep = '\t',row.names = T,quote = F)
geneset <- read.delim2('../02_DEERGs/gene.txt')

unmap <- geneset[!geneset$id%in%rownames(dat.pcg),]
unmap
