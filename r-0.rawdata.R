rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/BJTC-299_modIfied/")
if (! dir.exists("./00_rawdata")){
  dir.create("./00_rawdata")
}
setwd("./00_rawdata")

library(readr)
library(tidyverse)
library(org.Hs.eg.db)
dat.fpkm <- read_tsv('fpkm.tsv')%>%column_to_rownames(var = 'gene')
dat.count <- read_tsv('count.tsv')%>%column_to_rownames(var = 'gene')


grps = ifelse(colnames(dat.count) %>% str_ends(".1.."),"normal","tumor")
table(grps)
df.meta = data.frame(sample = colnames(dat.count), group = grps)
df.meta <- df.meta[order(df.meta$group),]
dat.count <- dat.count[,df.meta$sample]
dat.fpkm <- dat.fpkm[,df.meta$sample]
gene.transform = AnnotationDbi::select(org.Hs.eg.db, rownames(dat.count), "SYMBOL", "ENTREZID")
class(gene.transform$SYMBOL)
dat.count <- dat.count[!duplicated(dat.count$symbol),]%>%na.omit()
dat.fpkm <- dat.fpkm[rownames(dat.count),]
rownames(dat.fpkm) <- dat.count$symbol
rownames(dat.count) <- dat.count$symbol
dat.count <- dplyr::select(dat.count,-'symbol')
write.table(dat.count,file = 'dat.count.xls',sep = '\t',row.names = T,quote = F)
write.table(dat.fpkm,file = 'dat.fpkm.xls',sep = '\t',row.names = T,quote = F)

# 02-GTEx数据的获取 -------------------------------------------------------------

GTEx <- read.table('/data/nas1/luchunlin/project/BJTC-399-11/00_rawdata/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz',
                   header = T,
                   sep = '\t',
                   skip = 2)

# 上传读取注释文件
GTEx_anno <- read.table('/data/nas1/luchunlin/project/BJTC-399-11/00_rawdata/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt',
                        header = T,sep = '\t',quote = '')
table(GTEx_anno$SMTS)
### 从注释文件GTEx_anno中提取的样本ID
GTEx_Bladder <- GTEx_anno[which(GTEx_anno$SMTS == "Bladder"),] 
### 将GTEx_muscle文件中SAMPID列的'-'替换成'.'
GTEx_Bladder$SAMPID <- gsub('-','.',GTEx_Bladder$SAMPID,fixed = T)
head(GTEx)
### 提取矩阵
#GTEx_dat <- GTEx[colnames(GTEx)%in%GTEx_Bladder$SAMPID,]
cols <- intersect(GTEx_Bladder$SAMPID, colnames(GTEx))
GTEx_dat <- subset(GTEx, select = cols)
Name <- GTEx$Name
GTEx_dat <- cbind(Name, GTEx_dat)
### 调整格式:去除GTEx_dat文件Name列Ensembl号后面的点号
# GTEx_dat <- separate(GTEx_dat,
#                     Name,
#                     into=c('Name'),
#                     sep='\\.')
# ucsc <- read.delim2('gtex_RSEM_gene_fpkm.gz')
# ucsc <- column_to_rownames(ucsc,var = 'sample')

# GTEx_dat <- column_to_rownames(GTEx_dat,var = 'Name')
# GTEx_dat <- GTEx_dat[,colnames(GTEx_dat)%in%colnames(ucsc)]


### 表达量不为0的矩阵(去除基因表达量为0的行),即:56200-2162=54038
#GTEx_dat <- GTEx_dat[which(rowSums(GTEx_dat) > 0),] 

genecode<-read.table(file = '/data/nas1/luchunlin/project/BJTC-399-11/00_rawdata/probeMap_gencode.v23.annotation.gene.probemap')
probe2symbol<-genecode[,(1:2)]
colnames(probe2symbol)<-c('ID','symbol')
probe2symbol<-probe2symbol[-1,]
dat.GTEx<-GTEx_dat
colnames(dat.GTEx)[1] <- 'ID'
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
# 51211  21
unmap <- GTEx_dat[GTEx_dat$Name%in%'ENSG00000115380.19',]
rownames(unmap) <- 'EFEMP1'
unmap <- unmap[,-1]
dat.GTEx <- rbind(unmap,dat.GTEx)

# 03-3去除批次效应 --------------------------------------------------------------
library(sva)
# BiocManager::install("bladderbatch")
library(bladderbatch)
## 1.样本分组
group_GTEx <- data.frame(sample = colnames(dat.GTEx),
                         group = 'normal')
table(df.meta$group)
group_TCGA <- data.frame(sample = colnames(dat.count),
                         group = c(rep('normal',19),rep('tumor',411)))
group_list<- rbind(group_GTEx, group_TCGA)
write.table(group_list,file = 'group.merge.xls',sep = '\t',row.names = F,quote = F)
## 2.设置批次
table(group_list$group)
batch <- group_list
batch$batch <- (c(rep('1',21),rep('2',430)))
class(batch$batch)
batch$group <- factor(batch$group,levels = c('normal','tumor'))

library(lance)
dat.GTEx1 <- dat.GTEx%>%rownames_to_column(var = 'symbol')
dat.final1 <- dat.count%>%rownames_to_column(var = 'symbol')
merge <- merge(dat.GTEx1,dat.final1,by='symbol')%>%column_to_rownames(var = 'symbol')
### 注:batch含样本名、分组和批次,共3列信息
# BiocManager::install('sva')

## 3.去除批次效应
library(sva)
Combat_merge <- ComBat_seq(as.matrix(merge),
                           batch = batch$batch,group = batch$group)
Combat_merge <- as.data.frame(Combat_merge)
dat.all <- Combat_merge
write.table(dat.all,file = 'dat.merge.xls',sep = '\t',row.names = T,quote = F)

###配对样本
group.match <- df.meta
group.match$match <- substr(group.match$sample,1,12)
tumor <- group.match[which(group.match$group=='tumor'),]
normal <- group.match[which(group.match$group=='normal'),]
tumor.match <- tumor[tumor$match%in%normal$match,]
tumor.match <- tumor.match[!duplicated(tumor.match$match),]
group.match <- rbind(normal,tumor.match)%>%dplyr::select(-'match')
write.table(group.match,file = 'group.match.xls',sep = '\t',row.names = F,quote = F)
