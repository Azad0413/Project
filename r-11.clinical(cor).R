rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-386-10/")
if (! dir.exists("./11_clinical(cor)")){
  dir.create("./11_clinical(cor)")
}
setwd("./11_clinical(cor)")

##lncRNA与临床特征相关性

## 训练集-------
dat<-read.delim2('../00_rawdata/dat.fpkm.xls',row.names = 1)%>%lc.tableToNum()
colnames(dat)<-gsub('.','-',colnames(dat),fixed = T)
dat <- dat[,-c(1:72)]
hubgene <- read.delim2('../07_hubgene/hubgene.xls')
dat <- t(dat[hubgene$symbol,])%>%as.data.frame()
phenotype<-read.delim2('/data/nas1/luchunlin/TCGA_phenotype/TCGA-KIRC.GDC_phenotype.tsv.gz')
train_phenotype<-data.frame(sample=phenotype$submitter_id.samples,
                            age=phenotype$age_at_initial_pathologic_diagnosis,
                            gender=phenotype$gender.demographic,
                            grade=phenotype$neoplasm_histologic_grade,
                            stage=phenotype$tumor_stage.diagnoses,
                            T.stage=phenotype$pathologic_T,
                            N.stage=phenotype$pathologic_N,
                            M.stage=phenotype$pathologic_M)
dat$sample <- rownames(dat)
train_phenotype<-merge(train_phenotype,dat,by='sample')
#write.table(train_phenotype,file = 'phenotype.all.xls',sep = '\t',row.names = F,quote = F)

train_phenotype2<-train_phenotype
table(train_phenotype2$stage)
train_phenotype2$stage<-gsub('not reported',NA,train_phenotype2$stage)

train_phenotype2$stage<-gsub('stage iv','4',train_phenotype2$stage)
train_phenotype2$stage<-gsub('stage iii','3',train_phenotype2$stage)
train_phenotype2$stage<-gsub('stage ii','2',train_phenotype2$stage)
train_phenotype2$stage<-gsub('stage i','1',train_phenotype2$stage)

table(train_phenotype2$age)
#train_phenotype2$age<-ifelse(train_phenotype2$age>60,'1','0')

table(train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('a','',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('b','',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('c','',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('T4','3',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('T3','3',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('T2','2',train_phenotype2$T.stage)
train_phenotype2$T.stage<-gsub('T1','1',train_phenotype2$T.stage)

table(train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('NX',NA,train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('N0','0',train_phenotype2$N.stage)
train_phenotype2$N.stage<-gsub('N1','1',train_phenotype2$N.stage)

table(train_phenotype2$M.stage)
train_phenotype2$M.stage<-gsub('MX',NA,train_phenotype2$M.stage)
train_phenotype2$M.stage<-gsub('M0','0',train_phenotype2$M.stage)
train_phenotype2$M.stage<-gsub('M1','1',train_phenotype2$M.stage)
train_phenotype2$M.stage[train_phenotype2$M.stage==''] <- NA

table(train_phenotype2$grade)
train_phenotype2$grade<-gsub('GX',NA,train_phenotype2$grade)
train_phenotype2$grade<-gsub('G4','4',train_phenotype2$grade)
train_phenotype2$grade<-gsub('G3','3',train_phenotype2$grade)
train_phenotype2$grade<-gsub('G2','2',train_phenotype2$grade)
train_phenotype2$grade<-gsub('G1','1',train_phenotype2$grade)
train_phenotype2$grade[train_phenotype2$grade==''] <- NA

table(train_phenotype2$gender)
train_phenotype2$gender <- ifelse(train_phenotype2$gender=='male',1,0)
colnames(train_phenotype2)
colnames(train_phenotype2)<-c('id','Age','Gender','Grade','Stage','T.stage','N.stage','M.stage','CTD-2626G11.2','AP000696.2','RP11-528A4.2','LINC00645','RP4-655J12.4','RP11-321G12.1','RP11-195B3.1','PTCSC3')

## 相关性分析
library(Hmisc)
nc <- train_phenotype2%>%column_to_rownames(var = 'id')%>%as.data.frame()%>%lc.tableToNum()
nc <- na.omit(nc)
nc <- as.matrix(nc)
exp <- t(nc[,c(8:15)])
clinical <- t(nc[,c(1:7)])%>%lc.tableToNum()

m=rcorr(nc)$r[1:nrow(clinical),(ncol(nc)-length(hubgene$symbol)+1):ncol(nc)]
m<-t(m)
p=rcorr(nc)$P[1:nrow(clinical),(ncol(nc)-length(hubgene$symbol)+1):ncol(nc)]
p<-t(p)
library(dplyr)

# tmp = matrix(case_when(p<0.0001~"****",
#                        p<0.0001~"***",
#                        p<0.01~"**",
#                        p<0.05~"*",
#                        T~""),nrow = nrow(p))
cor <- m
cor <- signif(cor,3)
#cor[abs(cor)<0.21] <- ''
p <- signif(p,digits = 3)
textMatrix = paste(cor,"\n(",
                   p, ")",sep = "")
dim(textMatrix) = dim(m)
colnames(textMatrix)<-colnames(m)
rownames(textMatrix)<-rownames(m)
pdf(file = '01.lncRNA-Clinical_relationships(Train).pdf',w=7,h=6)
par(mar = c(9, 8, 3, 1));
labeledHeatmap(Matrix = m, 
               xLabels = colnames(m), 
               yLabels = rownames(m), 
               cex.lab = 1, 
               # ySymbols = colnames(MEs_col)[-20], 
               colorLabels = FALSE, 
               colors = blueWhiteRed(50, 0.8),
               textMatrix = textMatrix, 
               setStdMargins = FALSE, 
               cex.text = 0.8, 
               zlim = c(-1,1),
               main = paste("lncRNA-Clinical relationships"))
dev.off()

png(file = '01.lncRNA-Clinical_relationships(Train).png',w=600,h=500)
par(mar = c(9, 8, 3, 1));
labeledHeatmap(Matrix = m, 
               xLabels = colnames(m), 
               yLabels = rownames(m), 
               cex.lab = 1, 
               # ySymbols = colnames(MEs_col)[-20], 
               colorLabels = FALSE, 
               colors = blueWhiteRed(50, 0.8),
               textMatrix = textMatrix, 
               setStdMargins = FALSE, 
               cex.text = 0.8, 
               zlim = c(-1,1),
               main = paste("lncRNA-Clinical relationships"))
dev.off()

## validation---------
dat<-read.delim2('../00_rawdata/dat.va.xls')%>%lc.tableToNum()
dat <- na.omit(dat)
group <- read.delim2('../00_rawdata/group.va.xls')
tumor.sample <- group$sample[which(group$group=='Tumor')]
dat <- t(dat[,tumor.sample])%>%as.data.frame()
phenotype<-read.delim2('../00_rawdata/Merge_clinical.txt')
validation_phenotype<-data.frame(sample=phenotype$icgc_sample_id,
                            age=phenotype$donor_age_at_diagnosis,
                            gender=phenotype$donor_sex,
                            T.stage=phenotype$donor_tumour_stage_at_diagnosis,
                            N.stage=phenotype$donor_tumour_stage_at_diagnosis,
                            M.stage=phenotype$donor_tumour_stage_at_diagnosis)
dat$sample <- rownames(dat)
validation_phenotype<-merge(validation_phenotype,dat,by='sample')
write.table(validation_phenotype,file = 'phenotype.validation.xls',sep = '\t',row.names = F,quote = F)

validation_phenotype2<-validation_phenotype

table(validation_phenotype2$age)
#validation_phenotype2$age<-ifelse(validation_phenotype2$age>60,'1','0')

table(validation_phenotype2$T.stage)
validation_phenotype2$T.stage <- substr(validation_phenotype2$T.stage,1,2)

validation_phenotype2$T.stage<-gsub('T4','3',validation_phenotype2$T.stage)
validation_phenotype2$T.stage<-gsub('T3','3',validation_phenotype2$T.stage)
validation_phenotype2$T.stage<-gsub('T2','2',validation_phenotype2$T.stage)
validation_phenotype2$T.stage<-gsub('T1','1',validation_phenotype2$T.stage)

table(validation_phenotype2$N.stage)
validation_phenotype2$N.stage <- substr(validation_phenotype2$N.stage,3,4)
validation_phenotype2$N.stage<-gsub('NX',NA,validation_phenotype2$N.stage)
validation_phenotype2$N.stage<-gsub('N0','0',validation_phenotype2$N.stage)
validation_phenotype2$N.stage<-gsub('N1','1',validation_phenotype2$N.stage)

table(validation_phenotype2$M.stage)
validation_phenotype2$M.stage <- substr(validation_phenotype2$M.stage,5,6)
validation_phenotype2$M.stage<-gsub('MX',NA,validation_phenotype2$M.stage)
validation_phenotype2$M.stage<-gsub('M0','0',validation_phenotype2$M.stage)
validation_phenotype2$M.stage<-gsub('M1','1',validation_phenotype2$M.stage)

table(validation_phenotype2$gender)
validation_phenotype2$gender <- ifelse(validation_phenotype2$gender=='male',1,0)
colnames(validation_phenotype2)
colnames(validation_phenotype2)<-c('id','Age','Gender','T.stage','N.stage','M.stage','RP4-655J12.4','RP11-528A4.2','RP11-195B3.1','LINC00645','PTCSC3','RP11-321G12.1','AP000696.2')

## 相关性分析
library(Hmisc)
nc <- validation_phenotype2%>%column_to_rownames(var = 'id')%>%as.data.frame()%>%lc.tableToNum()
nc <- na.omit(nc)
nc <- as.matrix(nc)
exp <- t(nc[,c(6:12)])
clinical <- t(nc[,c(1:5)])%>%lc.tableToNum()

m=rcorr(nc)$r[1:nrow(clinical),(ncol(nc)-length(hubgene$symbol)+2):ncol(nc)]
m<-t(m)
p=rcorr(nc)$P[1:nrow(clinical),(ncol(nc)-length(hubgene$symbol)+2):ncol(nc)]
p<-t(p)
library(dplyr)

# tmp = matrix(case_when(p<0.0001~"****",
#                        p<0.0001~"***",
#                        p<0.01~"**",
#                        p<0.05~"*",
#                        T~""),nrow = nrow(p))
cor <- m
cor <- signif(cor,3)
#cor[abs(cor)<0.21] <- ''
p <- signif(p,digits = 3)
textMatrix = paste(cor,"\n(",
                   p, ")",sep = "")
dim(textMatrix) = dim(m)
colnames(textMatrix)<-colnames(m)
rownames(textMatrix)<-rownames(m)
pdf(file = '02.lncRNA-Clinical_relationships(validation).pdf',w=7,h=6)
par(mar = c(9, 8, 3, 1));
labeledHeatmap(Matrix = m, 
               xLabels = colnames(m), 
               yLabels = rownames(m), 
               cex.lab = 1, 
               # ySymbols = colnames(MEs_col)[-20], 
               colorLabels = FALSE, 
               colors = blueWhiteRed(50, 0.8),
               textMatrix = textMatrix, 
               setStdMargins = FALSE, 
               cex.text = 0.8, 
               zlim = c(-1,1),
               main = paste("lncRNA-Clinical relationships"))
dev.off()

png(file = '02.lncRNA-Clinical_relationships(validation).png',w=600,h=500)
par(mar = c(9, 8, 3, 1));
labeledHeatmap(Matrix = m, 
               xLabels = colnames(m), 
               yLabels = rownames(m), 
               cex.lab = 1, 
               # ySymbols = colnames(MEs_col)[-20], 
               colorLabels = FALSE, 
               colors = blueWhiteRed(50, 0.8),
               textMatrix = textMatrix, 
               setStdMargins = FALSE, 
               cex.text = 0.8, 
               zlim = c(-1,1),
               main = paste("lncRNA-Clinical relationships"))
dev.off()
