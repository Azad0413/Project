## hsa-486-3p  hsa-miR-941 调控基因， 与脂质代谢相关基因取交集后富集分析。

setwd("/data/nas1/luchunlin/project/BJZK-107/")
library(readr)
library(readxl)

## mirBase --------
## 486-3p
mirbase<-read_xlsx('工作簿3.xlsx')
mRNA<-mirbase$mRNA
gene_transform <- bitr(mRNA,
                       fromType = "REFSEQ",
                       toType = "SYMBOL",
                       OrgDb = "org.Hs.eg.db")
mirbase<-gene_transform$SYMBOL%>%as.data.frame()
mirbase_mir_486_3p<-read_xlsx('mirBase_miR-486-3p.xlsx')
mirbase_mir_486_3p<-mirbase_mir_486_3p$`Gene Symbol`%>%
  as.data.frame()
mirbase_mir_486_3p<-rbind(mirbase,mirbase_mir_486_3p)
mirbase_mir_486_3p<-mirbase_mir_486_3p[!duplicated(mirbase_mir_486_3p$.),]
mirbase_mir_486_3p<-as.data.frame(mirbase_mir_486_3p)
colnames(mirbase_mir_486_3p)<-'target'
## 2036
write.table(mirbase_mir_486_3p,file = 'mirbase_mir_486_3p.xls',
            row.names = F,
            sep = '\t')
## 941没有

## targetscan------
library(tidyverse)
## 486-3p
targetscan_mir_485_3p<-read_xlsx('TargetScan8.0__miR-486-3p.predicted_targets.xlsx')
## 5264
colnames(targetscan_mir_485_3p)
targetscan_mir_485_3p<-select(targetscan_mir_485_3p,"Target gene")

write.table(targetscan_mir_485_3p,file = 'targetscan_mir_485_3p.xls',
            row.names = F,
            sep = '\t')
## mir-941
targetscan_mir_941<-read_xlsx('TargetScan8.0__miR-941.predicted_targets.xlsx')
##1056
targetscan_mir_941<-select(targetscan_mir_941,"Target gene")
write.table(targetscan_mir_941,file = 'targetscan_mir_941.xls',
            row.names = F,
            sep = '\t')

## mirwalk-----
## 486-3p
mirwalk_mir_486_3p<-read.csv(file = 'miRWalk_miRNA_Targets_miR-486-3p.csv')
mirwalk_mir_486_3p<-mirwalk_mir_486_3p[!duplicated(mirwalk_mir_486_3p$genesymbol),]
##3998
mirwalk_mir_486_3p<-select(mirwalk_mir_486_3p,'genesymbol')
write.table(mirwalk_mir_486_3p,file = 'miRWalk_mir_486_3p.xls',
            row.names = F,
            sep = '\t')
## mir-941
mirwalk_mir_941<-read.csv(file = 'miRWalk_miRNA_Targets.csv')
mirwalk_mir_941<-mirwalk_mir_941[!duplicated(mirwalk_mir_941$genesymbol),]
## 2563
mirwalk_mir_941<-select(mirwalk_mir_941,'genesymbol')
write.table(mirwalk_mir_941,file = 'miRWalk_mir_941.xls',
            row.names = F,
            sep = '\t')
## 取交集------
## mir-486-3P
mir_486_3p<-mirbase_mir_486_3p[mirbase_mir_486_3p$target%in%targetscan_mir_485_3p$`Target gene`&mirbase_mir_486_3p$target%in%mirwalk_mir_486_3p$genesymbol,]%>%
  as.data.frame()
colnames(mir_486_3p)<-'target'
write.table(mir_486_3p,file = 'mir_486_3p_target.xls',
            row.names = F,
            sep = '\t')

library(ggvenn)
mydata1<-list('miRbase'=mirbase_mir_486_3p$target,'TargetScan'=targetscan_mir_485_3p$`Target gene`,'miRWalk'=mirwalk_mir_486_3p$genesymbol)
ggvenn(mydata1,c('miRbase','TargetScan','miRWalk'),
       fill_color = c("#ffb2b2","#b2e7cb","#F0E442"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb","#F0E442"),
       text_color = 'black')


## 454
## mir-941
mir_941<-mirwalk_mir_941[mirwalk_mir_941$genesymbol%in%targetscan_mir_941$`Target gene`,]%>%as.data.frame()
colnames(mir_941)<-'target'
## 220
write.table(mir_941,file = 'mir_941_target.xls',
            row.names = F,
            sep = '\t')

mydata2<-list('TargetScan'=targetscan_mir_941$`Target gene`,'miRWalk'=mirwalk_mir_941$genesymbol)
ggvenn(mydata2,c('TargetScan','miRWalk'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
## 脂肪代谢相关基因-----
geneset<-read_xlsx('GENESET.xlsx')
geneset<-geneset[!duplicated(geneset$Geneset),]
## 1045

### 脂肪代谢相关target--------
### 486-3p
lipid_target1<-geneset[geneset$Geneset%in%mir_486_3p$target,]
colnames(lipid_target1)<-'target'
### 19
write.table(lipid_target1,file = 'lipid_target1.xls',
            row.names = F,
            sep = '\t')

mydata3<-list("Target"=mir_486_3p$target,"Lipid metabolism"=geneset$Geneset)
ggvenn(mydata3,c('Target','Lipid metabolism'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
### mir-941
lipid_target2<-geneset[geneset$Geneset%in%mir_941$genesymbol,]
## 9
colnames(lipid_target2)<-'target'
write.table(lipid_target2,file = 'lipid_target2.xls',
            row.names = F,
            sep = '\t')
mydata4<-list("Target"=mir_941$genesymbol,"Lipid metabolism"=geneset$Geneset)
ggvenn(mydata4,c('Target','Lipid metabolism'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')

## 合在一起
lipid_target<-rbind(lipid_target1,lipid_target2)
## 2个重复 EMD22 DGKH
##GO/KEGG
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(GOplot)
gene_transform <- bitr(lipid_target$Geneset,
                       fromType = "SYMBOL",
                       toType = c("ENTREZID"),
                       OrgDb = "org.Hs.eg.db")

ego <- enrichGO(gene = gene_transform$ENTREZID,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                readable = TRUE)
write.table(ego,file = "GO.xls",sep = "\t",quote = F,row.names = F)
# 展示富集最显著的 GO term
go_bar <- barplot(ego, showCategory=5, split="ONTOLOGY") +
  facet_grid(ONTOLOGY ~ ., scales = "free")
go_bar
##05-2 KEGG富集分析（气泡图）
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "hsa",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk@result,file = "KEGG.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk, showCategory=15)
kk_dot
