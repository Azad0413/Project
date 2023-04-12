rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-385-10/")
if (! dir.exists("./09_ceRNA")){
  dir.create("./09_ceRNA")
}
setwd("./09_ceRNA")

##miRNA--------
library(GEOquery)
library(Biobase)
library(GEOquery)
library(Biobase)
gset<-getGEO("GSE41321",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[1]]))
gpl<-getGEO("GPL14767",destdir = '.')
a=gset[[1]]
gpl<-Table(gpl)    
colnames(gpl)
probe2symbol<-gpl %>%
  dplyr::select('ID','miRNA_ID')%>%
  filter('miRNA_ID'!='')%>%
  separate('miRNA_ID',c('symbol','drop'),sep = '///')%>%
  dplyr::select(-drop)
probe2symbol=probe2symbol[probe2symbol$symbol!='',]
dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symbol$ID<-as.character(probe2symbol$ID)
dat<-dat %>%
  inner_join(probe2symbol,by='ID')%>% 
  dplyr::select(-ID)%>%     ## 去除多余信息
  dplyr::select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
dat <- log2(dat+7)
pd <- pData(a)
group <- data.frame(sample=pd$geo_accession,group=pd$`cell type:ch1`)
table(group$group)
group$group <- ifelse(group$group=='Retinoblastoma','RB','control')
write.table(group,file = 'group(GSE41321).xls',sep = '\t',row.names = F,quote = F)
write.table(dat,file = 'dat(GSE41321).xls',sep = '\t',row.names = T,quote = F)

### DEG(miRNA)-------
df = dat
df.group = group[order(group$group),]
df = df[,df.group$sample]
table(df.group$group)
df.group$group = factor(df.group$group, levels = c("control", "RB"))
design.mat = cbind(control = ifelse(df.group$group == "control", 1, 0), 
                   RB = ifelse(df.group$group == "control", 0, 1))
contrast.mat = makeContrasts(contrasts="RB-control", levels=design.mat)

fit = lmFit(df, design.mat)
fit = contrasts.fit(fit, contrast.mat)
fit = eBayes(fit)
fit = topTable(fit, coef = 1, number = Inf, adjust.method = "fdr")
#fit = fit[c(1,4,5)]
DEG=na.omit(fit)
logFC_cutoff <- 1
DEG$change = as.factor(
  ifelse(DEG$adj.P.Val<0.05 & abs(DEG$logFC) > logFC_cutoff,
         ifelse(DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)
sig_diff <- subset(DEG,
                   DEG$adj.P.Val < 0.05 & abs(DEG$logFC) > logFC_cutoff)

dim(DEG)
dim(sig_diff)
## 54  6
summary(sig_diff$change)
# DOWN  NOT   UP 
#    0 54 
write.table(DEG,file = "DEG_all(GSE41321).xls",quote = F,sep = "\t",row.names = T)
write.table(sig_diff,file = "DEG_sig.xls(GSE41321)",quote = F,sep = "\t",row.names = T)

## 火山图------
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(Ipaper)
library(scales)
library(ggrepel)
dat_rep<-DEG[rownames(DEG)%in%
               rownames(head(sig_diff[order(sig_diff$logFC,decreasing = T),],10)),]
volcano_plot<- ggplot(data = DEG, 
                      aes(x = logFC,
                          y = -log10(adj.P.Val), 
                          color =change)) +
  scale_color_manual(values = c("darkgray","red")) +
  scale_x_continuous(breaks = c(-10,-6,-3,-1,0,1,3,6,10)) +
  scale_y_continuous(trans = "log1p",
                     breaks = c(0,1,5,10,20,50, 100,200)) +
  geom_point(size = 1.2, alpha = 0.4, na.rm=T) +
  theme_bw(base_size = 12, base_family = "Times") +
  geom_vline(xintercept = c(-1,1),
             lty = 4,
             col = "darkgray",
             lwd = 0.6)+
  geom_hline(yintercept = -log10(0.05),
             lty = 4,
             col = "darkgray",
             lwd = 0.6)+
  theme(legend.position = "right",
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(face="bold",
                                   color="black",
                                   family = "Times",
                                   size=13),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face = "bold",
                                   color = "black",
                                   size = 15),
        axis.text.y = element_text(face = "bold",
                                   color = "black",
                                   size = 15),
        axis.title.x = element_text(face = "bold",
                                    color = "black",
                                    size = 15),
        axis.title.y = element_text(face = "bold",
                                    color = "black",
                                    size = 15)) +
  geom_label_repel(
    data = dat_rep,
    aes(label = rownames(dat_rep)),
    max.overlaps = 20,
    size = 3,
    box.padding = unit(0.5, "lines"),
    min.segment.length = 0,
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )+
  labs(x = "log (Fold Change)",
       y = "-log10 (adj.P.Val)")
volcano_plot
ggsave('01.volcano(miRNAs).png', volcano_plot,width = 8, height = 7)
ggsave('01.volcano(miRNAs).pdf', volcano_plot,width = 8, height = 7)

## lncRNA-------
exp <- read_xlsx('GSE125903_All_Genes_FPKM.xlsx')
DEG <- read.csv('GSE125903_Differential_Expression_Genes_DESeq2.csv')%>%column_to_rownames(var = 'ID')
## 加载注释文件
library("rtracklayer")
#gtf_data = import('/data/nas1/luchunlin/pipeline/GENEANNO/gencode.v22.annotation.gtf.gz') #gtf的路径
gtf_data = import('/data/nas1/luchunlin/pipeline/GENEANNO/gencode.v41.long_noncoding_RNAs.gtf.gz') #gtf的路径
gtf_data = as.data.frame(gtf_data)
table(gtf_data$gene_type)
lncRNA=gtf_data%>%
  dplyr::filter(type=="gene",gene_type=="lncRNA")%>%
  dplyr::select(gene_id,gene_type,gene_name)
exp <- exp[exp$gene_short_name%in%lncRNA$gene_name,]
exp <- exp[!duplicated(exp$gene_short_name),]
write.table(dat,file = 'dat(GSE125903).xls',sep = '\t',row.names = T,quote = F)
##1675
DEG <- DEG[exp$gene_id,]
rownames(DEG) <- exp$gene_short_name
dat <- exp[,-2]%>%column_to_rownames(var = 'gene_short_name')
logFC_cutoff <- 1
DEG$change = as.factor(
  ifelse(DEG$padj<0.05 & abs(DEG$log2FoldChange) > logFC_cutoff,
         ifelse(DEG$log2FoldChange > logFC_cutoff ,'UP','DOWN'),'NOT')
)
sig_diff <- subset(DEG,
                   DEG$padj < 0.05 & abs(DEG$log2FoldChange) > logFC_cutoff)

dim(DEG)
dim(sig_diff)
##83
summary(sig_diff$change)
# DOWN  NOT   UP 
# 24    0   59
write.table(DEG,file = "DEG_all(GSE125903).xls",quote = F,sep = "\t",row.names = T)
write.table(sig_diff,file = "DEG_sig.xls(GSE125903)",quote = F,sep = "\t",row.names = T)
## 火山图------
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(Ipaper)
library(scales)
library(ggrepel)
dat_rep<-DEG[rownames(DEG)%in%
               rownames(rbind(head(sig_diff[order(sig_diff$log2FoldChange,decreasing = T),],10),
                              head(sig_diff[order(sig_diff$log2FoldChange,decreasing = F),],10))),]
volcano_plot<- ggplot(data = DEG, 
                      aes(x = log2FoldChange,
                          y = -log10(padj), 
                          color =change)) +
  scale_color_manual(values = c("blue","darkgray","red")) +
  scale_x_continuous(breaks = c(-10,-6,-3,-1,0,1,3,6,10)) +
  scale_y_continuous(trans = "log1p",
                     breaks = c(0,1,5,10,20,50, 100,200)) +
  geom_point(size = 1.2, alpha = 0.4, na.rm=T) +
  theme_bw(base_size = 12, base_family = "Times") +
  geom_vline(xintercept = c(-1,1),
             lty = 4,
             col = "darkgray",
             lwd = 0.6)+
  geom_hline(yintercept = -log10(0.05),
             lty = 4,
             col = "darkgray",
             lwd = 0.6)+
  theme(legend.position = "right",
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(face="bold",
                                   color="black",
                                   family = "Times",
                                   size=13),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face = "bold",
                                   color = "black",
                                   size = 15),
        axis.text.y = element_text(face = "bold",
                                   color = "black",
                                   size = 15),
        axis.title.x = element_text(face = "bold",
                                    color = "black",
                                    size = 15),
        axis.title.y = element_text(face = "bold",
                                    color = "black",
                                    size = 15)) +
  geom_label_repel(
    data = dat_rep,
    aes(label = rownames(dat_rep)),
    max.overlaps = 20,
    size = 3,
    box.padding = unit(0.5, "lines"),
    min.segment.length = 0,
    point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )+
  labs(x = "log (Fold Change)",
       y = "-log10 (adj.P.Val)")
volcano_plot
ggsave('02.volcano(lncRNAs).png', volcano_plot,width = 8, height = 7)
ggsave('02.volcano(lncRNAs).pdf', volcano_plot,width = 8, height = 7)


### miRNA
hubgene <- read.delim2('../06_model/hubgene.xls')

##mirWalk
PDE8B <- read_csv('PDE8B-miR.csv')
PDE8B <- PDE8B[which(PDE8B$bindingp>0.95 & PDE8B$energy < -15),]
PDE8B <- PDE8B[!duplicated(PDE8B$mirnaid),]
SPRY2 <- read_csv('SPRY2-miR.csv')
SPRY2 <- SPRY2[which(SPRY2$bindingp>0.95 & SPRY2$energy < -15),]
SPRY2 <- SPRY2[!duplicated(SPRY2$mirnaid),]
ESRRB <- read_csv('ESRRB-miR.csv')
ESRRB <- ESRRB[which(ESRRB$bindingp>0.95 & ESRRB$energy < -15),]
ESRRB <- ESRRB[!duplicated(ESRRB$mirnaid),]
gene2mir <- rbind(PDE8B,SPRY2,ESRRB)
gene2mir <- gene2mir[!duplicated(gene2mir$mirnaid),]

diffmir <- read.delim2('DEG_sig.xls(GSE41321)')
cemir <- data.frame(symbol=intersect(rownames(diffmir),gene2mir$mirnaid))
## 9
gene2mir2 <- rbind(PDE8B,SPRY2,ESRRB)
gene2mir2 <- gene2mir2[gene2mir2$mirnaid%in%cemir$symbol,]
gene2mir2 <- gene2mir2[,c(3,1)]
colnames(gene2mir2) <- c('symbol','miRNA')
write.table(gene2mir2,file = 'gene2mir.xls',sep = '\t',row.names = F,quote = F)


##celncRNA

mir1 <- read_xlsx('hsa-miR-146b-5p.xlsx')
mir1 <- mir1[!duplicated(mir1$geneName),]
mir2 <- read_xlsx('hsa-miR-188-5p.xlsx')
mir2 <- mir2[!duplicated(mir2$geneName),]
mir3 <- read_xlsx('hsa-miR-342-3p.xlsx')
mir3 <- mir3[!duplicated(mir3$geneName),]
mir4 <- read_xlsx('hsa-miR-665.xlsx')
mir4 <- mir4[!duplicated(mir4$geneName),]
mir2lnc <- rbind(mir1,mir2,mir3,mir4)
mir2lnc <- mir2lnc[!duplicated(mir2lnc$geneName),]
difflnc <- read.delim2('DEG_sig.xls(GSE125903)')
celnc <- data.frame(symbol=intersect(rownames(difflnc),mir2lnc$geneName))
##13
mir2lnc2 <- rbind(mir1,mir2,mir3,mir4)
mir2lnc2 <- mir2lnc2[mir2lnc2$geneName%in%celnc$symbol,]
mir2lnc2 <- mir2lnc2[,c(4,2)]
colnames(mir2lnc2) <- c('miRNA','lncRNA')
write.table(mir2lnc2,'mir2lnc.xls',sep = '\t',row.names = F,quote = F)
