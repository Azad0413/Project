rm(list = ls())
setwd("/data/nas1/luchunlin/project/CD-0601-2/")
if (! dir.exists("./11_IC50")){
  dir.create("./11_IC50")
}
setwd("./11_IC50")


library(lance)
library(tidyverse)
dat<-read.delim2('../00_rawdata/dat(GSE65682).xls',row.names = 1)%>%lc.tableToNum()
# colnames(dat)<-gsub('.','-',colnames(dat),fixed = T)
group <- read.delim2('../00_rawdata/group(GSE65682).xls')
control.sample <- group$sample[which(group$group=='control')]
dat<-dat[,group$sample]

hubgene <- read.delim2('../07_features/features.xls')
model_expr<-dat[hubgene$symbol, group$sample]

library(pRRophetic)
library(ggplot2)
set.seed(12345)

drug<-read.table(file = '/data/nas1/luchunlin/pipeline/Medicinal_Sensity/drugs.txt',sep='\t',header=F)
ic50<-data.frame(group$sample)
a<-data.frame(row.names=group$sample,group$group)
colnames(a)<-'group'
cnt<-1
while (cnt < 139) {
  predictedPtype <- pRRopheticPredict(as.matrix(model_expr), drug[cnt,],selection=1)
  Tipifarnib<-data.frame(predictedPtype)
  colnames(Tipifarnib)<-drug[cnt,]
  a<-cbind(a,Tipifarnib)
  cnt = cnt + 1
}
write.table(a,'IC50.xls',sep='\t',quote=F)
b<-a
# b[b<0]<-NA
# 先写成函数的形式，方便调用
removeRowsAllNa  <- function(x){x[apply(x, 1, function(y) any(!is.na(y))),]}
removeColsAllNa  <- function(x){x[, apply(x, 2, function(y) any(!is.na(y)))]}
c<-removeColsAllNa(b)
na_flag <- apply(is.na(c), 2, sum)
x <- c[, which(na_flag == 0)]
dim(x)
medicinal_result <- t(subset(x, select = -group)) 
table(group$group)
Sepsis_group <- group$sample[which(group$group=='Sepsis')]
control_group <- group$sample[which(group$group=='control')]
pvalue = padj = log2FoldChange <- matrix(0, nrow(medicinal_result), 1)
for (i in 1:nrow(medicinal_result)){
  pvalue[i, 1] = p.value = wilcox.test(medicinal_result[i, Sepsis_group],
                                       medicinal_result[i, control_group])$p.value
  log2FoldChange[i, 1] = mean(medicinal_result[i, Sepsis_group]) - 
    mean(medicinal_result[i, control_group])
}
padj <- p.adjust(as.vector(pvalue), "fdr", n = length(pvalue))
rTable <- data.frame(log2FoldChange, 
                     pvalue, 
                     padj,
                     row.names = rownames(medicinal_result))
Sepsis_group_res <- signif(apply(medicinal_result[rownames(rTable), Sepsis_group], 
                               1,
                               median), 4)
control_group_res <- signif(apply(medicinal_result[rownames(rTable), control_group], 
                              1, 
                              median), 4)
rTable <- data.frame(Sepsis_group_res, 
                     control_group_res,
                     rTable[, c("padj", "pvalue", "log2FoldChange")])
rTable$drugs <- rownames(rTable)
rTable$sig <- ifelse(rTable$pvalue < 0.05,
                     ifelse(rTable$pvalue < 0.01, 
                            ifelse(rTable$pvalue < 0.001,
                                   ifelse(rTable$pvalue < 0.0001,
                                          paste(rTable$drugs, "****",  sep = ""),
                                          paste(rTable$drugs, "***", sep = "")),
                                   paste(rTable$drugs, "**", sep = "")),
                            paste(rTable$drugs, "*",  sep = "")), 
                     rTable$drugs)

write.table(rTable,
            file = "drugs_wilcox_test.xls",
            quote = F,
            row.names = F)
DE.drug<-rTable[which(rTable$pvalue<0.05),]
DE.drug_up <- DE.drug[which(DE.drug$log2FoldChange>0),]
##66
DE.drug_down <- DE.drug[which(DE.drug$log2FoldChange<0),]
##64
library(tidyr)
library(ggplot2)
library(ggpubr)
library(Ipaper)
all(rownames(rTable) == rownames(medicinal_result))
drugs_res <- data.frame(drugs=rownames(medicinal_result), medicinal_result, pvalue=rTable$pvalue)
drugs_res <- drugs_res[rownames(drugs_res)%in%rownames(DE.drug_down),]
# drugs_res <- drugs_res[rownames(rTable2),]

drugs_res <- drugs_res[which(rownames(drugs_res) %in% rownames(DE.drug)),]
violin_dat <- gather(drugs_res, key=indivs, value=score, -c("drugs","pvalue"))
violin_dat$indivs <- ifelse(gsub("\\.","-",violin_dat$indivs) %in% Sepsis_group,
                            "Sepsis", "control") 
violin_dat$indivs <- factor(violin_dat$indivs, levels = c("Sepsis", "control"))
head(violin_dat)
drugs_hub_boxplot1 <- ggboxplot(violin_dat, x = "indivs", y = "score",
                                color = "indivs", palette = c("#A73030FF", "#0073C2FF"),
                                add = "jitter",
                                short.panel.labs = T,
                                ggtheme = theme_bw()) +
  stat_compare_means(label = "p.signif", label.x = 1.4, vjust = 0.5)
drugs_hub_boxplot <- facet(drugs_hub_boxplot1,
                           facet.by = "drugs",
                           short.panel.labs = T,
                           panel.labs.background = list(fill = "white"),
                           ncol = 7,
                           scales = "free_y") + xlab("") + ylab("IC(50)") +
  # geom_text(data=data_text,
  #           mapping=aes(x=x,y=y,label=label),nudge_x=0.1,nudge_y=0.1)+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        # strip.background = element_blank(),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 10, face = "bold",angle = 45,hjust = 1),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        text = element_text(size = 13, face = "bold"))
drugs_hub_boxplot
ggsave(filename = "01.drugs.plot(down).pdf", height = 12, width = 14)
ggsave(filename = "01.drugs.plot(down).png", height = 13, width = 16)


drugs_res <- data.frame(drugs=rownames(medicinal_result), medicinal_result, pvalue=rTable$pvalue)
drugs_res <- drugs_res[rownames(drugs_res)%in%rownames(DE.drug_up),]
# drugs_res <- drugs_res[rownames(rTable2),]

drugs_res <- drugs_res[which(rownames(drugs_res) %in% rownames(DE.drug)),]
violin_dat <- gather(drugs_res, key=indivs, value=score, -c("drugs","pvalue"))
violin_dat$indivs <- ifelse(gsub("\\.","-",violin_dat$indivs) %in% Sepsis_group,
                            "Sepsis", "control") 
violin_dat$indivs <- factor(violin_dat$indivs, levels = c("Sepsis", "control"))
head(violin_dat)
drugs_hub_boxplot1 <- ggboxplot(violin_dat, x = "indivs", y = "score",
                                color = "indivs", palette = c("#A73030FF", "#0073C2FF"),
                                add = "jitter",
                                short.panel.labs = T,
                                ggtheme = theme_bw()) +
  stat_compare_means(label = "p.signif", label.x = 1.4, vjust = 0.5)
drugs_hub_boxplot <- facet(drugs_hub_boxplot1,
                           facet.by = "drugs",
                           short.panel.labs = T,
                           panel.labs.background = list(fill = "white"),
                           ncol = 7,
                           scales = "free_y") + xlab("") + ylab("IC(50)") +
  # geom_text(data=data_text,
  #           mapping=aes(x=x,y=y,label=label),nudge_x=0.1,nudge_y=0.1)+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        # strip.background = element_blank(),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 10, face = "bold",angle = 45,hjust = 1),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        text = element_text(size = 13, face = "bold"))
drugs_hub_boxplot
ggsave(filename = "02.drugs.plot(up).pdf", height = 12, width = 14)
ggsave(filename = "02.drugs.plot(up).png", height = 13, width = 16)

# 
# 
# 
# library(oncoPredict)
# library(data.table)
# library(gtools)
# library(reshape2)
# library(ggpubr)
# library(reshape2)
# library(ggpubr)
# 
# th=theme(axis.text.x = element_text(angle = 45,vjust = 0.5))
# dir='/data/nas1/luchunlin/GDSC/DataFiles/Training Data/'
# GDSC2_Expr = readRDS(file=file.path(dir,'GDSC2_Expr (RMA Normalized and Log Transformed).rds'))
# GDSC2_Res = readRDS(file = file.path(dir,"GDSC2_Res.rds"))
# GDSC2_Res <- exp(GDSC2_Res) 
# 
# testExpr<- as.matrix(model_expr)
# dim(testExpr)  
# 
# calcPhenotype(trainingExprData = GDSC2_Expr,
#               trainingPtype = GDSC2_Res,
#               testExprData = testExpr,
#               batchCorrect = 'eb',  #   "eb" for ComBat  
#               powerTransformPhenotype = TRUE,
#               removeLowVaryingGenes = 0.2,
#               minNumSamples = 5, 
#               printOutput = TRUE, 
#               removeLowVaringGenesFrom = 'rawData' )
# result1 <- read.csv('calcPhenotype_Output/DrugPredictions.csv',row.names = 1)
# 
# drug <- data.frame(drug=colnames(result1))
# # drug_plot1 <- result1[c('LFM.A13_192','QL.XII.47_235','QL.XII.61_1203','QL.X.138_331','Lenalidomide_1020')]
# drug_plot <- log2(result1+1)
# # colnames(drug_plot1) <- c('LFM-A13','QL-XII-47','QL-XII-61','QL-X-138','Lenalidomide')
# #LFM.A13_192	QL.XII.47_235 QL.XII.61_1203  QL.X.138_331		Lenalidomide_1020
# drug <- separate(drug,col = drug,into = c('drug','drop'),sep = '_')%>%dplyr::select(-'drop')
# colnames(drug_plot) <- drug$drug
# 
# 
# ### 差异分析--------
# drug_plot <- drug_plot%>%as.data.frame()%>%
#   rownames_to_column(var = 'sample')
# 
# drug_plot <- merge(group,drug_plot,by='sample')
# drug_plot2 <- tidyr::gather(drug_plot,Drug,IC50,-c('sample','group'))
# 
# library(rstatix)
# library(ggplot2)
# library(ggpubr)
# res.drug <- drug_plot2%>%
#   group_by(Drug)%>%
#   wilcox_test(IC50 ~ group) %>%
#   adjust_pvalue(method = 'BH')%>%
#   add_significance('p')
# 
# write.table(res.drug,file = 'res.drug.xls',sep = '\t',row.names = F,quote = F)
# DE.drug<-res.drug[which(res.drug$p<0.05),]
# # write.table(DE.drug,file = 'DE.drug.xls',sep = '\t',row.names = F,quote = F)
# colnames(drug_plot2)
# violin<-drug_plot2[drug_plot2$Drug%in%res.drug$Drug[which(res.drug$p<0.05)],]
# drug_plot <- ggplot(violin, aes(x=risk,
#                                 y=IC50,
#                                 fill=risk)) +
#   geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
#   #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
#   stat_boxplot(geom="errorbar",
#                width=0.1,
#                position = position_dodge(0.9)) +
#   geom_boxplot(width=0.3,
#                position=position_dodge(0.9),
#                outlier.shape = NA,fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
#   scale_fill_manual(values= c("#FF6A6A", "#20B2AA"), name = "Group")+
#   labs(title="", x="", y = "IC50",size=20) +
#   stat_compare_means(data = violin,
#                      mapping = aes(group = risk),
#                      label ="p.signif",
#                      method = 'wilcox.test',
#                      paired = F,label.x = 1.4) +
#   theme_bw()+
#   theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18),
#         axis.text.x=element_text(angle=0,hjust=,colour="black",face="bold",size=10), 
#         axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
#         axis.title.x=element_text(size=16,face="bold"),
#         axis.title.y=element_text(size=16,face="bold"),
#         legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
#         legend.title = element_text(face = "bold", size = 12),
#         legend.position = "top",
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank())+facet_wrap(~Drug,scales = "free",nrow = 2) +
#   guides(fill='none')
# drug_plot
# ggsave(filename = '01.drug.plot.pdf',drug_plot,w=6,h=5)
# ggsave(filename = '01.drug.plot.png',drug_plot,w=6,h=5)