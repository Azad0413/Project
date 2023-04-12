# rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/BJTC-236")
if (! dir.exists("./00_rawdata")){
  dir.create("./00_rawdata")
}
setwd("./00_rawdata")
## 读取从xena下载的数据
expr<-read_tsv(file = 'TCGA-BRCA.htseq_counts.tsv')
expr<-as.data.frame(expr)
rownames(expr)<-expr[,1]
expr<-expr[,-1]
expr_fpkm<-read_tsv(file = 'TCGA-BRCA.htseq_fpkm.tsv')
expr_fpkm<-as.data.frame(expr_fpkm)
rownames(expr_fpkm)<-expr_fpkm[,1]
expr_fpkm<-expr_fpkm[,-1]
## xena下载的数据经过了log2+1转化，需要将其还原
expr<-2^expr-1
expr_fpkm<-2^expr_fpkm-1
## 对数据进行id转化
genecode<-read.table(file = 'gencode.v22.annotation.gene.probeMap')
probe2symbol<-genecode[,(1:2)]
colnames(probe2symbol)<-c('id','symbol')
expr$id<-rownames(expr)
dat<-merge(expr,probe2symbol,by='id')
dat<-dat[,c(1,1219,2:1218)]
dat<-dat[!duplicated(dat$symbol),]
rownames(dat)<-dat$symbol
dat<-dat[,-c(1,2)]

expr_fpkm$id<-rownames(expr_fpkm)
dat_fpkm<-merge(expr_fpkm,probe2symbol,by='id')
dat_fpkm<-dat_fpkm[,c(1,1219,2:1218)]
dat_fpkm<-dat_fpkm[!duplicated(dat_fpkm$symbol),]
rownames(dat_fpkm)<-dat_fpkm$symbol
dat_fpkm<-dat_fpkm[,-c(1,2)]
dat_fpkm<-dat_fpkm[rownames(dat_fpkm)%in%rownames(dat_final),]
dataNorm.BRCA <- TCGAanalyze_Normalization(tabDF = dat,
                                           geneInfo = geneInfo,
                                           method = "gcContent")
# 将标准化后的数据再过滤，得到最终的数据
dataFilt.BRCA.final <- TCGAanalyze_Filtering(tabDF = dataNorm.BRCA,
                                             method = "quantile", 
                                             qnt.cut =  0.25)
dim(dataFilt.BRCA.final)
# [1] 13627  1217
## 筛选癌症组织，去掉癌旁组织。01-09为肿瘤，10-19为正常对照
mete=data.frame(colnames(dataFilt.BRCA.final))  # 取第一行样本id
for (i in 1:length(mete[,1])) {
  num=as.numeric(as.character(substring(mete[i,1],14,15)))
  if(num %in% seq(1,9)){mete[i,2]="T"}
  if(num %in% seq(10,29)){mete[i,2]="N"}
}
names(mete)=c("id","group")
mete$group=as.factor(mete$group)
mete=subset(mete,mete$group=="T")
exp_tumor<-dataFilt.BRCA.final[,which(colnames(dataFilt.BRCA.final)%in%mete$id)]
exp_tumor<-as.data.frame(exp_tumor)
exp_tumor_fpkm<-dat_fpkm[,which(colnames(dat_fpkm)%in%mete$id)]
exp_tumor_fpkm<-exp_tumor_fpkm[rownames(exp_tumor_fpkm)%in%rownames(exp_tumor),]
# [1104]
exp_control<-dataFilt.BRCA.final[,which(!colnames(dataFilt.BRCA.final)%in%mete$id)]
exp_control<-as.data.frame(exp_control)
# [113]
write.table(exp_tumor,file = "tumor_exp.xls",
            quote = F,
            sep = "\t",
            row.names = T)
dat_final<-cbind(exp_control,exp_tumor)
write.table(dat_final,file = "exp.xls",
            quote = F,
            sep = "\t",
            row.names = T)
## 分组
group<-data.frame(sample=colnames(dat_final),
                  group=c(rep('Normal',113),rep('Tumor',1104)))
## 下载临床数据
clinical<-GDCquery_clinic(project = "TCGA-BRCA",type = "clinical")
# 02 差异分析----------
library(DESeq2)
setwd("/data/nas1/luchunlin/project/BJTC-236")
if (! dir.exists("./01_DEGs")){
  dir.create("./01_DEGs")
}
setwd("./01_DEGs")
##02-1 差异基因鉴定
# 构建dds矩阵
##样本分组信息
colData<-group
colData$group<-factor(colData$group,levels = c("Normal","Tumor"))

write.table(colData,
            file = "group.xls",
            quote = F,
            row.names = F)
dds<-DESeqDataSetFromMatrix(countData = dat_final,colData=colData,design = ~group)
dds = dds[rownames(counts(dds)) > 1,]
## 对原始dds进行normalize
dds<-DESeq(dds)
## 显示dds信息
dds
# 提取DESeq2分析结果
## 使用Result函数提取差异分析结果。
## 将提取的差异分析结果定义为变量“res”。
## contrast：定义谁和谁比较
res =results(dds, contrast = c("group","Normal","Tumor"))
## 对结果res利用order（）函数按照pvalue值进行排序。
res =res[order(res$pvalue),]
head(res)
summary(res)
## 保存所有输出结果
write.csv(res,file="All_results.csv")
# 显示显著差异基因数目
table(res$padj<0.05)
# 对显著性差异结果进行提取和保存
## 获取padj< 0.05，表达倍数取以2为对数后大于1或者小于-1的差异表达基因。
## 使用subset（）函数过滤需要结果至变量DEG中
## Usage:subset(x,...),x为objects，...为筛选的参数或者条件
DEG <- subset(res, padj < 0.05 & abs(log2FoldChange) >0.5 )
DEG<-as.data.frame(res)
DEG<-na.omit(DEG)
## 使用dim查看该结果的维度、规模
dim(DEG)
head(DEG)
## 添加change列
logFC_cutoff<-0.5
DEG$change=as.factor(
  ifelse(DEG$padj<0.05&abs(DEG$log2FoldChange)>logFC_cutoff,
         ifelse(DEG$log2FoldChange>logFC_cutoff,'UP','DOWN'),'NOT'))
table(DEG$change)
#  DOWN   NOT    UP  
#3921 7163 2543
sig_diff <- subset(DEG,
                   DEG$padj < 0.05 & abs(DEG$log2FoldChange) >= logFC_cutoff)
## 11422
DEG_write <- cbind(GeneSymbol=rownames(DEG), DEG)
write.table(DEG_write, file = "DEG_all.xls",
            quote = F,
            sep = "\t",
            row.names = F)
sig_diff_write <- cbind(GeneSymbol=rownames(sig_diff), sig_diff)
write.table(sig_diff_write, file = "DEG_sig.xls",
            quote = F,
            sep = "\t",
            row.names = F)
# 02-1 绘制火山图-----
library(ggplot2)
library(ggthemes)
library(RColorBrewer)
library(Ipaper)
library(scales)
library(ggrepel)
library(EnhancedVolcano)
dat_rep<-DEG[rownames(DEG)%in%
               c(head(rownames(subset(sig_diff,sig_diff$log2FoldChange>5)),10),
                 head(rownames(subset(sig_diff,sig_diff$log2FoldChange< -7)),10)),]

volcano_plot<-EnhancedVolcano(DEG,lab = rownames(DEG),
                              x='log2FoldChange',
                              y='padj',
                              pCutoff = 0.05,
                              FCcutoff = 0.5,
                              maxoverlapsConnectors = Inf,
                              legendLabels = c('NS','Log2FC','adj.P.Value',
                                               'adj.P.Value & Log2FC'),
                              selectLab = c(rownames(dat_rep)),
                              drawConnectors = T,
                              boxedLabels = T,
)+
  labs(y = "-log10 (adj.P.Value)")
volcano_plot
ggsave('Volcano.pdf',volcano_plot,w=7,h=6)
ggsave('Volcano.png',volcano_plot,w=7,h=6)
# 02-2 绘制热图-----
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library(gplots)
group_rt<-group$group
group_rt<-as.data.frame(group_rt)
rt<-dat_final
colnames(group_rt)<-'group'
rownames(group_rt)<-group$sample
heat<-rt[rownames(rt)%in%
           c(head(rownames(subset(sig_diff,sig_diff$log2FoldChange>5)),10),head(rownames(subset(sig_diff,sig_diff$log2FoldChange< -5)),10)),]
x<-log2(heat+1)
#x<-t(scale(t(heat)))
ann_colors<-list(
  Group = c(Normal="lightblue",Tumor="darkorange"))
pdf(file = 'heatmap.pdf',w=8,h=7)
pheatmap(mat=x,
         annotation_col = group_rt,
         color=bluered(50),
         scale = "row",
         annotation_colors = ann_colors,
         fontsize = 10,
         show_colnames = F,
         cluster_cols = F,
         cluster_rows = T)
dev.off()

png(file = 'heatmap.png',w=700,h=600)
pheatmap(mat=x,
         annotation_col = group_rt,
         color=bluered(50),
         scale = "row",
         annotation_colors = ann_colors,
         fontsize = 10,
         show_colnames = F,
         cluster_cols = F,
         cluster_rows = T)
dev.off()
# 03 差异坏死性凋亡基因necroptosis------
setwd("/data/nas1/luchunlin/project/BJTC-236")
if (! dir.exists("./02_necroptosis")){
  dir.create("./02_necroptosis")
}
setwd("./02_necroptosis")
necroptosis<-read_xlsx('/data/nas1/luchunlin/project/BJTC-236/02_necroptosis/necroptosis2.xlsx')
necroptosis<-necroptosis[!duplicated(necroptosis$Geneset),]
## 133
all_gene<-dat_final[rownames(dat_final)%in%necroptosis$Geneset,]
DENPs<-sig_diff[rownames(sig_diff)%in%necroptosis$Geneset,]
# 40
write.table(DENPs,file = 'diff_necroptosis.xls',
            sep = '\t',
            quote = F,
            row.names = T)
# 04 GO/KEGG富集-------

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(GOplot)
#PRKN
diff_gene_names <- rownames(DENPs)
gene_transform <- bitr(diff_gene_names,
                       fromType = "SYMBOL",
                       toType = c("ENTREZID", "ENSEMBL", "REFSEQ"),
                       OrgDb = "org.Hs.eg.db")
ego <- enrichGO(gene = gene_transform$ENTREZID,
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID",
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05,
                qvalueCutoff = 0.05,
                readable = TRUE)
write.table(ego,file = "GO.xls",sep = "\t",quote = F,row.names = F)
# 展示富集最显著的 GO term
go_bar <- barplot(ego, showCategory=5, split="ONTOLOGY",label_format = 50) +
  facet_grid(ONTOLOGY ~ ., scales = "free")
go_bar
ggsave(filename = 'GO_plot.pdf',go_bar,w=9,h=7)
ggsave(filename = 'GO_plot.png',go_bar,w=9,h=7)
## KEGG富集分析（气泡图）
kk <- enrichKEGG(gene = gene_transform$ENTREZID,
                 keyType = "kegg",
                 organism = "hsa",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.05,
                 qvalueCutoff = 0.05)
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write.table(kk@result,file = "KEGG.xls",sep = "\t",quote = F,row.names = F)
kk_dot <- dotplot(kk,showCategory=25,label_format=50)
kk_dot
ggsave(filename = 'KEGG_plot.pdf',kk_dot,w=9,h=7)
ggsave(filename = 'KEGG_plot.png',kk_dot,w=9,h=7)
# 05 风险模型构建-----
## 05-1 单因素cox回归----
setwd("/data/nas1/luchunlin/project/BJTC-236")
if (! dir.exists("./03_univariate_cox")){
  dir.create("./03_univariate_cox")
}
setwd("./03_univariate_cox")
survival<-read_tsv(file = 'TCGA-BRCA.survival.tsv')
survival_cancer<-survival[survival$sample%in%mete$id,]
survival_dat<-t(exp_tumor_fpkm)
survival_dat<-survival_dat[rownames(survival_dat)%in%survival_cancer$sample,]
## 1082个匹配的临床信息
## 合并生存数据
train_dat<-survival_dat[,colnames(survival_dat)%in%rownames(DENPs)]
#train_dat<-log2(train_dat)+1
train_dat<-t(scale(t(train_dat)))
train_dat<-as.data.frame(train_dat)
train_dat$sample<-rownames(train_dat)
train_dat<-merge(survival_cancer,train_dat,by='sample')
rownames(train_dat)<-train_dat$sample
train_dat<-train_dat[,-c(1,3)]
colnames(train_dat)
### 按照生存状态（OS）进行三七分组。七分作为训练集，三分作为验证集。
library(caret)
set.seed(27)
###27  9 ###16 6 ## 13 11 #28  #17  
expr<-train_dat
expr$sample<-rownames(expr)
sam<-createDataPartition(expr$OS,p= .7,list = F)
train_sample<-expr$sample[sam]
train_sample<-as.data.frame(train_sample)
test_sample<-expr$sample[-sam]
test_sample<-as.data.frame(test_sample)
train_data<-train_dat[rownames(train_dat)%in%train_sample$train_sample,]
test_data<-train_dat[rownames(train_dat)%in%test_sample$test_sample,]

### 单因素cox
library(survival)
library(survminer)
colnames_sum <- colnames(train_data)
colnames_sum <- gsub("-","_",colnames_sum)
colnames_sum <- gsub(" ","_",colnames_sum)
colnames(train_data) <- colnames_sum
covariates <- colnames_sum[-which(colnames_sum %in% c("OS", "OS.time"))]
#Surv()函数产生一个生存对象  生存时间对生存的影响 对每一个变量构建生存分析公式
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste("Surv(OS.time, OS)~", x)))  #as.formula(). 将字符串转换成公式。构建formula对象
# coxph函数用于计算cox模型 循环对每一个特征做cox回归分析
univ_models <- lapply(univ_formulas,
                      function(x) {coxph(x, data = train_data)})

univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         #获取p值
                         p.value<-signif(x$wald["pvalue"], digits=3)
                         #获取HR
                         HR <-signif(x$coef[2], digits=3);
                         #获取95%置信区间
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"],3)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],3)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", 
                                      HR.confint.upper, ")")
                         res<-c(p.value,HR)
                         names(res)<-c("p.value","HR (95% CI for HR)")
                         return(res)
                       })

## coef是公式中的回归系数b（有时候也叫beta值）。exp(coef)是cox模型中的风险比（HR）
## z代表wald统计量，是coef除以其标准误se(coef)。ower .95 upper .95则是exp(coef)的95%置信区间，可信区间越窄，可信度越高，你的实验越精确，越是真理。
res_mod <- t(as.data.frame(univ_results, check.names = FALSE))
res_mod <- as.data.frame(res_mod)
res_results_0.05 <- res_mod[which(as.numeric(res_mod$p.value) < 0.05),]
res_results_0.05 <- na.omit(res_results_0.05)
write.table(res_results_0.05,
            file = "univariate_cox_result_0.05.xls",
            quote = F,
            row.names = T,
            sep = '\t')
dim(res_results_0.05)
library(tidyr)
res_results_0.05_2 <- separate(res_results_0.05, "HR (95% CI for HR)",
                               into = c("HR", "HR.95L", "HR.95H"),
                               sep = " ")
res_results_0.05_2 <- separate(res_results_0.05_2, "HR.95L",
                               into = c("HR.95L", "HR.95H"),
                               sep = "\\-")
res_results_0.05_2$HR.95L <- gsub("\\(", "", res_results_0.05_2$HR.95L)
res_results_0.05_2$HR.95H <- gsub("\\)", "", res_results_0.05_2$HR.95H)

res_results_0.05_2[,1:ncol(res_results_0.05_2)] <- as.numeric(unlist(res_results_0.05_2[,1:ncol(res_results_0.05_2)]))
res_results_0.05_2 <- res_results_0.05_2[order(res_results_0.05_2$HR),]
hz <- paste(round(res_results_0.05_2$HR,3),
            "(",round(res_results_0.05_2$HR.95L,3),
            "-",round(res_results_0.05_2$HR.95H,3),")",sep = "")


tabletext <- cbind(c(NA,"Gene",rownames(res_results_0.05_2)),
                   c(NA,"P value",ifelse(res_results_0.05_2$p.value<0.001,
                                         "< 0.001",
                                         round(res_results_0.05_2$p.value,3))),
                   c(NA,"Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1
library(forestplot)
pdf(file = "univariate_cox_forest.pdf", height = 6, width = 10, onefile = F)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE, TRUE,rep(FALSE, 57)),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,NA,res_results_0.05_2$HR),
           lower=c(NA,NA,res_results_0.05_2$HR.95L), #95%置信区间下限
           upper=c(NA,NA,res_results_0.05_2$HR.95H), #95%置信区间上限
           boxsize=0.2,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(0, 1, 100, 400), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(0.8,"cm"), #固定行高
           graphwidth = unit(.5,"npc"), #图在表中的宽度比例
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines = list("3" = gpar(col = "black", lty = 1, lwd = 2)),
           # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
           #                 "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
           #                 "59" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
           # mar=unit(rep(0.01, times = 4), "cm"),#图形页边距
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Hazard Ratio",
           grid = T) # 垂直于x轴的网格线，对应每个刻度
dev.off()
png(filename = "univariate_cox_forest.png", height = 600, width = 800)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE, TRUE,rep(FALSE, 57)),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,NA,res_results_0.05_2$HR),
           lower=c(NA,NA,res_results_0.05_2$HR.95L), #95%置信区间下限
           upper=c(NA,NA,res_results_0.05_2$HR.95H), #95%置信区间上限
           boxsize=0.2,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           # xticks = c(0, 1, 2, 4, 6, 8), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(0.8,"cm"), #固定行高
           graphwidth = unit(.5,"npc"), #图在表中的宽度比例
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines = list("3" = gpar(col = "black", lty = 1, lwd = 2)),
           # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
           #                 "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
           #                 "59" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
           # mar=unit(rep(0.01, times = 4), "cm"),#图形页边距
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Hazard Ratio",
           grid = T) # 垂直于x轴的网格线，对应每个刻度

dev.off()
## 05-2 Lasso回归-----
setwd("/data/nas1/luchunlin/project/BJTC-236")
if (! dir.exists("./04_Lasso")){
  dir.create("./04_Lasso")
}
setwd("./04_Lasso")
library(glmnet)
x_all <- subset(train_data, select = -c(OS, OS.time))
x_all <- x_all[,rownames(res_results_0.05)]
y_all <- subset(train_data, select = c(OS, OS.time))
# 拟合模型
fit <- glmnet(as.matrix(x_all), Surv(y_all$OS.time,y_all$OS), 
              family = "cox") 
#dev.new()
#png(filename = "lasso_model.png", height = 450, width = 600)
plot(fit, xvar = "lambda",label = TRUE, las=1)
# dev.off()
pdf(file = "lasso_model.pdf", height = 6)
plot(fit, xvar = "lambda",label = TRUE, las=1)
dev.off()
# 交叉验证拟合模型
set.seed(6)
cvfit = cv.glmnet(as.matrix(x_all),
                  Surv(y_all$OS.time,y_all$OS),nfold=50,
                  family = "cox") 

#png(filename = "lasso_verify.png", height = 450, width = 600)
plot(cvfit, las =1)
#dev.off()
#pdf(file = "lasso_verify.pdf", height = 5)
#plot(cvfit, las =1)
#dev.off()

# 提取指定lambda时特征的系数
coef.min = coef(cvfit, s = "lambda.min")  ## lambda.min & lambda.1se 取一个
cvfit$lambda.min
# [1] 0.003887636
# 找出那些回归系数没有被惩罚为0的
active.min = which(coef.min@i != 0)

# 提取基因名称
lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i+1]
lasso_geneids
# [1] "BAK1"     "HEBP2"    "HSP90AA1" "TRAF2"    "TNIP1"    "JAK1"     "MAP3K5"   "IRF9"   
write(lasso_geneids, "lasso_genes.csv")
write.csv(x_all,file = "Lasso_x.csv",quote = F)
write.csv(y_all,file = "Lasso_y.csv",quote = F)
## 05-3 多因素cox回归--------
setwd("/data/nas1/luchunlin/project/BJTC-236")
if (! dir.exists("./05_Multivariate_cox")){
  dir.create("./05_Multivariate_cox")
}
setwd("./05_Multivariate_cox")
library(survival)
library(survminer)

#gene_list<-rownames(res_results_0.05)
gene_list <- lasso_geneids
cox_data <- as.formula(paste0('Surv(OS.time, OS)~',
                              paste(gene_list,
                                    sep = '',
                                    collapse = '+')))

cox_more <- coxph(cox_data,
                  data = train_data)
cox_zph <- cox.zph(cox_more)
cox_table <- cox_zph$table[-nrow(cox_zph$table),]
cox_formula <- as.formula(paste("Surv(OS.time, OS)~",
                                paste(rownames(cox_table)[cox_table[,3]>0.05],
                                      collapse = "+")))
cox_more_2 <- coxph(cox_formula, data = train_data)
# 计算共线性
cox_correlation <- cor(train_data[, rownames(cox_table)[cox_table[,3]>0.05]],
                       method = "pearson")
library(GGally)
cox_corr <- ggpairs(train_data[, rownames(cox_table)[cox_table[,3]>0.05]],
                    axisLabels = "show") +
  theme_bw() +
  theme(panel.background = element_rect(color = "black",
                                        size = 1,
                                        fill = "white"),
        panel.grid = element_blank())
cox_corr
# 评估共线性
library(rms)
vif <- vif(cox_more_2)
#some people said if the square root of VIF >2, they might be co-linear
sqrt(vif) < 2

# 多因素森林图
x <- summary(cox_more_2)
#获取p值
p.value<-signif(as.matrix(x$coefficients)[,5],2)
#获取HR
HR <-signif(as.matrix(x$coefficients)[,2],2)
#获取95%置信区间
HR.confint.lower <- signif(x$conf.int[,3],2)
HR.confint.upper <- signif(x$conf.int[,4],2)
multi_res=data.frame(p.value=p.value,
                     'HR (95% CI for HR)'=paste(HR," (",HR.confint.lower,"-",HR.confint.upper,")",sep=""),
                     stringsAsFactors = F,
                     check.names = F)
multi_res
write.csv(multi_res,
          file = "multivariate_cox_result.csv",
          quote = F,
          row.names = T)

library(tidyr)
res_results_0.05_2 <- multi_res
res_results_0.05_2 <- separate(multi_res, "HR (95% CI for HR)",
                               into = c("HR", "HR.95L", "HR.95H"),
                               sep = " ")
res_results_0.05_2 <- separate(res_results_0.05_2, "HR.95L",
                               into = c("HR.95L", "HR.95H"),
                               sep = "\\-")
res_results_0.05_2$HR.95L <- gsub("\\(", "", res_results_0.05_2$HR.95L)
res_results_0.05_2$HR.95H <- gsub("\\)", "", res_results_0.05_2$HR.95H)

res_results_0.05_2[,1:ncol(res_results_0.05_2)] <- as.numeric(unlist(res_results_0.05_2[,1:ncol(res_results_0.05_2)]))

res_results_0.05_2 <- res_results_0.05_2[order(res_results_0.05_2$HR),]

hz <- paste(round(res_results_0.05_2$HR,3),
            "(",round(res_results_0.05_2$HR.95L,3),
            "-",round(res_results_0.05_2$HR.95H,3),")",sep = "")


tabletext <- cbind(c(NA, "Gene",rownames(res_results_0.05_2)),
                   c(NA,"P value",ifelse(res_results_0.05_2$p.value<0.001,
                                         "< 0.001",
                                         round(res_results_0.05_2$p.value,3))),
                   c(NA,"Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1
library(forestplot)
#pdf(file = "multivariate_cox_forest.pdf", height = 4.5, width = 10, onefile = F)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE, TRUE,rep(FALSE, 57)),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,NA,res_results_0.05_2$HR),
           lower=c(NA,NA,res_results_0.05_2$HR.95L), #95%置信区间下限
           upper=c(NA,NA,res_results_0.05_2$HR.95H), #95%置信区间上限
           boxsize=0.1,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(0,50,100,150,200,250), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(1.2,"cm"), #固定行高
           graphwidth = unit(.5,"npc"), #图在表中的宽度比例
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines = list("3" = gpar(col = "black", lty = 1, lwd = 2)),
           # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
           #                 "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
           #                 "59" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
           # mar=unit(rep(0.01, times = 4), "cm"),#图形页边距
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Hazard Ratio",
           grid = T) # 垂直于x轴的网格线，对应每个刻度
#C指数即一致性指数，用来评价模型的预测能力。c指数是指所有病人对子中预测结果与实际结果一致的对子所占的比例。
C_index <- cox_more_2$concordance['concordance']
if(C_index >= 0.9){
  print("High accuracy")
}else{
  if(C_index < 0.9 & C_index >= 0.7){
    print("Medium accuracy")
  }else{
    print("Low accuracy")
  }
}
sum.surv<-summary(cox_more_2)
c_index<-sum.surv$concordance
c_index
# 06 风险模型的构建与验证------
setwd("/data/nas1/luchunlin/project/BJTC-236")
if (! dir.exists("./06_risk")){
  dir.create("./06_risk")
}
setwd("./06_risk")

## 06-1 计算每个患者的风险评分，展示生存状态分布-----
cox_more_2$coefficients
riskScore=predict(cox_more_2,type="lp",newdata=train_data)
coxGene=rownames(as.data.frame(cox_more_2$coefficients))
outCol=c("OS","OS.time",coxGene)
risk=as.vector(ifelse(riskScore>median(riskScore),0,1))
risk <- as.data.frame(c(cbind(id=rownames(cbind(train_data[,outCol],
                                                riskScore,
                                                risk)),
                              cbind(train_data[,outCol],
                                    riskScore,
                                    risk))))
table(risk$risk)
#0  1 
#379 379 
library(ggplot2)
library(ggthemes)
library(Ipaper)
median(riskScore)
# [1] 0.03281169
risk_dis <- ggplot(risk, aes(x=reorder(id, riskScore), 
                             y=riskScore, 
                             color = factor(risk, 
                                            levels = c(0, 1), 
                                            labels = c("High Risk", "Low Risk")))) +
  geom_point() +
  scale_color_manual(values = c("#A73030FF", "#0073C2FF")) + 
  scale_x_discrete(breaks = risk[order(risk$riskScore),]$id[c(1,100,200,300,400,500,600,700)],
                   labels = c(1,100,200,300,400,500,600,700),
                   expand = c(0.02,0)) +
  geom_vline(xintercept = nrow(risk[which(risk$risk==1),]) + 0.5,
             lty = 2) +
  geom_hline(yintercept = median(riskScore),
             lty =2) +
  labs(x = "Patients(increasing risk score)",
       y = "Risk Score",
       title = "Train Risk Score Distribution") + 
  theme_base() +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.justification = c(0,1),
        legend.position = c(0,1),
        #        legend.margin = margin(c(-5,4,4,3)),
        legend.background = element_rect(color = "black", size = .3),
        plot.title = element_text(size = 15, hjust = 0.5))
risk_dis
surv_stat <- ggplot(risk, aes(x=reorder(id, riskScore),
                              y=OS.time/365,
                              color = factor(OS,
                                             levels = c(0,1),
                                             labels = c("Alive", "Dead")))) +
  geom_point() +
  scale_color_manual(values = c("#0073C2FF","#A73030FF")) +
  scale_x_discrete(breaks = risk[order(risk$riskScore),]$id[c(1,100,200,300,400,500,600,700)],
                   labels = c(1,100,200,300,400,500,600,700),
                   expand = c(0.02,0)) +
  ylim(c(0,40))+
  geom_vline(xintercept = nrow(risk[which(risk$risk==1),]) + 0.5,
             lty = 2) +
  labs(x = "Patients(increasing risk score)",
       y = "Progression-free Interval (days)",
       title = "Train survival state distribution") + 
  theme_base() +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.justification = c(0,1),
        legend.position = c(0,1),
        #        legend.margin = margin(c(-5,4,4,3)),
        legend.background = element_rect(color = "black", size = .3),
        plot.title = element_text(size = 15, hjust = 0.5))

surv_stat
library(survival)
library(survminer)
kmfit<-survfit(Surv(OS.time, OS) ~ risk, data =  risk)
train_survival_median <- ggsurvplot(kmfit,
                                    pval = TRUE, 
                                    conf.int = F,
                                    legend.labs=c("High risk","Low risk" ),
                                    legend.title="Risk score",
                                    title="Train KM",
                                    font.main = c(15,"bold"),
                                    risk.table = TRUE, 
                                    risk.table.col = "strata", 
                                    linetype = "strata", 
                                    surv.median.line = "hv", 
                                    ggtheme = theme_bw(), 
                                    palette = c("#A73030FF", "#0073C2FF"))
train_survival_median
## 06-2 Kaplan-Meier生存分析，KM生存曲线图，1、3、5年绘制ROC曲线-------------
riskscore <- function(survival_cancer_df,
                      candidate_genes_for_cox,
                      cox_report){
  library("dplyr")
  risk_score_table <- survival_cancer_df[, candidate_genes_for_cox]
  for (each_sig_gene in colnames(risk_score_table)){
    risk_score_table$each_sig_gene <- risk_score_table[,each_sig_gene]*
      (summary(cox_report)$coefficients[each_sig_gene,1])
  }
  risk_score_table <- cbind(risk_score_table,
                            "total_risk_score"=exp(rowSums(risk_score_table))) %>%
    cbind(survival_cancer_df[,c("OS.time", "OS")])
  risk_score_table <- risk_score_table[,c("OS.time",
                                          "OS",
                                          candidate_genes_for_cox,
                                          "total_risk_score")]
}
candidate_genes_for_cox2 <- c(rownames(cox_table)[cox_table[,3]>0.05])
risk_score_table_multi_cox2 <- riskscore(train_data,
                                         candidate_genes_for_cox2,
                                         cox_more_2)

multi_ROC <- function(time_vector, risk_score_table){
  library(survivalROC)
  single_ROC <- function(single_time){
    for_ROC <- survivalROC(Stime=risk_score_table$OS.time,
                           status=risk_score_table$OS,
                           marker=risk_score_table$riskScore,
                           predict.time=single_time,method = 'KM')
    data.frame('True_positive'=for_ROC$TP, 'False_positive'=for_ROC$FP,
               'Cut_values'=for_ROC$cut.values, 'Time'=rep(single_time, length(for_ROC$TP)),
               'AUC'=rep(for_ROC$AUC, length(for_ROC$TP)))
  }
  multi_ROC_list <- lapply(time_vector, single_ROC)
  do.call(rbind, multi_ROC_list)
}
for_multi_ROC <- multi_ROC(time_vector = c(365*seq(1,5,2)), 
                           risk_score_table = risk)
for_multi_ROC$Time <- factor(for_multi_ROC$Time)
# 画ROC曲线 
library(scales)
library(geomROC)
library(plotROC)
auc_y1 <- round(for_multi_ROC[which(for_multi_ROC$Time==365),5][1],2)
auc_y3 <- round(for_multi_ROC[which(for_multi_ROC$Time==1095),5][1],2)
auc_y5 <- round(for_multi_ROC[which(for_multi_ROC$Time==1825),5][1],2)

ROC <- ggplot(for_multi_ROC, aes(x=False_positive,
                                 y=True_positive, 
                                 label=Cut_values, 
                                 color=Time)) + 
  scale_color_manual(breaks = c("365", "1095", "1825"),
                     labels = c("1 years", "3 years", "5 years"),
                     values = c("#00468b", "#A73030FF", "#42b540")) +
  geom_roc(labels = F, stat = 'identity') + 
  style_roc() + 
  geom_abline(slope = 1, intercept = 0, color = 'gray', linetype=2) +
  theme_bw() +
  labs(title = "Train ROC") +
  theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
        panel.grid = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold")) +
  # annotate('text', x=.75, y=.25, label=paste('AUC of 1 years =', round(auc_y1,2))) + 
  # annotate('text', x=.75, y=.15, label=paste('AUC of 3 years =', round(auc_y3,2))) + 
  # annotate('text', x=.75, y=.05, label=paste('AUC of 5 years =', round(auc_y5,2))) +
  annotate("text", x=0.75, y=c(0.25, 0.15, 0.05),
           label = c(paste('AUC of 1 years =', format(auc_y1,nsmall=2)),
                     paste('AUC of 3 years =', format(auc_y3,nsmall=2)),
                     paste('AUC of 5 years =', format(auc_y5,nsmall=2))))
ROC
ggsave(filename = 'Train_roc.pdf',ROC,w=6,h=5)
ggsave(filename = 'Train_roc.png',ROC,w=6,h=5)
## 06-3 外部数据库验证（KM，1、3、5年ROC）---------
setwd("/data/nas1/luchunlin/project/BJTC-236")
if (! dir.exists("./07_Validation")){
  dir.create("./07_Validation")
}
setwd("./07_Validation")
#ICGC_dat<-read_xlsx('/data/nas1/luchunlin/project/BJTC-236/08_External_va/ICGC_dat.xlsx')
#ICGC_dat<-as.data.frame(ICGC_dat)
#rownames(ICGC_dat)<-ICGC_dat$Symbol
#ICGC_dat<-ICGC_dat[,-1]
#ICGC_survival<-read_xlsx('/data/nas1/luchunlin/project/BJTC-236/08_External_va/ICGC_survival.xlsx')
#test_dat<-t(scale(ICGC_dat))
#test_dat<-as.data.frame(test_dat)
#test_dat$sample<-rownames(test_dat)
#test_dat<-merge(ICGC_survival,test_dat,by='sample')
#rownames(test_dat)<-test_dat$sample
#test_dat<-test_dat[,-1]

# 开始验证

riskScore_out = predict(cox_more_2,type="lp",newdata=test_data)
risk_out=as.vector(ifelse(riskScore_out>median(riskScore_out),0,1))
risk_out <- as.data.frame(c(cbind(id=rownames(cbind(test_data[,outCol],
                                                    riskScore_out,
                                                    risk_out)),
                                  cbind(test_data[,outCol],
                                        riskScore_out,
                                        risk_out))))

library(ggplot2)
library(ggthemes)
median(riskScore_out)

risk_dis_out <- ggplot(risk_out, aes(x=reorder(id, riskScore_out), 
                                     y=riskScore_out, 
                                     color = factor(risk_out, 
                                                    levels = c(0, 1), 
                                                    labels = c("High Risk", "Low Risk")))) +
  geom_point() +
  scale_color_manual(values = c("#A73030FF", "#0073C2FF")) + 
  scale_x_discrete(breaks = risk_out[order(risk_out$riskScore_out),]$id[c(1,100,200,300,400,500,600,700,800)],
                   labels = c(1,100,200,300,400,500,600,700,800),
                   expand = c(0.02,0)) +
  geom_vline(xintercept = nrow(risk_out[which(risk_out$risk_out==1),]) + 0.5,
             lty = 2) +
  geom_hline(yintercept = median(riskScore_out),
             lty =2) +
  labs(x = "Patients(increasing risk score)",
       y = "Risk Score",
       title = "Validation Risk Score Distribution") + 
  theme_base() +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.justification = c(0,1),
        legend.position = c(0,1),
        #        legend.margin = margin(c(-5,4,4,3)),
        legend.background = element_rect(color = "black", size = .3),
        plot.title = element_text(size = 15, hjust = 0.5))
risk_dis_out
surv_stat_out <- ggplot(risk_out, aes(x=reorder(id, riskScore_out),
                                      y=OS.time/365,
                                      color = factor(OS,
                                                     levels = c(0,1),
                                                     labels = c("Alive", "Dead")))) +
  geom_point() +
  scale_color_manual(values = c("#0073C2FF", "#A73030FF")) +
  scale_x_discrete(breaks = risk_out[order(risk_out$riskScore_out),]$id[c(1,100,200,300,400,500,600,700,800)],
                   labels = c(1,100,200,300,400,500,600,700,800),
                   expand = c(0.02,0)) +
  ylim(x=c(0,40)) +
  geom_vline(xintercept = nrow(risk_out[which(risk_out$risk_out==1),]) + 0.5,
             lty = 2) +
  labs(x = "Patients(increasing risk score)",
       y = "Survival time (years)",
       title = "Validation Survival State Distribution") + 
  theme_base() +
  theme(legend.title = element_blank(),
        legend.key = element_blank(),
        legend.justification = c(0,1),
        legend.position = c(0,1),
        #        legend.margin = margin(c(-5,4,4,3)),
        legend.background = element_rect(color = "black", size = .3),
        plot.title = element_text(size = 15, hjust = 0.5))
surv_stat_out

library(survival)
library(survminer)
kmfit_out <- survfit(Surv(OS.time, OS) ~ risk_out, data =  risk_out)
verify_survival_median <- ggsurvplot(kmfit_out,
                                     pval = TRUE, 
                                     conf.int = F,
                                     legend.labs=c("High risk","Low risk" ),
                                     legend.title="Risk score",
                                     title="Validition KM",
                                     font.main = c(15,"bold"),
                                     risk.table = TRUE, 
                                     risk.table.col = "strata", 
                                     linetype = "strata", 
                                     surv.median.line = "hv", 
                                     ggtheme = theme_bw(), 
                                     palette = c("#A73030FF", "#0073C2FF"))
verify_survival_median
riskscore <- function(survival_cancer_df,
                      candidate_genes_for_cox,
                      cox_report){
  library("dplyr")
  risk_score_table <- survival_cancer_df[, candidate_genes_for_cox]
  for (each_sig_gene in colnames(risk_score_table)){
    risk_score_table$each_sig_gene <- risk_score_table[,each_sig_gene]*
      (summary(cox_report)$coefficients[each_sig_gene,1])
  }
  risk_score_table <- cbind(risk_score_table,
                            "total_risk_score"=exp(rowSums(risk_score_table))) %>%
    cbind(survival_cancer_df[,c("OS.time", "OS")])
  risk_score_table <- risk_score_table[,c("OS.time",
                                          "OS",
                                          candidate_genes_for_cox,
                                          "total_risk_score")]
}
candidate_genes_for_cox2 <- c(rownames(cox_table)[cox_table[,3]>0.05])
risk_score_table_multi_cox2 <- riskscore(test_data,
                                         candidate_genes_for_cox2,
                                         cox_more_2)

multi_ROC <- function(time_vector, risk_score_table){
  library(survivalROC)   
  single_ROC <- function(single_time){
    for_ROC <- survivalROC(Stime=risk_score_table$OS.time,
                           status=risk_score_table$OS,
                           marker=risk_score_table$riskScore_out,
                           predict.time=single_time,method = 'KM')
    data.frame('True_positive'=for_ROC$TP, 'False_positive'=for_ROC$FP,
               'Cut_values'=for_ROC$cut.values, 'Time'=rep(single_time, length(for_ROC$TP)),
               'AUC'=rep(for_ROC$AUC, length(for_ROC$TP)))
  }
  multi_ROC_list <- lapply(time_vector, single_ROC)
  do.call(rbind, multi_ROC_list)
}
for_multi_ROC <- multi_ROC(time_vector = c(365*seq(1,5,2)), 
                           risk_score_table = risk_out)
for_multi_ROC$Time <- factor(for_multi_ROC$Time)

# 画ROC曲线 
library(scales)
library(geomROC)
library(plotROC)
auc_y1 <- round(for_multi_ROC[which(for_multi_ROC$Time==365),5][1],2)
auc_y3 <- round(for_multi_ROC[which(for_multi_ROC$Time==1095),5][1],2)
auc_y5 <- round(for_multi_ROC[which(for_multi_ROC$Time==1825),5][1],2)

ROC <- ggplot(for_multi_ROC, aes(x=False_positive,
                                 y=True_positive, 
                                 label=Cut_values, 
                                 color=Time)) + 
  scale_color_manual(breaks = c("365", "1095", "1825"),
                     labels = c("1 years", "3 years", "5 years"),
                     values = c("#00468b", "#A73030FF", "#42b540")) +
  geom_roc(labels = F, stat = 'identity') + 
  style_roc() + 
  geom_abline(slope = 1, intercept = 0, color = 'gray', linetype=2) +
  theme_bw() +
  labs(title = "Validation ROC") +
  theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
        panel.grid = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold")) +
  # annotate('text', x=.75, y=.25, label=paste('AUC of 1 years =', round(auc_y1,2))) + 
  # annotate('text', x=.75, y=.15, label=paste('AUC of 3 years =', round(auc_y3,2))) + 
  # annotate('text', x=.75, y=.05, label=paste('AUC of 5 years =', round(auc_y5,2))) +
  annotate("text", x=0.75, y=c(0.25, 0.15, 0.05),
           label = c(paste('AUC of 1 years =', format(auc_y1,nsmall=2)),
                     paste('AUC of 3 years =', format(auc_y3,nsmall=2)),
                     paste('AUC of 5 years =', format(auc_y5,nsmall=2))))
ROC
ggsave(filename = 'Validation_roc.pdf',ROC,w=6,h=5)
ggsave(filename = 'Validation_roc.png',ROC,w=6,h=5)
# 07 预后模型的构建与评价--------
setwd("/data/nas1/luchunlin/project/BJTC-236")
if (! dir.exists("./08_prog_model")){
  dir.create("./08_prog_model")
}
setwd("./08_prog_model")
## 07-1 单因素Cox----------
train_phenotype<-clinical
survival_cancer2<-survival_cancer[survival_cancer$sample%in%rownames(train_data),]
train_phenotype<-train_phenotype[train_phenotype$submitter_id%in%survival_cancer2$'_PATIENT',]
colnames(survival_cancer2)<-c('sample','OS','submitter_id','OS.time')
train_phenotype<-merge(survival_cancer2,train_phenotype,by='submitter_id')
train_phenotype<-data.frame(id=train_phenotype$sample,
                            age=train_phenotype$age_at_index,
                            #                           gender=train_phenotype$gender,
                            #                           stage=train_phenotype$ajcc_pathologic_stage,
                            T_stage=train_phenotype$ajcc_pathologic_t,
                            N_stage=train_phenotype$ajcc_pathologic_n,
                            M_stage=train_phenotype$ajcc_pathologic_m,
                            OS=as.numeric(train_phenotype$OS),
                            OS.time=as.numeric(train_phenotype$OS.time))
train_phenotype[train_phenotype == "N/A"] <- NA
## 将不常见分期替换为NA
#train_phenotype$stage<-gsub('Stage X',NA,train_phenotype$stage)
train_phenotype$T_stage<-gsub('TX',NA,train_phenotype$T_stage)
train_phenotype$N_stage<-gsub('NX',NA,train_phenotype$N_stage)
train_phenotype$M_stage<-gsub('MX',NA,train_phenotype$M_stage)

## 改成数值型
train_phenotype$age <- as.numeric(train_phenotype$age)
#train_phenotype$gender <- ifelse(train_phenotype$gender == "male", 1, 2)
#train_phenotype$stage<- gsub("Stage", "", train_phenotype$stage)
#train_phenotype$stage<- gsub("A", "", train_phenotype$stage)
#train_phenotype$stage<- gsub("B", "", train_phenotype$stage)
#train_phenotype$stage<- gsub("C", "", train_phenotype$stage)
#train_phenotype$stage<- gsub("IV", "4", train_phenotype$stage)
#train_phenotype$stage<- gsub("III", "3", train_phenotype$stage)
#train_phenotype$stage<- gsub("II", "2", train_phenotype$stage)
#train_phenotype$stage<- gsub("I", "1", train_phenotype$stage)
#train_phenotype$stage<-as.numeric(train_phenotype$stage)

train_phenotype$T_stage <- gsub("T", "", train_phenotype$T_stage)
train_phenotype$T_stage<-gsub("a", "", train_phenotype$T_stage)
train_phenotype$T_stage<-gsub("b", "", train_phenotype$T_stage)
train_phenotype$T_stage<-gsub("c", "", train_phenotype$T_stage)
train_phenotype$T_stage<-gsub("d", "", train_phenotype$T_stage)
train_phenotype$T_stage <- as.numeric(train_phenotype$T_stage)

train_phenotype$M_stage <- gsub("M", "", train_phenotype$M_stage)
train_phenotype$M_stage <- as.numeric(train_phenotype$M_stage)

train_phenotype$N_stage <- gsub("N", "", train_phenotype$N_stage)
train_phenotype$N_stage <- gsub("a", "", train_phenotype$N_stage)
train_phenotype$N_stage <- gsub("b", "", train_phenotype$N_stage)
train_phenotype$N_stage <- gsub("c", "", train_phenotype$N_stage)
train_phenotype$N_stage <- gsub("(i-)", "", train_phenotype$N_stage,fixed = T)
train_phenotype$N_stage <- gsub("(i+)", "", train_phenotype$N_stage,fixed = T)
train_phenotype$N_stage <- gsub("mi", "", train_phenotype$N_stage)
train_phenotype$N_stage <- as.numeric(train_phenotype$N_stage)

dim(train_phenotype)

# train_phenotype <- train_phenotype[-which(train_phenotype$OS.time == 0),]
all(rownames(train_data) %in% train_phenotype$id)
dim(train_phenotype)

sub_risk <- subset(risk, select = c(id, riskScore))
train_risk_clinical <- merge(train_phenotype,
                             sub_risk,
                             by = "id")
rownames(train_risk_clinical) <- train_risk_clinical$id
train_risk_clinical = subset(train_risk_clinical, select = -c(id))
dim(train_risk_clinical)
colnames_train <- colnames(train_risk_clinical)
covariates_train <- colnames_train[-which(colnames_train %in% c("OS", "OS.time"))]

train_risk_clinical$T_stage <- factor(train_risk_clinical$T_stage)
train_risk_clinical$N_stage <- factor(train_risk_clinical$N_stage)
train_risk_clinical$M_stage <- factor(train_risk_clinical$M_stage)
#train_risk_clinical$stage <- factor(train_risk_clinical$stage)
library(survival)
res.risk = coxph(Surv(time = OS.time, event = OS) ~ riskScore, data = train_risk_clinical) %>% summary
res.risk = c(res.risk$conf.int[-2], res.risk$coefficients[5])
res.age = coxph(Surv(time = OS.time, event = OS) ~ age, data = train_risk_clinical) %>% summary
res.age = c(res.age$conf.int[-2], res.age$coefficients[5])

#res.gender = coxph(Surv(time = OS.time, event = OS) ~ gender, data = train_risk_clinical) %>% summary
#res.gender = c(res.gender$conf.int[-2], res.gender$coefficients[5])


res.T_stage = coxph(Surv(time = OS.time, event = OS) ~ T_stage, data = train_risk_clinical) %>% summary
res.T_stage = cbind(res.T_stage$conf.int[,-2], res.T_stage$coefficients[,5])

res.M_stage = coxph(Surv(time = OS.time, event = OS) ~ M_stage, data = train_risk_clinical) %>% summary
res.M_stage = c(res.M_stage$conf.int[-2], res.M_stage$coefficients[5])

res.N_stage = coxph(Surv(time = OS.time, event = OS) ~ N_stage, data = train_risk_clinical) %>% summary
res.N_stage = cbind(res.N_stage$conf.int[,-2], res.N_stage$coefficients[,5])

#res.stage = coxph(Surv(time = OS.time, event = OS) ~ stage, data = train_risk_clinical) %>% summary
#res.stage = cbind(res.stage$conf.int[,-2], res.stage$coefficients[,5])

res.ref.T1 = c(1,1,1,NA)
res.ref.N0 = c(1,1,1,NA)
#res.ref.1 = c(1,1,1,NA)
res = rbind(res.risk, res.age, res.ref.T1, res.T_stage,res.ref.N0, res.N_stage, res.M_stage) %>% as.data.frame()
rownames(res)
res$Indicators = c("riskScore", "Age","T1(Reference)","T2","T3", "T4","N0(Reference)","N1", "N2",'N3',
                   'M1 vs M0')
colnames(res) = c("hr","low","up","pv","Indicator")
res$p = signif(res$pv, 2) %>% paste0("p = ", .)
res$p[is.na(res$pv)] = NA
res$Indicator = factor(res$Indicator, levels = rev(res$Indicator))
rownames(res) <- res$Indicator
res2 <- data.frame(p.value=res$pv,
                   HR=res$hr,
                   HR.95L=res$low,
                   HR.95H=res$up,
                   Indicator=res$Indicator)
rownames(res2) <- res2$Indicator
write.table(res2, file = "univariate_cox_prog_forest.xls", sep = "\t", quote = F)
res2 <- subset(res2, select = -c(Indicator))
library(tidyr)
hz <- paste(round(res2$HR,3),
            "(",round(res2$HR.95L,3),
            "-",round(res2$HR.95H,3),")",sep = "")
hz
hz[c(3,7)] <- ""

tabletext <- cbind(c(NA,rownames(res2)),
                   c("P value",ifelse(res2$p.value<0.0001,
                                      "< 0.0001",
                                      round(res2$p.value,4))),
                   c("Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1
library(forestplot)
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,rep(FALSE, 7)),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,res2$HR),
           lower=c(NA,res2$HR.95L), #95%置信区间下限
           upper=c(NA,res2$HR.95H), #95%置信区间上限
           boxsize=0.1,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(0, 1, 2, 3,4,5, 6,7, 8,9,10), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(1.2,"cm"), #固定行高
           graphwidth = unit(.6,"npc"), #图在表中的宽度比例
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),
           # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
           #                 "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
           #                 "59" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
           # mar=unit(rep(0.01, times = 4), "cm"),#图形页边距
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Hazard Ratio",
           grid = T,
           title = "Univariate") # 垂直于x轴的网格线，对应每个刻度
## 07-2 多因素Cox----------

res.mul = coxph(Surv(time = OS.time, event = OS) ~ riskScore + age + T_stage + N_stage + M_stage, data = train_risk_clinical)%>% summary
res.mul = cbind(res.mul$conf.int[,-2], res.mul$coefficients[,5]) %>% as.data.frame()
res.mul = rbind(res.mul[1:2,], res.ref.T1, res.mul[c(3:5),],res.ref.N0,res.mul[6:9,])
res.mul$Indicators = c("riskScore","Age","T1(Reference)", "T2","T3","T4", "N0(Reference)","N1",'N2','N3', "M1 vs M0")
colnames(res.mul) = c("hr","low","up","pv","Indicator")
res.mul$p = signif(res.mul$pv, 2) %>% paste0("p = ", .)
res.mul$p[is.na(res.mul$pv)] = NA
res.mul$Indicator = factor(res.mul$Indicator, levels = rev(res.mul$Indicator))
rownames(res.mul) <- res.mul$Indicator

multi_res <- data.frame(p.value=res.mul$pv,
                        HR=res.mul$hr,
                        HR.95L=res.mul$low,
                        HR.95H=res.mul$up,
                        Indicator=res.mul$Indicator)
rownames(multi_res) <- multi_res$Indicator
multi_res
write.csv(multi_res,
          file = "multivariate_cox_prog_result.csv",
          quote = F,
          row.names = T)
multi_res <- subset(multi_res, select = -c(Indicator))

library(tidyr)
hz <- paste(round(multi_res$HR,3),
            "(",round(multi_res$HR.95L,3),
            "-",round(multi_res$HR.95H,3),")",sep = "")
hz

hz[c(3,7)] <- ""
tabletext <- cbind(c(NA,rownames(multi_res)),
                   c("P value",ifelse(multi_res$p.value<0.0001,
                                      "< 0.0001",
                                      round(multi_res$p.value,4))),
                   c("Hazard Ratio(95% CI)",hz))
nrow(tabletext) + 1
forestplot(labeltext=tabletext, 
           graph.pos=4,  #为Pvalue箱线图所在的位置
           is.summary = c(TRUE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE,FALSE,FALSE,TRUE,FALSE,FALSE),
           col=fpColors(box="red", lines="royalblue", zero = "gray50"),
           mean=c(NA,multi_res$HR),
           lower=c(NA,multi_res$HR.95L), #95%置信区间下限
           upper=c(NA,multi_res$HR.95H), #95%置信区间上限
           boxsize=0.1,lwd.ci=3,   #箱子大小，线的宽度
           ci.vertices.height = 0.08,ci.vertices=TRUE, #置信区间用线宽、高、型
           zero=1,lwd.zero=0.5,    #zero线宽 基准线的位置
           colgap=unit(5,"mm"),    #列间隙
           xticks = c(0, 1, 2,3, 4,5, 6,7, 8), #横坐标刻度
           lwd.xaxis=2,            #X轴线宽
           lineheight = unit(1.2,"cm"), #固定行高
           graphwidth = unit(.6,"npc"), #图在表中的宽度比例
           cex=0.9, fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           hrzl_lines = list("2" = gpar(col = "black", lty = 1, lwd = 2)),
           # hrzl_lines=list("2" = gpar(lwd=2, col="black"),
           #                 "3" = gpar(lwd=2, col="black"), #第三行顶部加黑线，引号内数字标记行位置
           #                 "59" = gpar(lwd=2, col="black")),#最后一行底部加黑线,"16"中数字为nrow(tabletext)+1
           # mar=unit(rep(0.01, times = 4), "cm"),#图形页边距
           #fpTxtGp函数中的cex参数设置各个组件的大小
           txt_gp=fpTxtGp(label=gpar(cex=1),
                          ticks=gpar(cex=0.8, fontface = "bold"),
                          xlab=gpar(cex = 1, fontface = "bold"),
                          title=gpar(cex = 1.25, fontface = "bold")),
           xlab="Hazard Ratio",
           grid = T,
           title = "Multivariate") # 垂直于x轴的网格线，对应每个刻度
##07-3 构建COX模型，绘制列线图---------
multi_cov<-c('riskScore','age',"T_stage","N_stage","M_stage")
cox_data_prog <- as.formula(paste0('Surv(OS.time, OS)~',
                                   paste(multi_cov,
                                         sep = '',
                                         collapse = '+')))
cox_more_prog <- coxph(cox_data_prog,
                       data = as.data.frame(train_risk_clinical))

# Nomogram
library(rms)
ddist <- datadist(train_risk_clinical)
options(datadist='ddist')

cox_data_prog


# 构建COX模型，绘制列线图
res.cox <- psm(cox_data_prog,
               data = train_risk_clinical, dist = 'lognormal')
surv <- Survival(res.cox) # 构建生存概率函数
function(x) surv(365, x) # 1年事件发生概率
function(x) surv(1095, x) # 2年事件发生概率
function(x) surv(1825, x) # 3年事件发生概率

nom.cox <- nomogram(res.cox,
                    fun = list(function(x) surv(365, x),
                               function(x) surv(1095, x),
                               function(x) surv(1825, x)),
                    funlabel=c("1-year Survival Probability", "3-year Survival Probability", "5-year Survival Probability"),
                    maxscale = 10,
                    fun.at = c(0.01,seq(0.1,0.9,by=0.2),0.95,0.99),
                    lp=F)

plot(nom.cox, cex.axis  = 1.5, cex.var = 1.7)
png(filename = "nomogram_line_points.png", height = 800, width = 1200)
plot(nom.cox, cex.axis  = 1.5, cex.var = 1.7)
dev.off()
pdf(file = "nomogram_line_points.pdf", height = 10, width = 17)
plot(nom.cox, cex.axis  = 1.5, cex.var = 1.7)
dev.off()
##07-4 构建校准曲线---------

##绘制1年生存期校准曲线
coxm_1 <- cph(cox_data_prog,
              data=train_risk_clinical,
              surv=T,
              x=T,
              y=T,
              time.inc = 365)

cal_1<-calibrate(coxm_1,u=365,cmethod='KM',m=500,B=1000)

##绘制2年生存期校曲线
##time.in 和 u 要是一样的，都是要评价的时间节点
coxm_3 <- cph(cox_data_prog,
              data=train_risk_clinical,
              surv=T,
              x=T,y=T,
              time.inc = 3*365)
cal_3 <-calibrate(coxm_3,u=3*365,cmethod='KM',m=500,B=1000)

coxm_5 <- cph(cox_data_prog,
              data=train_risk_clinical,
              surv=T,
              x=T,y=T,
              time.inc = 5*365)
cal_5 <-calibrate(coxm_5,u=5*365,cmethod='KM',m=500,B=1000)

par(mar=c(7,4,4,3),cex=1.5)
plot(cal_1,
     subtitles = F,
     lwd=2,lty=1, ##设置线条形状和尺寸
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 1-5 year Progression-free Interval',#便签
     ylab='Actual 1-5 year Progression-free Interval(Proportion)',#标签
     col="#00468b",#设置一个颜色
     xlim = c(0.5,1),ylim = c(0.3,1)) ##x轴和y轴范围
plot(cal_3,
     add = T,
     subtitles = F,
     lwd=2,lty=1,  ##设置线条宽度和线条类型
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 1-5 year Progression-free Interval',#便签
     ylab='Actual 1-5 year Progression-free Interval(Proportion)',#标签
     col="#ed0000",#设置一个颜色
     xlim = c(0.5,1),ylim = c(0.3,1)) ##x轴和y轴范围
plot(cal_5,
     add = T,
     subtitles = F,
     lwd=2,lty=1, ##设置线条形状和尺寸
     errbar.col=c(rgb(0,118,192,maxColorValue = 255)), ##设置一个颜色
     xlab='Nomogram-Predicted Probability of 1-5year Progression-free Interval',#便签
     ylab='Actual 1-year Progression-free Interval(Proportion)',#标签
     col="#42b540",#设置一个颜色
     xlim = c(0.5,1),ylim = c(0.3,1)) ##x轴和y轴范围



#加上图例
legend("bottomright", legend=c("1-year", "3-year", "5-year"), 
       col=c("#00468b", "#ed0000", "#42b540"), 
       lwd=2)
#调整对角线
abline(0,1,lty=5,lwd=2,col="grey")

## 07-5 KM----------
res.mul = coxph(Surv(time = OS.time, event = OS)~ riskScore + age + T_stage + N_stage + M_stage, data = train_risk_clinical)
train_risk_clinical$riskScore_prog = predict(res.mul, newdata = train_risk_clinical, type = "lp")

train_risk_clinical$risk_prog = ifelse(train_risk_clinical$riskScore_prog > median(train_risk_clinical$riskScore_prog, na.rm = T), "high", "low")
train_risk_clinical$risk_prog = factor(train_risk_clinical$risk_prog, levels = c("high", "low"), labels = c("High risk", "Low risk"))
surv.fit = survfit(Surv(time = OS.time,event = OS) ~ risk_prog, data = train_risk_clinical)
prog_survival_median <- ggsurvplot(surv.fit,
                                   pval = TRUE, 
                                   conf.int = F,
                                   legend.labs=c("High risk","Low risk" ),
                                   legend.title="Risk score",
                                   title="Prognosis Model KM",
                                   font.main = c(15,"bold"),
                                   risk.table = TRUE, 
                                   risk.table.col = "strata", 
                                   linetype = "strata", 
                                   surv.median.line = "hv",
                                   ggtheme = theme_bw(), 
                                   palette = c("#A73030FF", "#0073C2FF"))
prog_survival_median


## 07-6 ROC--------
library(survivalROC)
library(tidyverse)
###计算每个患者的风险评分，展示生存状态分布

# 开始验证
train_risk_clinical2 <- train_risk_clinical

outCol2 <- c("OS", "OS.time", rownames(as.data.frame(cox_more_prog$coefficients)))

library(survival)
library(survminer)

riskscore <- function(survival_cancer_df,
                      candidate_genes_for_cox,
                      cox_report){
  library("dplyr")
  risk_score_table <- survival_cancer_df[, candidate_genes_for_cox]
  for (each_sig_gene in colnames(risk_score_table)){
    risk_score_table$each_sig_gene <- risk_score_table[,each_sig_gene]*
      (summary(cox_report)$coefficients[each_sig_gene,1])
  }
  risk_score_table <- cbind(risk_score_table,
                            "total_risk_score"=exp(rowSums(risk_score_table))) %>%
    cbind(survival_cancer_df[,c("OS.time", "OS")])
  risk_score_table <- risk_score_table[,c("OS.time",
                                          "OS",
                                          candidate_genes_for_cox,
                                          "total_risk_score")]
}


# candidate_genes_for_cox2 <- c(rownames(cox_table)[cox_table[,3]>0.05])

risk_score_table_multi_cox2<-train_risk_clinical2[,-10]
multi_ROC_out <- function(time_vector, risk_score_table){
  library(survivalROC)
  single_ROC <- function(single_time){
    for_ROC <- survivalROC(Stime=risk_score_table$OS.time,
                           status=risk_score_table$OS,
                           marker=risk_score_table$riskScore_prog,
                           predict.time=single_time,method = 'KM')
    data.frame('True_positive'=for_ROC$TP, 'False_positive'=for_ROC$FP,
               'Cut_values'=for_ROC$cut.values, 'Time'=rep(single_time, length(for_ROC$TP)),
               'AUC'=rep(for_ROC$AUC, length(for_ROC$TP)))
  }
  multi_ROC_list <- lapply(time_vector, single_ROC)
  do.call(rbind, multi_ROC_list)
}

for_multi_ROC <- multi_ROC_out(time_vector = c(365*seq(1,5,2)), 
                               risk_score_table = train_risk_clinical2)
for_multi_ROC$Time <- factor(for_multi_ROC$Time)

# 画ROC曲线 
library(scales)
library(geomROC)
library(plotROC)
auc_y1 <- round(for_multi_ROC[which(for_multi_ROC$Time==365),5][1],2)
auc_y3 <- round(for_multi_ROC[which(for_multi_ROC$Time==1095),5][1],2)
auc_y5 <- round(for_multi_ROC[which(for_multi_ROC$Time==1825),5][1],2)

ROC <- ggplot(for_multi_ROC, aes(x=False_positive,
                                 y=True_positive, 
                                 label=Cut_values, 
                                 color=Time)) + 
  scale_color_manual(breaks = c("365", "1095", "1825"),
                     labels = c("1 years", "3 years", "5 years"),
                     values = c("#00468b", "#A73030FF", "#42b540")) +
  geom_roc(labels = F, stat = 'identity') + 
  style_roc() + 
  geom_abline(slope = 1, intercept = 0, color = 'gray', linetype=2) +
  theme_bw() +
  labs(title = "ROC for Prognosis Model") +
  theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
        panel.grid = element_blank(),
        plot.title = element_text(size = 15, hjust = 0.5, face = "bold")) +
  # annotate('text', x=.75, y=.25, label=paste('AUC of 1 years =', round(auc_y1,2))) + 
  # annotate('text', x=.75, y=.15, label=paste('AUC of 2 years =', round(auc_y3,2))) + 
  # annotate('text', x=.75, y=.05, label=paste('AUC of 3 years =', round(auc_y5,2))) +
  annotate("text", x=0.75, y=c(0.25, 0.15, 0.05),
           label = c(paste('AUC of 1 years =', format(auc_y1,nsmall=2)),
                     paste('AUC of 3 years =', format(auc_y3,nsmall=2)),
                     paste('AUC of 5 years =', format(auc_y5,nsmall=2))))
ROC
ggsave(filename = 'ROC for Prognosis Model.pdf',ROC,w=6,h=5)
ggsave(filename = 'ROC for Prognosis Model',ROC,w=6,h=5)
#08 风险评分与临床指标相关性分析-------------
setwd("/data/nas1/luchunlin/project/BJTC-236")
if (! dir.exists("./09_clinical_index")){
  dir.create("./09_clinical_index")
}
setwd("./09_clinical_index")

train_phenotype2<-clinical
survival_cancer2<-survival_cancer[survival_cancer$sample%in%rownames(train_data),]
train_phenotype2<-train_phenotype2[train_phenotype2$submitter_id%in%survival_cancer2$'_PATIENT',]
colnames(survival_cancer2)<-c('sample','OS','submitter_id','OS.time')
train_phenotype2<-merge(survival_cancer2,train_phenotype2,by='submitter_id')
train_phenotype2<-data.frame(id=train_phenotype2$sample,
                             age=train_phenotype2$age_at_index,
                             #gender=train_phenotype2$gender,
                             stage=train_phenotype2$ajcc_pathologic_stage,
                             T_stage=train_phenotype2$ajcc_pathologic_t,
                             N_stage=train_phenotype2$ajcc_pathologic_n,
                             M_stage=train_phenotype2$ajcc_pathologic_m,
                             OS=as.numeric(train_phenotype2$OS),
                             OS.time=as.numeric(train_phenotype2$OS.time))
train_phenotype2[train_phenotype2 == "N/A"] <- NA
## 将不常见分期替换为NA
train_phenotype2$stage<-gsub('Stage X',NA,train_phenotype2$stage)
train_phenotype2$T_stage<-gsub('TX',NA,train_phenotype2$T_stage)
train_phenotype2$N_stage<-gsub('NX',NA,train_phenotype2$N_stage)
train_phenotype2$M_stage<-gsub('MX',NA,train_phenotype2$M_stage)

## 改成数值型
train_phenotype2$age <- as.numeric(train_phenotype2$age)
#train_phenotype2$gender <- ifelse(train_phenotype2$gender == "male", 1, 2)
train_phenotype2$stage<- gsub("A", "", train_phenotype2$stage)
train_phenotype2$stage<- gsub("B", "", train_phenotype2$stage)
train_phenotype2$stage<- gsub("C", "", train_phenotype2$stage)
train_phenotype2$stage<- gsub("IV", "4", train_phenotype2$stage)
train_phenotype2$stage<- gsub("III", "3", train_phenotype2$stage)
train_phenotype2$stage<- gsub("II", "2", train_phenotype2$stage)
train_phenotype2$stage<- gsub("I", "1", train_phenotype2$stage)

train_phenotype2$T_stage<-gsub("a", "", train_phenotype2$T_stage)
train_phenotype2$T_stage<-gsub("b", "", train_phenotype2$T_stage)
train_phenotype2$T_stage<-gsub("c", "", train_phenotype2$T_stage)
train_phenotype2$T_stage<-gsub("d", "", train_phenotype2$T_stage)

train_phenotype2$N_stage <- gsub("a", "", train_phenotype2$N_stage)
train_phenotype2$N_stage <- gsub("b", "", train_phenotype2$N_stage)
train_phenotype2$N_stage <- gsub("c", "", train_phenotype2$N_stage)
train_phenotype2$N_stage <- gsub("(i-)", "", train_phenotype2$N_stage,fixed = T)
train_phenotype2$N_stage <- gsub("(i+)", "", train_phenotype2$N_stage,fixed = T)
train_phenotype2$N_stage <- gsub("mi", "", train_phenotype2$N_stage)
train_phenotype3 <- merge(sub_risk,
                          train_phenotype2,
                          by = "id")

train_phenotype3$stage <- gsub("^Stage 3$", "Stage 3/4", train_phenotype3$stage)
train_phenotype3$stage <- gsub("^Stage 4$", "Stage 3/4", train_phenotype3$stage)
train_phenotype3$N_stage <- gsub("N0", "N0/N1", train_phenotype3$N_stage)
train_phenotype3$N_stage <- gsub("^N1$", "N0/N1", train_phenotype3$N_stage)
train_phenotype3$N_stage <- gsub("^N2$", "N2/N3", train_phenotype3$N_stage)
train_phenotype3$N_stage <- gsub("^N3$", "N2/N3", train_phenotype3$N_stage)
write.table(train_phenotype3,
            file = "clinical_risk.csv",
            row.names = T,
            sep = "\t",
            quote = F)

library(ggpubr)
library(Ipaper)
library(ggthemes)
## 08-1 stage-----
my_comparisons <- list(c("Stage 1","Stage 2"),c("Stage 2","Stage 3/4"),c("Stage 1","Stage 3/4"))
stage_data <- data.frame(riskScore = train_phenotype3$riskScore,
                         stage = factor(train_phenotype3$stage,
                                        levels = c("Stage 1","Stage 2", "Stage 3/4")))
stage_data <- na.omit(stage_data)

stage<-ggplot(stage_data,aes(x = stage, y = riskScore, fill = stage)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "riskScore",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("Stage") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              test = t.test,
              map_signif_level = T,
              y_position = c(2,3,4))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 18),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13))+
  guides(fill='none')
stage

## 08-2 T stage-----
my_comparisons <- list(c("T1","T2"),c("T2","T3"),c("T3","T4"))
t_stage_data <- data.frame(riskScore = train_phenotype3$riskScore,
                           T_stage = factor(train_phenotype3$T_stage,
                                            levels = c("T1","T2","T3","T4")))
t_stage_data <- na.omit(t_stage_data)

t_stage<-ggplot(t_stage_data,aes(x = T_stage, y = riskScore, fill = T_stage)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "riskScore",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("T Stage") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              test = t.test,
              map_signif_level = T,
              y_position = c(2,3,4,4,5))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 18),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13))+
  guides(fill='none')
t_stage
## 08-3 M stage-----
my_comparisons <- list(c("M0","M1"))
m_stage_data <- data.frame(riskScore = train_phenotype3$riskScore,
                           M_stage = factor(train_phenotype3$M_stage,
                                            levels = c("M0","M1")))
m_stage_data <- na.omit(m_stage_data)

m_stage<-ggplot(m_stage_data,aes(x = M_stage, y = riskScore, fill = M_stage)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "riskScore",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("M Stage") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              test = t.test,
              map_signif_level = T,
              y_position = c(2,3,4,5))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 18),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13))+
  guides(fill='none')
m_stage
## 08-4 N stage-----
my_comparisons <- list(c("N0/N1","N2/N3"))
n_stage_data <- data.frame(riskScore = train_phenotype3$riskScore,
                           N_stage = factor(train_phenotype3$N_stage,
                                            levels = c("N0/N1","N2/N3")))
n_stage_data <- na.omit(n_stage_data)

n_stage<-ggplot(n_stage_data,aes(x = N_stage, y = riskScore, fill = N_stage)) +
  geom_boxplot(alpha=0.7) +
  scale_y_continuous(name = "riskScore",expand = c(0.1,0.1))+
  scale_x_discrete(name = "") +
  ggtitle("N Stage") +
  theme_bw() +
  geom_signif(comparisons = my_comparisons,
              test = t.test,
              map_signif_level = T,
              y_position = c(4,5))+
  theme(plot.title = element_text(size = 16, face =  "bold"),
        text = element_text(size = 18),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 13))+
  guides(fill='none')
n_stage
library(patchwork)
all_clinical_index <-  t_stage + stage + m_stage +n_stage +
  plot_layout(ncol = 2) & 
  theme_bw() & 
  theme(legend.title = element_blank(),
        legend.position = "top",
        legend.text = element_text(face = "bold"),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid = element_blank(),
        axis.title = element_text(face = "bold"),
        axis.text = element_text(face = "bold"))
all_clinical_index
# 09 GSEA-------
setwd("/data/nas1/luchunlin/project/BJTC-236")
if (! dir.exists("./10_GSEA")){
  dir.create("./10_GSEA")
}
setwd("./10_GSEA")
risk2 <- risk
risk2$risk_label <- ifelse(risk$risk == 0, "High", "Low")
gsea_exp<-exp_tumor_fpkm[,risk$id]
all(colnames(gsea_exp) == risk$id)
dim(gsea_exp)
## 13627   758
write.table(gsea_exp,
            file = "gsea_exp.xls",
            quote = F,
            sep = "\t",
            row.names = T)
# 分组
group_risk <- risk2$risk_label %>% as.factor()
design_risk <- model.matrix(~0 + group_risk)
rownames(design_risk) <- colnames(gsea_exp)
colnames(design_risk) <- levels(group_risk)
compare_risk <- makeContrasts("High-Low", levels = design_risk)

fit <- lmFit(gsea_exp, design_risk)
fit2 <- contrasts.fit(fit ,compare_risk)
fit3 <- eBayes(fit2)
allGeneSets <- topTable(fit3, coef = 1, number = Inf)

logFCcutoff <- 0
allGeneSets$change = as.factor(
  ifelse(allGeneSets$adj.P.Val < 0.05 & abs(allGeneSets$logFC) > logFCcutoff,
         ifelse(allGeneSets$logFC > logFCcutoff,'UP','DOWN'),'NOT')
)

DEGeneSets <- subset(allGeneSets,
                     allGeneSets$adj.P.Val < 0.05 & abs(allGeneSets$logFC) > 0)
write.table(allGeneSets,
            file = "DEG_for_riskscore.xls",
            quote = F,
            sep = "\t",
            row.names = T)

dim(DEGeneSets)
summary(DEGeneSets$change)
## DOWN  NOT   UP 
## 3475    0 2657 
genelist <- allGeneSets$logFC
names(genelist) <- rownames(allGeneSets)

geneList <- sort(genelist, decreasing = T)

DEGeneSets <- DEGeneSets[order(DEGeneSets$adj.P.Val),]
dim(DEGeneSets)
## GSEA-kegg----
kegg_set<- read.gmt("c2.cp.kegg.v7.5.1.symbols.gmt")
set.seed(1)
kegg_gsea <- GSEA(geneList, TERM2GENE = kegg_set, pvalueCutoff = 0.05)
kegg_result <- kegg_gsea@result
dim(kegg_result)

gseaplot2(kegg_gsea,1:5,color = c('#7B68EE','#CD3333','#20B2AA','#FF8C00','#FF6666'),
          title = 'KEGG GSEA',
          base_size = 11,
          rel_heights = c(1.5, 0.3, 0.5))

write.table(kegg_result,file = "KEGG_GSEA.xls",sep = "\t",quote = F,row.names = F)

# 10 TIICs---------
setwd("/data/nas1/luchunlin/project/BJTC-236")
if (! dir.exists("./11_TIICs")){
  dir.create("./11_TIICs")
}
setwd("./11_TIICs")
## CIBERSORT------

design_risk <- as.data.frame(design_risk)
high_sample <- rownames(design_risk)[which(design_risk$High == 1)]
low_sample <- rownames(design_risk)[which(design_risk$Low == 1)]
length(high_sample)
# [1] 379
length(low_sample)
# [1] 379

setwd('/data/nas1/luchunlin/pipeline/CIBERSORT')
write.table(gsea_exp,
            file = "gsea_exp.txt",
            quote = F,
            sep = "\t",
            row.names = T)
{ 
  source("Cibersort.R")
  result <- CIBERSORT('/data/nas1/luchunlin/pipeline/CIBERSORT/LM22.txt',
                      'gsea_exp.txt', 
                      perm = 1000, ##Permutations for significance analysis是用来计算单个样本估算免疫浸润的p值，大多数文章会采用1000次。数值越大，运行时间越久，
                      QN = F)
  cibersort_raw <- read.table("CIBERSORT-Results.txt",
                              header = T,
                              sep = "\t",
                              row.names = 1,
                              check.names = F)
  cibersort_result <- t(cibersort_raw[,-c(23,24,25)])
}
{
  tiics_result <- cibersort_result
  pvalue = padj = log2FoldChange <- matrix(0, nrow(tiics_result), 1)
  for (i in 1:nrow(tiics_result)){
    pvalue[i, 1] = p.value = wilcox.test(tiics_result[i, high_sample],
                                         tiics_result[i, low_sample])$p.value
    log2FoldChange[i, 1] = mean(tiics_result[i, high_sample]) - 
      mean(tiics_result[i, low_sample])
  }
  padj <- p.adjust(as.vector(pvalue), "fdr", n = length(pvalue))
  rTable <- data.frame(log2FoldChange, 
                       pvalue, 
                       padj,
                       row.names = rownames(tiics_result))
  high <- signif(apply(tiics_result[rownames(rTable), high_sample], 
                       1,
                       mean), 4)
  low <- signif(apply(tiics_result[rownames(rTable), low_sample], 
                      1, 
                      mean), 4)
  rTable <- data.frame(high, 
                       low,
                       rTable[, c("padj", "pvalue", "log2FoldChange")])
  rTable$immune_cell <- rownames(rTable)
  rTable$sig <- ifelse(rTable$padj < 0.05,
                       ifelse(rTable$padj < 0.01, 
                              ifelse(rTable$padj < 0.001,
                                     ifelse(rTable$padj < 0.0001,
                                            paste(rTable$immune_cell, "****",  sep = ""),
                                            paste(rTable$immune_cell, "***", sep = "")),
                                     paste(rTable$immune_cell, "**", sep = "")),
                              paste(rTable$immune_cell, "*",  sep = "")), 
                       rTable$immune_cell)
  
  write.table(rTable,
              file = "cibersort_tiics_wilcox_test.xls",
              quote = F,
              row.names = F,
              sep = '\t')
}
diff_cibersort_Table<-rTable[which(rTable$padj<0.05),]
write.table(diff_cibersort_Table,
            file = 'diff_cibersort_Table.xls',
            quote = F,
            sep = '\t')
design_risk <- as.data.frame(design_risk)
high_sample <- rownames(design_risk)[which(design_risk$High == 1)]
low_sample <- rownames(design_risk)[which(design_risk$Low == 1)]
length(high_sample)
# [1] 379
length(low_sample)
# [1] 379

# devtools::install_github("ricardo-bion/ggradar",
#                          dependencies = TRUE)

library(ggradar)
cibersort_result2 <- cibersort_result
cibersort_result2 <- cibersort_result2[rTable$immune_cell,]
rownames(cibersort_result2) <- rTable$sig
cibersort_result_t <- as.data.frame(t(cibersort_result2))
cibersort_result_t <- cbind(group = ifelse(rownames(cibersort_result_t) %in% high_sample,
                                           "High", "Low"),
                            cibersort_result_t)     
cibersort_result_t <- aggregate(cibersort_result_t[,2:ncol(cibersort_result_t)],
                                by = list(group = cibersort_result_t$group), mean)
tiics_radar <- ggradar(cibersort_result_t,
                       values.radar = c("0", "0.12", "0.24"),
                       grid.min = 0,
                       grid.mid = 0.12,
                       grid.max = 0.24,
                       group.line.width = 0.8, 
                       group.point.size = 1.5,
                       plot.extent.x.sf = 2,
                       plot.extent.y.sf = 1.25,
                       axis.label.size = 4.5,
                       grid.label.size = 4,
                       background.circle.colour = "white",
                       gridline.mid.colour = "grey",
                       gridline.min.colour = "grey",
                       axis.line.colour = "grey",
                       group.colours = c("#A73030FF", "#0073C2FF"),
                       plot.title = "TIICs Distribution") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = c(0.8,1))
tiics_radar
# 11 免疫评分、基质评分----------
setwd("/data/nas1/luchunlin/project/BJTC-236")
if (! dir.exists("./12_estimate")){
  dir.create("./12_estimate")
}
setwd("./12_estimate")

## 11-1 estimate-------
library(estimate)
expr_train <- log2(gsea_exp+1)
write.table(expr_train, 
            'expr_sample405_log2.txt', 
            col.names = T, 
            row.names = T, 
            quote = F, sep="\t")
# 生成expr_train.gct
filterCommonGenes(input.f = './expr_sample405_log2.txt', 
                  output.f = 'expr_train.gct', 
                  id = 'GeneSymbol')
# [1] "Merged dataset includes 8810 genes (1602 mismatched)."
# 生成train_purity.gct
estimateScore('expr_train.gct', 'train_purity.gct', platform="affymetrix")
# [1]"1 gene set: StromalSignature  overlap= 136"
# [1] "1 gene set: StromalSignature  overlap= 136"
es_score <- read.table('train_purity.gct', skip = 2, header = T)
immu_score <- es_score[,3:length(es_score)]
rownames(immu_score) <- es_score$NAME
immu_score<-t(immu_score)
immu_score<-as.data.frame(immu_score)
immu_score$sample<-colnames(expr_train)
rownames(immu_score)<-immu_score$sample
write.table(es_score,
            file = "es_score.xls",
            sep = "\t",
            quote = F,
            row.names = F)
violin_dat <- immu_score
violin_dat$group <- ifelse(violin_dat$sample %in% high_sample,
                           "High", "Low")
library(tidyr)
library(ggplot2)
library(ggpubr)
library(Ipaper)

{
  p1 <- ggplot(violin_dat, aes(x=group,y=StromalScore, fill=group))+
    geom_violin() + #绘制小提琴图
    stat_boxplot(geom="errorbar",
                 width=0.1,
                 position = position_dodge(0.9)) +
    geom_boxplot(width=0.5,
                 position=position_dodge(0.9),
                 outlier.shape = NA,
                 fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
    geom_point(aes(fill = group),
               size = 0.5,
               position = position_dodge(0.9))+
    scale_fill_manual(values = c("#CD3333", "#4682B4"), name = "Group") + 
    stat_compare_means(aes(group = group),method = 'wilcox')+ 
    theme_bw()+ #背景变为白色
    ylim(-3000,3000) +
    theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
          panel.grid = element_blank()) + 
    theme(axis.text.x=element_text(family="Times",size=12,face="bold"), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
          axis.text.y=element_text(family="Times",size=12,face="bold"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
          axis.title.y=element_text(family="Times",size = 15,face="bold"), #设置y轴标题的字体属性
          axis.title.x=element_text(family="Times",size = 10,face="bold"),
          plot.title = element_text(hjust = 0.5,family="Times",size = 15,face="bold"),
          legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold"))+
    labs(title="Stromal Score", x="", y="Stromal Score")
  p1
  p2 <- ggplot(violin_dat, aes(x=group,y=ImmuneScore, fill=group))+
    geom_violin() + #绘制小提琴图
    stat_boxplot(geom="errorbar",
                 width=0.1,
                 position = position_dodge(0.9)) +
    geom_boxplot(width=0.5,
                 position=position_dodge(0.9),
                 outlier.shape = NA,
                 fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
    geom_point(aes(fill = group),
               size = 0.5,
               position = position_dodge(0.9))+
    scale_fill_manual(values = c("#CD3333", "#4682B4"), name = "Group") +
    stat_compare_means(aes(group = group),method = 'wilcox')+ 
    theme_bw()+ #背景变为白色
    ylim(-2000,4000) +
    theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
          panel.grid = element_blank()) + 
    theme(axis.text.x=element_text(family="Times",size=12,face="bold"), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
          axis.text.y=element_text(family="Times",size=12,face="bold"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
          axis.title.y=element_text(family="Times",size = 15,face="bold"), #设置y轴标题的字体属性
          axis.title.x=element_text(family="Times",size = 10,face="bold"),
          plot.title = element_text(hjust = 0.5,family="Times",size = 15,face="bold"),     
          legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold"))+
    labs(title="Immune Score", x="", y="Immune Score")
  p2
  p3 <- ggplot(violin_dat, aes(x=group, y=ESTIMATEScore, fill=group))+
    geom_violin() + #绘制小提琴图
    stat_boxplot(geom="errorbar",
                 width=0.1,
                 position = position_dodge(0.9)) +
    geom_boxplot(width=0.5,
                 position=position_dodge(0.9),
                 outlier.shape = NA,
                 fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
    geom_point(aes(fill = group),
               size = 0.5,
               position = position_dodge(0.9))+
    scale_fill_manual(values = c("#CD3333", "#4682B4"), name = "Group") +
    stat_compare_means(aes(group = group),method = 'wilcox')+ 
    theme_bw()+ #背景变为白色
    theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
          panel.grid = element_blank()) + 
    ylim(-5000, 7000) +
    theme(axis.text.x=element_text(family="Times",size=12,face="bold"), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
          axis.text.y=element_text(family="Times",size=12,face="bold"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
          axis.title.y=element_text(family="Times",size = 15,face="bold"), #设置y轴标题的字体属性
          axis.title.x=element_text(family="Times",size = 10,face="bold"),
          plot.title = element_text(hjust = 0.5,family="Times",size = 15,face="bold"),
          legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold"))+
    labs(title="ESTIMATE Score", x="", y="ESTIMATE Score")
  p3
  p4 <- ggplot(violin_dat, aes(x=group, y=TumorPurity, fill=group))+
    geom_violin() + #绘制小提琴图
    stat_boxplot(geom="errorbar",
                 width=0.1,
                 position = position_dodge(0.9)) +
    geom_boxplot(width=0.5,
                 position=position_dodge(0.9),
                 outlier.shape = NA,
                 fill='white')+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
    geom_point(aes(fill = group),
               size = 0.5,
               position = position_dodge(0.9))+
    scale_fill_manual(values = c("#CD3333", "#4682B4"), name = "Group") +
    stat_compare_means(aes(group = group),method = 'wilcox')+ 
    theme_bw()+ #背景变为白色
    theme(panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
          panel.grid = element_blank()) + 
    ylim(0, 1.5) +
    theme(axis.text.x=element_text(family="Times",size=12,face="bold"), #设置x轴刻度标签的字体显示倾斜角度为15度，并向下调整1(hjust = 1)，字体簇为Times大小为20
          axis.text.y=element_text(family="Times",size=12,face="bold"), #设置y轴刻度标签的字体簇，字体大小，字体样式为plain
          axis.title.y=element_text(family="Times",size = 15,face="bold"), #设置y轴标题的字体属性
          axis.title.x=element_text(family="Times",size = 10,face="bold"),
          plot.title = element_text(hjust = 0.5,family="Times",size = 15,face="bold"),
          legend.text = element_text(face = "bold"),
          legend.title = element_text(face = "bold"))+
    labs(title="Tumor Purity", x="", y="Tumor Purity")
  p4
  p5 <- cowplot::plot_grid(p1,p2,p3,p4,
                           nrow = 2, 
                           align = 'h', 
                           vjust = -0.3)
  p5
}
# 12 免疫检验点分子-----
setwd("/data/nas1/luchunlin/project/BJTC-236")
if (! dir.exists("./13_checkpoint")){
  dir.create("./13_checkpoint")
}
setwd("./13_checkpoint")
checkpoint <- read.table("checkpoint.txt",
                         header = F)
checkpoint <- checkpoint$V1
length(checkpoint)
checkpoint_DEG <- allGeneSets[which(rownames(allGeneSets)%in%checkpoint),]
dim(checkpoint_DEG)
# 34 7
logFC_cutoff <- 0
checkpoint_DEG$change = as.factor(
  ifelse(checkpoint_DEG$adj.P.Val < 0.05 & abs(checkpoint_DEG$logFC) > logFC_cutoff,
         ifelse(checkpoint_DEG$logFC > logFC_cutoff ,'UP','DOWN'),'NOT')
)

sig_checkpoint <- rownames(checkpoint_DEG[which(checkpoint_DEG$adj.P.Val < 0.05),])
length(sig_checkpoint)
#25
sig_checkpoint

write.table(sig_checkpoint, 
            'sig_checkpoint.xls',
            quote = F, sep="\t")
sig_expr <- as.data.frame(gsea_exp[sig_checkpoint,])
sig_expr <- log2(sig_expr + 1)
sig_expr$gene <- rownames(sig_expr)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(Ipaper)
head(sig_expr[,1:3])
violin_dat <- gather(sig_expr, key=sample, value='log2(expr+1)', -c("gene"))
head(violin_dat)
violin_dat$group <- ifelse(violin_dat$sample %in% high_sample,
                           "High", "Low") 
head(violin_dat)

violin_plot <- ggplot(violin_dat, aes(x=gene, 
                                      y=`log2(expr+1)`,
                                      fill=group)) +
  #  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#A73030FF", "#0073C2FF"), name = "Group")+
  labs(title="Immune Checkpoint", x="", y = "log2(expr+1)",size=20) +
  stat_compare_means(data = violin_dat,
                     mapping = aes(group = group),
                     label ="p.signif",
                     method = 't.test',
                     paired = TRUE) +
  #  geom_signif(comparisons = my_comparisons,
  #              test = t.test,
  #              map_signif_level = T)+
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=8),
        legend.title = element_text(face = "bold", size = 10),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
violin_plot

# 13 肿瘤突变负荷----------
setwd("/data/nas1/luchunlin/project/BJTC-236")
if (! dir.exists("./14_TMB")){
  dir.create("./14_TMB")
}
setwd("./14_TMB")
library(TCGAbiolinks) 
mut_BRCA <- GDCquery_Maf("BRCA", save.csv = FALSE, directory = "GDCdata", pipelines = "varscan2")
dim(mut_BRCA)
## 93612   120
colnames(mut_BRCA) ## 查看所需要的列数
mut_BRCA<-mut_BRCA[,c(1,5,9,16)]
## 去掉正常样本，保留肿瘤样本
mut_BRCA2<-mut_BRCA[,substr(mut_BRCA[,4],14,15)!='11']
## 查看突变类型
table(mut_BRCA2$Variant_Classification)
## 去掉GenVisR包不识别的"Splice_Region"
mut_BRCA3 <- mut_BRCA2[mut_BRCA2$Variant_Classification!='Splice_Region',] 

library(GenVisR)
## 选取突变数量最高的
pdf(file='BRCA_waterfall.pdf',height=14,width=18)
waterfall(mut_BRCA3,mainRecurCutoff = 0.05, ## 筛选基因突变频率大于0.05的基因进行绘制
          mainGrid = F,     ## 是否绘制网格
          section_heights=c(3,10,3),plot_proportions=TRUE,   ## 绘制突变类型比例柱状图
          main_geneLabSize = 14)   
dev.off()
png(file='BRCA_waterfall.png',height=1400,width=1800)
waterfall(mut_BRCA3,mainRecurCutoff = 0.05, ## 筛选基因突变频率大于0.05的基因进行绘制
          mainGrid = F,     ## 是否绘制网格
          section_heights=c(3,10,3),plot_proportions=TRUE,   ## 绘制突变类型比例柱状图
          main_geneLabSize = 14)   
dev.off()

