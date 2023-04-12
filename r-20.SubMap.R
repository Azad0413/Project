rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-300-8/")
if (! dir.exists("./20_SubMap")){
  dir.create("./20_SubMap")
}
setwd("./20_SubMap")

##############################################################submap分析
######一般只用训练集
train <- read.delim2('../08_risk/risk.xls')
order <- train[order(train$risk,decreasing = T),]

generateInputFileForSubMap <- function(in_gct, gct_file, cls_file, sam_info, type_name = "type"){
  in_gct <- data.frame(GeneID=rownames(in_gct),
                       description="na",
                       in_gct, 
                       stringsAsFactors = F,
                       check.names = F)
  cat("#1.2\n", file = gct_file)
  cat(nrow(in_gct),"\t",ncol(in_gct)-2,"\n", file = gct_file, append = T)
  cat(paste(colnames(in_gct), collapse = "\t"),"\n", file = gct_file, append = T)
  for(i in 1:nrow(in_gct)) cat(paste(in_gct[i,], collapse = "\t"),"\n", file = gct_file, append = T)
  
  cat(nrow(sam_info),length(levels(factor(sam_info$rank))),1, "\n", file = cls_file )
  cat("#", paste0(levels(factor(sam_info[, type_name])), collapse = " " ), "\n", file = cls_file, sep = "", append = T)
  cat(as.numeric(factor(sam_info[, type_name])), file = cls_file, append = T)
}



# 创建submap需要的数据格式 (SKCM)
skcm.immunotherapy.logNC <- read.table("skcm.immunotherapy.47samples.log2CountsNorm.txt",
                                       sep = "\t",row.names = 1,header = T,
                                       check.names = F,stringsAsFactors = F) #原文提供的log2转化的标准化count值
rownames(skcm.immunotherapy.logNC) <- rownames(skcm.immunotherapy.logNC)
skcm.immunotherapy.info <- read.table("skcm.immunotherapy.47sampleInfo.txt",
                                      sep = "\t",row.names = 1,header = T,
                                      check.names = F,stringsAsFactors = F)

skcm.immunotherapy.info <- skcm.immunotherapy.info[order(skcm.immunotherapy.info$label),]
skcm.immunotherapy.info$rank <- rep(c(1,2,3,4),times=as.character(table(skcm.immunotherapy.info$label))) #1: CTLA4_noR 2: CTLA4_R 3:PD1_noR 4:PD1_R

dat<-read.delim2('../00_rawdata/dat.fpkm.xls',row.names = 1)%>%lc.tableToNum()
colnames(dat)<-gsub('.','-',colnames(dat),fixed = T)

# 创建submap需要的数据格式 (TCGA)
tmp <- dat[,train$id] # submap不允许出现flat value, 因此最好选取过滤掉低表达的表达谱，这里使用的数据过滤了超过90%样本表达值均<1的基因
# tmp <- tmp[rowSums(tmp>1)> ncol(tmp)*0.9,]#可以先不过滤基因
GENELIST <- intersect(rownames(tmp),rownames(skcm.immunotherapy.logNC)) # 取交集后的基因列表


sam_info <- skcm.immunotherapy.info
in_gct <- skcm.immunotherapy.logNC[GENELIST,rownames(skcm.immunotherapy.info)]

# 产生参考组的输出数据的文件名
gct_file <- "skcm.immunotherapy.for.SubMap.gct"
cls_file <- "skcm.immunotherapy.for.SubMap.cls"
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")



#提取出高低风险组的样本，顺序排列
ann <- order
samples.C1 <- ann[which(ann$risk == 1),'id']
samples.C2 <- ann[which(ann$risk == 0),'id']

sam_info <- data.frame("risk"=c(samples.C1,samples.C2),row.names = c(samples.C1,samples.C2))
sam_info$rank <- rep(c(1,2),times=c(length(samples.C1),length(samples.C2))) #1: C1,即HPV16-IMM 2: C2,即HPV16-KRT

# 产生输出数据的文件名
gct_file <- "Immune2.for.SubMap.gct"
cls_file <- "Immune2.for.SubMap.cls"

in_gct <- log2(tmp[GENELIST,rownames(sam_info)] + 1) # 产生和示例数据类似的形式，log2转化的标准化count值
generateInputFileForSubMap(in_gct = in_gct, gct_file = gct_file, cls_file = cls_file, sam_info = sam_info, type_name = "rank")

#--------------------------------------#
# 请参照文件夹submap使用教程完成该部分 #
#--------------------------------------#

# 参照文件夹中submap使用教程得到输出文件SubMap_SubMapResult.txt




#结果可视化
# 把submap结果/130452/SubMap_SubMapResult.txt文件中的值填入相应的位置
# 输入文件中的名义p值($nominal.p.matrix.Fisher)和校正p值(p-value correction method (for Fisher's statistics):  Bonferroni)绘制热图
heatmap.YlGnPe <- c("#440259","#345F8C","#228C8A","#78CE51","#FAE71F")
cherry    <- "#700353"
lightgrey <- "#dcddde"
library(pheatmap)
tmp <- matrix(c(1,0.5194805,1,0.01598402,1,1,1,1,
                0.7672328,0.06493506,0.9940060,0.001998002,0.7662338,0.40359640,0.1688312,0.844155844), # Bonferroni校正p值
              nrow = 4,byrow = T,dimnames = list(c("low_p","high_p","low_b","high_b"),
                                                 c("CTLA4-noR","CTLA4-R","PD1-noR","PD1-R")))

pheatmap(tmp, cellwidth = 40, cellheight = 40,
         cluster_rows = F,cluster_cols = F,
         color = heatmap.YlGnPe[5:1],
         display_numbers = matrix(ifelse(tmp < 0.05,'p < 0.05',''),nrow(tmp)),number_format = "%.3f",
         gaps_row = 2,
         annotation_row = data.frame(pvalue=c("Nominal p value","Nominal p value","Bonferroni corrected","Bonferroni corrected"),row.names = rownames(tmp)),
         annotation_colors = list(pvalue=c("Nominal p value"=lightgrey,"Bonferroni corrected"=cherry)),
         filename = "heatmap_submap.pdf")

pheatmap(tmp, cellwidth = 40, cellheight = 40,
         cluster_rows = F,cluster_cols = F,
         color = heatmap.YlGnPe[5:1],
         display_numbers = matrix(ifelse(tmp < 0.05,'p < 0.05',''),nrow(tmp)),number_format = "%.3f",
         gaps_row = 2,
         annotation_row = data.frame(pvalue=c("Nominal p value","Nominal p value","Bonferroni corrected","Bonferroni corrected"),row.names = rownames(tmp)),
         annotation_colors = list(pvalue=c("Nominal p value"=lightgrey,"Bonferroni corrected"=cherry)),
         filename = "heatmap_submap.png")
