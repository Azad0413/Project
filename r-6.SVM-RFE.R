rm(list = ls())
setwd("/data/nas1/luchunlin/project/HF-0106-2/")
if (! dir.exists("./06_SVM_RFE")){
  dir.create("./06_SVM_RFE")
}
setwd("./06_SVM_RFE")


library(magrittr)
library(tibble)
# 加载库和包
library(mlbench)
library(caret)
filtedat <- read.delim2('../00_rawdata/dat(GSE169568).xls',row.names = 1)%>%lc.tableToNum()
dat_group <- read.delim2('../00_rawdata/group(GSE169568).xls')%>%column_to_rownames(var = 'sample')
hubgene <- read.delim2('../02_intersect/DECFRGs.xls')

hubgene_exp <- filtedat[match(hubgene$symbol,rownames(filtedat)),]%>% 
  t %>% as.data.frame %>% rownames_to_column(var = 'sample')

hubgene_exp$group <- dat_group[match(hubgene_exp$sample,rownames(dat_group)),1]
hubgene_exp <- hubgene_exp[,-1]

# SVM-RFE-------------------------------------------------------------------------
set.seed(8)
control <- rfeControl(functions = caretFuncs, method = "cv", number = 5)
group <- as.factor(hubgene_exp$group)
# 执行SVM-RFE算法
num <- ncol(hubgene_exp)-1

results <- rfe(hubgene_exp[,1:num], 
               group, 
               sizes = c(1:num), 
               rfeControl = control,
               method = "svmRadial"
)
# saveRDS(results,"../00_rawdat/00_temp/GSE4_SVM_RFE.rds")
# 结果分析
print(results)
svmrfe_result <- predictors(results)
svmrfe_result

write.table(svmrfe_result,"svmrfe_result.txt",sep = "\t")
# saveRDS(svmrfe_result,"../00_rawdat/00_temp/svmrfe_GSE1result.rds")
# 绘制结果
pdf(file = paste0("01.SVM_RFE_Accuracy.pdf"),width = 5,height = 4)
a <- dev.cur()   #记录pdf设备
png(file = paste0("01.SVM_RFE_Accuracy.png"),width = 5, height=4, units="in", res=300) 
dev.control("enable")
par(mar = c(2,2,2,2));#下、左、上、右
plot(results, type=c("o"),
     xgap.axis = 1
)

dev.copy(which = a)  #复制来自png设备的图片到pdf
dev.off()
dev.off()
