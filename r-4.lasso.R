rm(list = ls())
setwd("/data/nas1/luchunlin/project/HF-0106-2/")
if (! dir.exists("./04_lasso")){
  dir.create("./04_lasso")
}
setwd("./04_lasso")

library(magrittr)
library(ggplot2)
hub_gene <- read.delim2('../02_intersect/DECFRGs.xls')
train.dat<-read.csv("../00_rawdata/dat(GSE169568).xls", sep='\t')  #%>% lc.tableToNum
group = read.delim2("../00_rawdata/group(GSE169568).xls")
dat<-train.dat[hub_gene$symbol,group$sample]%>%t%>%as.data.frame()
dat<-merge(dat,group,by.x='row.names',by.y='sample')%>%tibble::column_to_rownames(var = 'Row.names')
table(group$group)
dat$group<-factor(dat$group,levels = c('UC','control'))
#LASSO----------
set.seed(8)
##4  8
library(glmnet)
res.lasso <- cv.glmnet(as.matrix(dat[-ncol(dat)]), dat$group, family = "binomial", 
                       type.measure = "auc")

plot(res.lasso)
plot(res.lasso$glmnet.fit, xvar = 'lambda')
ggsave("01.lasso.CV.png", plot(res.lasso), width = 6, height = 5, dpi = 300, units = "in", bg = "white")
ggsave("01.lasso.CV.pdf", plot(res.lasso), width = 6, height = 5, dpi = 300, units = "in", bg = "white")
ggsave("02.lasso.Coef.png", plot(res.lasso$glmnet.fit, xvar = 'lambda'), width = 6, height = 5, dpi = 300, units = "in", bg = "white")
ggsave("02.lasso.Coef.pdf", plot(res.lasso$glmnet.fit, xvar = 'lambda'), width = 6, height = 5, dpi = 300, units = "in", bg = "white")
l.coef<-coef(res.lasso$glmnet.fit,s=res.lasso$lambda.min,exact= F)
l.coef
coef.min = coef(res.lasso, s = "lambda.min")  ## lambda.min & lambda.1se 取一个
res.lasso$lambda.min
# 找出那些回归系数没有被惩罚为0的
active.min = which(coef.min@i != 0)
# 提取基因名称
lasso_geneids <- coef.min@Dimnames[[1]][coef.min@i+1]
lasso_geneids <- lasso_geneids[-1]
lasso_geneids
#  [1] "MAP2K1" "CREBBP" "TAF1"   "HP"     "CXCL10" "MMP9"   "ITGA2B"  
write.table(lasso_geneids,file = 'lasso.gene.xls',sep = '\t',row.names = F,quote = F)
