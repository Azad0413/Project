rm(list = ls())
setwd("/data/nas1/luchunlin/project/CD-0601-2/")
if (! dir.exists("./05_lasso")){
  dir.create("./05_lasso")
}
setwd("./05_lasso")

library(magrittr)
library(ggplot2)
hub_gene <- read.delim2('../03_intersect/intersect.xls')
train.dat<-read.delim2("../00_rawdata/dat(GSE65682).xls", row.names = 1)  %>% lc.tableToNum
group = read.delim2("../00_rawdata/group(GSE65682).xls")
dat<-train.dat[hub_gene$symbol,group$sample]%>%t%>%as.data.frame()
dat<-merge(dat,group,by.x='row.names',by.y='sample')%>%tibble::column_to_rownames(var = 'Row.names')
table(group$group)
dat$group<-factor(dat$group,levels = c('Sepsis','control'))
#LASSO----------
set.seed(3)
##4
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
#  [1] "BCKDHB"   "SLC7A6"   "SLC22A4"  "SCLY"     "MCCC2"    "KCTD7"    "WDR59"    "MRI1"     "QDPR"     "QRSL1"   
# [11] "LRRC47"   "ALDH18A1" "PYCR2"    "BLMH"   
write.table(lasso_geneids,file = 'lasso.gene.xls',sep = '\t',row.names = F,quote = F)
