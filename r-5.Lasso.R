rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-370-8/")
if (! dir.exists("./05_lasso")){
  dir.create("./05_lasso")
}
setwd("./05_lasso")
library(magrittr)
library(ggplot2)
hub_gene <- read.delim2('../03_DEGs/intersect.xls')
train.dat<-read.delim2("../03_DEGs/dat_final.xls", row.names = 1)  %>% lc.tableToNum
#train.dat <- log2(train.dat+1)
group = read.delim2("../03_DEGs/group.xls")
table(group$group)
dat<-train.dat[hub_gene$.,group$sample]%>%t%>%as.data.frame()
dat<-merge(dat,group,by.x='row.names',by.y='sample')%>%tibble::column_to_rownames(var = 'Row.names')
dat$group<-factor(dat$group,levels = c('RB','control'))
set.seed(123)
library(glmnet)
res.lasso <- cv.glmnet(as.matrix(dat[-ncol(dat)]), dat$group, family = "binomial", 
                       type.measure = "auc",nfolds = 10)

plot(res.lasso)
plot(res.lasso$glmnet.fit, xvar = 'lambda')
ggsave("01.lasso.CV.png", plot(res.lasso), width = 6, height = 5, dpi = 300, units = "in", bg = "white")
ggsave("01.lasso.CV.pdf", plot(res.lasso), width = 6, height = 5, dpi = 300, units = "in", bg = "white")
ggsave("02.lasso.Coef.png", plot(res.lasso$glmnet.fit, xvar = 'lambda'), width = 6, height = 5, dpi = 300, units = "in", bg = "white")
ggsave("02.lasso.Coef.pdf", plot(res.lasso$glmnet.fit, xvar = 'lambda'), width = 6, height = 5, dpi = 300, units = "in", bg = "white")
l.coef<-coef(res.lasso$glmnet.fit,s=res.lasso$lambda.min,exact= F)
l.coef
# 提取基因名称
lasso_geneids <- l.coef@Dimnames[[1]][l.coef@i+1]
lasso_geneids <- lasso_geneids[-1]
lasso_geneids
write(lasso_geneids, "lasso_genes.csv")
