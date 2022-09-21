rm(list = ls())
### 验证集----------
setwd("/data/nas1/luchunlin/project/BJTC-258")
if (! dir.exists("./06_External_va")){
  dir.create("./06_External_va")
}
setwd("./06_External_va")
# 开始验证
library(GEOquery)
library(Biobase)
gset<-getGEO("GSE32062",
             destdir = '.',
             GSEMatrix = T,
             getGPL = F)
expr<-as.data.frame(exprs(gset[[2]]))

gpl<-getGEO("GPL6480",destdir = '.')
a=gset[[2]]
gpl<-Table(gpl)    
colnames(gpl)
probe2symobl<-gpl %>%
  dplyr::select('ID','GENE_SYMBOL')%>%
  filter('GENE_SYMBOL'!='')%>%
  separate('GENE_SYMBOL',c('symbol','drop'),sep = '///')%>%
  dplyr::select(-drop)
probe2symobl=probe2symobl[probe2symobl$symbol!='',]
dat<-expr
dat$ID<-rownames(dat)
dat$ID<-as.character(dat$ID)
probe2symobl$ID<-as.character(probe2symobl$ID)
dat<-dat %>%
  inner_join(probe2symobl,by='ID')%>% 
  dplyr::select(-ID)%>%     ## 去除多余信息
  dplyr::select(symbol,everything())%>%     ## 重新排列
  mutate(rowMean=rowMeans(.[grep('GSM',names(.))]))%>%    ## 求出平均数
  arrange(desc(rowMean))%>%       ## 把表达量的平均值从大到小排序
  distinct(symbol,.keep_all = T)%>%      ## symbol留下第一个
  dplyr::select(-rowMean)%>%     ## 反向选择去除rowMean这一列
  tibble::column_to_rownames(colnames(.)[1])   ## 把第一列变成行名并删除
pd<-pData(a)
setwd("/data/nas1/luchunlin/project/BJTC-258/00_rawdata/")
write.table(dat,file = 'dat.va.xls',sep = '\t',row.names = T,quote = F)
setwd("/data/nas1/luchunlin/project/BJTC-258/06_External_va/")
# pca-----
uncoxgene<-read.delim2('/data/nas1/luchunlin/project/BJTC-258/04_PCA/univariate_cox_result_0.05.xls')
uncoxgene<-data.frame(symbol=rownames(uncoxgene))
#dat.pca<-log2(dat+1)
dat.pca<-dat
dat.pca <- dat.pca[uncoxgene$symbol,]
dat.pca <- na.omit(dat.pca)
#mads <- apply(dat.pca,1,mad)
#dat.pca <- dat.pca[rev(order(mads))[1:800],]
dat.pca <- scale(dat.pca,center = T,scale = T)
dat.pca <- as.data.frame(t(dat.pca))
pca<-prcomp(dat.pca
            #,center = T,scale. = T
)
screeplot(pca,type="lines")
df<-pca$x  ## 提取PC score
df <- as.data.frame(df)
# 提取主成分的方差贡献率,生成坐标轴标题
summ <- summary(pca)
summ
xlab <- paste0("PC1(",round(summ$importance[2,1]*100,2),"%)")
ylab <- paste0("PC2(",round(summ$importance[2,2]*100,2),"%)")
#pcaScore<-data.frame(sample=rownames(df),scores=df$PC1)
## 生存信息
y <- eigen(cor(dat.pca))
# 根据累计贡献率大于60%确定主成分
sum(y$values[1:10])/sum(y$values)
s <- df[,1:10]
scores = 0.0
for (i in 1:10)
  scores=(y$values[i]*s[,i])/(sum(y$values[1:10]))+scores
pcaScore <- cbind(s,scores) # 输出综合得分信息
pcaScore <- data.frame(sample=rownames(pcaScore),scores=pcaScore$scores)
survival<-data.frame(sample=pd$geo_accession,OS.time=pd$`os (m):ch1`,OS=pd$`death (1):ch1`)
table(survival$OS)
#survival$OS<-ifelse(survival$OS=='Dead',1,0)
#dat<-dat[,survival$sample]
dat <- merge(pcaScore,survival,by = "sample")
##计算最佳截断值
dat$scores<-as.numeric(dat$scores)
dat$OS.time <- as.numeric(dat$OS.time)*30
dat$OS <- as.numeric(dat$OS)
library(survminer)
# res.cut<-surv_cutpoint(dat,time = 'OS.time',event = 'OS',variables = 'scores')
# summary(res.cut)
# cutpoint<-'1.478678'
dat$group<-ifelse(dat$scores>median(dat$scores),'High','Low')
table(dat$group)
# 去除掉"score"列
dat <- dat[,-2]
group<-dat[,-c(2,3)]
write.table(group,file = 'group.va.xls',sep = '\t',row.names = F,quote = F)
dat$group <- factor(dat$group)
# 进行KM生存分析
library(survival)
kmfit <- survfit(Surv(OS.time,OS) ~ group, data = dat)
# 绘制KM曲线
KM <- ggsurvplot(kmfit,
                 pval = TRUE, 
                 conf.int = F,
                 legend.labs=c("High","Low" ),
                 legend.title="group",
                 title="KM",
                 font.main = c(15,"bold"),
                 risk.table = TRUE, 
                 risk.table.col = "strata", 
                 linetype = "strata", 
                 surv.median.line = "hv", 
                 ggtheme = theme_bw(), 
                 palette = c("#A73030FF", "#0073C2FF"))
KM

dat <- merge(pcaScore,survival,by = "sample")
##计算最佳截断值
dat$scores<-as.numeric(dat$scores)

dat$OS.time <- as.numeric(dat$OS.time)
dat$OS <- as.numeric(dat$OS)
# res.cut<-surv_cutpoint(dat,time = 'OS.time',event = 'OS',variables = 'scores')
# summary(res.cut)
# cutpoint<-'-1.848056'
dat$group<-ifelse(dat$scores>median(dat$scores),'High','Low')
library(survivalROC)
library(tidyverse)
rt = subset(dat, select = c(OS, OS.time, scores))
rt$OS.time <- rt$OS.time / 365
survivalROC_helper <- function(t) {
  
  survivalROC(Stime=rt$OS.time, status=rt$OS, marker = rt$scores, 
              
              predict.time =t, method="KM")
}
survivalROC_data <- data_frame(t = c(5,7,9)) %>%
  
  mutate(survivalROC = map(t, survivalROC_helper),
         
         ## Extract scalar AUC
         
         auc = map_dbl(survivalROC, magrittr::extract2, "AUC"),
         
         ## Put cut off dependent values in a data_frame
         
         df_survivalROC = map(survivalROC, function(obj) {
           
           as_data_frame(obj[c("cut.values","TP","FP")])
           
         })) %>%
  
  dplyr::select(-survivalROC) %>%
  
  unnest() %>%
  
  arrange(t, FP, TP)
survivalROC_data1 <- survivalROC_data %>% 
  
  mutate(auc =sprintf("%.2f",auc))%>% 
  
  unite(year, t,auc,sep = " year AUC: ")


year =factor(survivalROC_data1$year)
survivalROC_data1 %>%
  
  ggplot(mapping = aes(x = FP, y = TP)) +
  scale_color_manual(
    values = c("#BC392F", "#076DB2", "#DD822D")) +
  
  geom_path(aes(color= year))+
  
  geom_abline(intercept = 0, slope = 1, linetype = "dashed")+
  
  facet_wrap( ~ year) +
  
  theme_bw() +
  labs(x = "False positive fraction", y = "True positive fraction", title = "Test ROC") +
  theme(axis.text.x = element_text(vjust = 0.5),
        
        legend.key = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 16),
        legend.position = "none",
        panel.background = element_rect(colour = 'black', size = 1, fill = 'white'),
        panel.grid = element_blank(),
        plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
        legend.text = element_text(size = 12),
        axis.title = element_text(size = 15),
        text = element_text(face = "bold"))

ggsave(filename = "02.test_ROC.pdf", height = 4, width = 8)
ggsave(filename = "02.test_ROC.png", height = 4, width = 8)
