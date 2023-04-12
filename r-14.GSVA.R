rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-300-8/")
if (! dir.exists("./14_GSVA")){
  dir.create("./14_GSVA")
}
setwd("./14_GSVA")

dat<-read.delim2('../00_rawdata/dat.fpkm.xls',row.names = 1)%>%lc.tableToNum()
colnames(dat)<-gsub('.','-',colnames(dat),fixed = T)
group<-read.delim2('../08_risk/risk.xls')%>%dplyr::select(c('id','riskScore'))
colnames(group)<-c('sample','group')
group$group <- ifelse(group$group>median(group$group),'High_risk','Low_risk')
High.sample <- group$sample[which(group$group=='High_risk')]
group2 <- group
colnames(group)
colnames(group2)<-c('id','label')
gsva_exp<-dat[,group2$id]
all(colnames(gsva_exp) == group2$id)
dim(gsva_exp)
## 19712    85
# 分组
library(limma)
library(GSEABase)
library(GSVA)
group<- group2$label %>% as.factor()
design <- model.matrix(~0 + group)
rownames(design) <- colnames(gsva_exp)
colnames(design) <- levels(group)
compare<- makeContrasts("High_risk-Low_risk", levels = design)
KEGG_ref <- getGmt("/data/nas1/luchunlin/pipeline/GSVA/c2.cp.kegg.v7.4.symbols.gmt")
es_KEGG <- gsva(as.matrix(gsva_exp), KEGG_ref,
                min.sz=10, max.sz=500, verbose=TRUE)
es_KEGG <- as.data.frame(es_KEGG)
fit <- lmFit(es_KEGG, design)
fit2 <- contrasts.fit(fit ,compare)
fit3 <- eBayes(fit2)
allGeneSets <- topTable(fit3, coef = 1, number = Inf)
logFCcutoff <- 0
allGeneSets$change = as.factor(
  ifelse(allGeneSets$adj.P.Val < 0.05 & abs(allGeneSets$logFC) > logFCcutoff,
         ifelse(allGeneSets$logFC > logFCcutoff,'UP','DOWN'),'NOT')
)

DEGeneSets <- subset(allGeneSets,
                     allGeneSets$adj.P.Val < 0.05 & abs(allGeneSets$t) > logFCcutoff )

write.table(allGeneSets, 
            file = "KEGG.xls",
            quote = F,
            sep = "\t",
            row.names = T)
DEGeneSets <- DEGeneSets[order(DEGeneSets$adj.P.Val),]
dim(DEGeneSets)
# [1]31
write.table(DEGeneSets,
            file = 'diff_KEGG.xls',
            quote = F,
            sep = "\t",
            row.names = T)
### 发散条形图绘制

#install.packages('ggprism')
library(ggprism)
## barplot
dat_plot<-data.frame(id=rownames(allGeneSets),
                     t=allGeneSets$t)
allGeneSets <- allGeneSets[dat_plot$id,]
dat_plot$threshold = factor(allGeneSets$change)
# dat_plot$threshold = factor(ifelse(dat_plot$p  >-1, ifelse(dat_plot$t >= 1 ,'Up','NoSignifi'),'Down'),levels=c('Up','Down','NoSignifi'))
table(dat_plot$threshold)
dat_plot$threshold <- ifelse(dat_plot$threshold=='UP','Up',
                             ifelse(dat_plot$threshold=='DOWN','Down','NoSignifi'))
dat_plot<-dat_plot%>%arrange(t)

dat_plot$id<-factor(dat_plot$id,levels=dat_plot$id)

p <- ggplot(data = dat_plot,aes(x = id,y = t,fill = threshold)) +
  geom_col()+
  coord_flip() +
  scale_fill_manual(values = c('Up'= '#36638a','NoSignifi'='#cccccc','Down'='#7bcd7b')) +
   geom_hline(yintercept = c(-1,1),color = 'white',size = 0.5,lty='dashed') +
  scale_y_continuous(breaks = c(-6,-5,-4,-3,-2,-1,0,1,2,3,4,5)) +
  xlab('') + 
  ylab('t value of GSVA score, tumour versus non-malignant') + #注意坐标轴旋转了
  guides(fill=F)+ # 不显示图例
  theme_prism(border = T) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )
p
# 小于-2的数量
low1 <- dat_plot %>% filter(t < -1) %>% nrow()
# 小于0总数量
low0 <- dat_plot %>% filter( t < 0) %>% nrow()
# 小于2总数量
high0 <- dat_plot %>% filter(t < 1) %>% nrow()
# 总的柱子数量
high1 <- nrow(dat_plot)

# 依次从下到上添加标签
p <- p + geom_text(data = dat_plot[1:low1,],aes(x = id,y = 0.1,label = id),
                   hjust = 0,color = 'black')  + # 小于-1的为黑色标签
  geom_text(data = dat_plot[(low1 +1):low0,],aes(x = id,y = 0.1,label = id),
            hjust = 0,color = 'grey') + # 灰色标签
  geom_text(data = dat_plot[(low0 + 1):high0,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'grey') + # 灰色标签
  geom_text(data = dat_plot[(high0 +1):high1,],aes(x = id,y = -0.1,label = id),
            hjust = 1,color = 'black') # 大于1的为黑色标签
p
ggsave(filename = '01.GSVA.pdf',p,w=15,h=24)
ggsave(filename = '01.GSVA.png',p,w=15,h=26)
