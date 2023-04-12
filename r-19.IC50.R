rm(list = ls())
setwd("/data/nas1/luchunlin/project/SJZZK-428-10/")
if (! dir.exists("./19_IC50")){
  dir.create("./19_IC50")
}
setwd("./19_IC50")
dat<-read.delim2('../00_rawdata/dat.fpkm.xls',row.names = 1)%>%lc.tableToNum()
colnames(dat)<-gsub('.','-',colnames(dat),fixed = T)
risk<-read.delim2('../09_risk/risk.xls')
high.sample<-risk$id[which(risk$risk==0)]
low.sample<-risk$id[which(risk$risk==1)]
dat<-log2(dat[,risk$id]+1)
library(pRRophetic)
library(ggplot2)
set.seed(12345)
lasso_geneids<-read.delim2('../08_Lasso/lasso_genes.csv',header = F)
model_expr<-dat[lasso_geneids$V1, risk$id]
riskscore<-data.frame(risk$id,risk$risk)
colnames(riskscore)<-c('sample','risk')
riskscore$risk[which(riskscore$risk==1)] <-'Low risk'
riskscore$risk[which(riskscore$risk==0)] <-'High risk'
head(riskscore)
drug<-read.table(file = '/data/nas1/luchunlin/pipeline/Medicinal_Sensity/drugs.txt',sep='\t',header=F)
ic50<-data.frame(riskscore$sample)
a<-data.frame(row.names=riskscore$sample,riskscore$risk)
colnames(a)<-'risk'
cnt<-1
while (cnt < 139) {
  predictedPtype <- pRRopheticPredict(as.matrix(model_expr), drug[cnt,],selection=1)
  Tipifarnib<-data.frame(predictedPtype)
  colnames(Tipifarnib)<-drug[cnt,]
  a<-cbind(a,Tipifarnib)
  cnt = cnt + 1
}
# a <- a[,c("Paclitaxel", "Docetaxel",
#           "Cisplatin","Gemcitabine","Mitomycin.C","Doxorubicin")]
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
medicinal_result <- t(subset(x, select = -risk)) 
high_group <- risk$id[which(risk$risk==0)]
low_group <- risk$id[which(risk$risk==1)]
pvalue = padj = log2FoldChange <- matrix(0, nrow(medicinal_result), 1)
for (i in 1:nrow(medicinal_result)){
  pvalue[i, 1] = p.value = wilcox.test(medicinal_result[i, high_group],
                                       medicinal_result[i, low_group])$p.value
  log2FoldChange[i, 1] = mean(medicinal_result[i, high_group]) - 
    mean(medicinal_result[i, low_group])
}
padj <- p.adjust(as.vector(pvalue), "fdr", n = length(pvalue))
rTable <- data.frame(log2FoldChange, 
                     pvalue, 
                     padj,
                     row.names = rownames(medicinal_result))
high_group_res <- signif(apply(medicinal_result[rownames(rTable), high_group], 
                               1,
                               median), 4)
low_group_res <- signif(apply(medicinal_result[rownames(rTable), low_group], 
                              1, 
                              median), 4)
rTable <- data.frame(high_group_res, 
                     low_group_res,
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
# rTable <- rTable[c("Paclitaxel", "Docetaxel",
#                    "Cisplatin","Gemcitabine","Mitomycin.C","Doxorubicin"),]
write.table(rTable,
            file = "drugs_wilcox_test.xls",
            quote = F,
            row.names = F)
DE.drug<-rTable[which(rTable$pvalue<0.05),]

##117
## 发散条形图绘制
#install.packages('ggprism')
library(ggprism)
## 横坐标药物，纵坐标：IC50(H)/IC50(L)-1
dat_plot<-data.frame(drug=rownames(rTable),
                     'IC50(H)/IC50(L)-1'=(rTable$high_group_res/rTable$low_group_res-1),
                     pvalue=rTable$pvalue,
                     padj=rTable$padj)
dat_plot$threshold=factor(ifelse(dat_plot$padj<0.05&dat_plot$pvalue<0.05,'P<0.05 & FDR<0.05',ifelse(dat_plot$pvalue<0.05&dat_plot$padj>0.05,'P<0.05 & FDR>0.05','P>0.05 & FDR>0.05')))
dat_plot<-dat_plot%>%arrange(desc(dat_plot$IC50.H..IC50.L..1))

dat_plot$drug<-factor(dat_plot$drug,levels=dat_plot$drug)

p <- ggplot(data = dat_plot,aes(x = drug,y = IC50.H..IC50.L..1,fill = threshold)) +
  geom_col()+
  scale_fill_manual(values = c('P<0.05 & FDR<0.05'= '#5F9EA0','P>0.05 & FDR>0.05'='#cccccc','P<0.05 & FDR>0.05'='#FFD700')) +
  xlab('') +
  ylab('Exp(Median IC50(H))/Exp(Median IC50(L))-1') +
  theme_prism(border = T) +
  theme(
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(angle = 65,size = 7,
                               hjust = 1,vjust = 1),
    axis.text.y = element_text(size = 13),
    legend.position = c(0.85,0.85),
    legend.text = element_text(size = 8,face = 'bold')
  )
p
ggsave(filename = '01.all_medicinal.pdf',w=15,h=8)
ggsave(filename = '01.all_medicinal.png',w=15,h=8)


library(tidyr)
library(ggplot2)
library(ggpubr)
library(Ipaper)
all(rownames(rTable) == rownames(medicinal_result))
drugs_res <- data.frame(drugs=rownames(medicinal_result), medicinal_result, padj=rTable$padj)
# drugs_res <- drugs_res[rownames(rTable2),]
drugs_res <- drugs_res[which(rownames(drugs_res) %in% 
                               c("Paclitaxel", "Docetaxel",
                                 "Cisplatin","Gemcitabine","Mitomycin.C","Doxorubicin")),]
violin_dat <- gather(drugs_res, key=indivs, value=risk, -c("drugs","padj"))
violin_dat$indivs <- ifelse(gsub("\\.","-",violin_dat$indivs) %in% high_group,
                            "High risk", "Low risk") 
violin_dat$indivs <- factor(violin_dat$indivs, levels = c("High risk", "Low risk"))
head(violin_dat)
drugs_hub_boxplot1 <- ggboxplot(violin_dat, x = "indivs", y = "risk",
                                color = "indivs", palette = c("#A73030FF", "#0073C2FF"),
                                add = "jitter",
                                short.panel.labs = T,
                                ggtheme = theme_bw()) +
  stat_compare_means(label = "p.signif", label.x = 1.4, vjust = 0.5)
drugs_hub_boxplot <- facet(drugs_hub_boxplot1,
                           facet.by = "drugs",
                           short.panel.labs = T,
                           panel.labs.background = list(fill = "white"),
                           ncol = 3,
                           scales = "free_y") + xlab("") + ylab("Estimated ln(IC50)") +
  # geom_text(data=data_text,
  #           mapping=aes(x=x,y=y,label=label),nudge_x=0.1,nudge_y=0.1)+
  theme(panel.grid = element_blank(),
        legend.position = "none",
        # strip.background = element_blank(),
        axis.title.y = element_text(size = 18, face = "bold"),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 10, face = "bold"),
        axis.text.x = element_text(size = 13, face = "bold"),
        plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
        text = element_text(size = 13, face = "bold"))
drugs_hub_boxplot
ggsave(filename = "02.drugs.plot.pdf", height = 6, width = 10)
ggsave(filename = "02.drugs.plot.png", height = 6, width = 10)


## 药物靶基因的表达水平--------
target <- read_xlsx('target.xlsx')%>%select('Target')
target <- target[!duplicated(target$Target),]

target.dat <- dat[target$Target,]
hub_exp2<-target.dat
hub_exp2$Symbol<-rownames(hub_exp2)
hub_exp2<-gather(hub_exp2,
                 key = sample,
                 value = expr,
                 -c("Symbol"))

## 样本分组
hub_exp2$Group<-ifelse(hub_exp2$sample%in%low.sample,'Low risk','High risk')
##分面图形
library(ggplot2)
library(tidyverse)
library(ggpubr)
library(ggsci)
library(rstatix)
stat.test<-hub_exp2%>%
  group_by(Symbol)%>%
  t_test(expr ~ Group)%>%
  adjust_pvalue(method = 'fdr')
# stat.test$p.adj<-ifelse(stat.test$p.adj<0.001,"***",
#                         ifelse(stat.test$p.adj<0.05,"**",
#                                ifelse(stat.test$p.adj<0.05,"*",'ns')))
exp_plot <- ggplot(hub_exp2, aes(x=Symbol, 
                                 y=expr,
                                 fill=Group)) +
  #  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  stat_boxplot(geom="errorbar", 
               width=0.1,
               position = position_dodge(0.9)) +
  geom_boxplot(width=0.7,
               position=position_dodge(0.9),
               outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#A73030FF", "#0073C2FF"), name = "Group")+
  labs(title="Target Gene Expression Level between High-Low risk Group", x="", y = "expression level",size=20) +
  stat_compare_means(data = hub_exp2,
                     mapping = aes(group = Group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F) +
  #  geom_signif(comparisons = my_comparisons,
  #              test = t.test,
  #              map_signif_level = T)+
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=15),
        axis.text.x=element_text(angle=45,hjust=1,colour="black",face="bold",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=8),
        legend.title = element_text(face = "bold", size = 10),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
exp_plot
ggsave('03.expression.pdf',exp_plot,width = 7,height = 5)
ggsave('03.expression.png',exp_plot,width = 7,height = 5)
