rm(list = ls())
setwd("/data/nas1/luchunlin/project/BJTC-258/")
if (! dir.exists("./10_immune_infiltration")){
  dir.create("./10_immune_infiltration")
}
setwd("./10_immune_infiltration")

dat.tcga<-read.delim2("/data/nas1/luchunlin/project/BJTC-258/00_rawdata/dat.fpkm.xls", row.names = 1)%>% lc.tableToNum
colname<-data.frame(sample=colnames(dat.tcga))
colname$sample<-gsub('.','-',colname$sample,fixed = T)
colnames(dat.tcga)<-colname$sample
group<-read.delim2('/data/nas1/luchunlin/project/BJTC-258/05_survival/clinical.xls')%>%
  dplyr::select(c('sample','group'))
high.sample<-group$sample[which(group$group=='High')]

library(immunedeconv)
library(RColorBrewer)
library(tidyverse)
mypalette <- colorRampPalette(brewer.pal(7,"Paired"))
# PART A------
## 01 quantiseq-----
res.quantiseq<-deconvolute(as.matrix(dat.tcga),method = 'quantiseq')
write.table(res.quantiseq,"quantiseq.txt",sep = '\t',col.names = T,row.names = F,quote = F)
save(file = 'quantiseq.Rdata',res.quantiseq)
##画图
pdf('01.quantiseq.box.pdf',w=8,h=6)
res.quantiseq %>%
  gather(sample, fraction, -cell_type) %>%
  merge(group,by='sample')%>%
  # 绘制堆积条形图
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
  geom_bar(position = 'stack',stat = 'identity')+
  scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  labs(x='',
       y='Relative Percent',
       fill='')+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'top') +
  scale_fill_manual(values = mypalette(14))+
  facet_grid(~group,scales= "free",space= "free")
dev.off()
png('01.quantiseq.box.png',w=600,h=500)
res.quantiseq %>%
  gather(sample, fraction, -cell_type) %>%
  merge(group,by='sample')%>%
  # 绘制堆积条形图
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
  geom_bar(position = 'stack',stat = 'identity')+
  scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  labs(x='',
       y='Relative Percent',
       fill='')+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'top') +
  scale_fill_manual(values = mypalette(14))+
facet_grid(~group,scales= "free",space= "free")
dev.off()
## 差异-------
dat.quantiseq <- res.quantiseq %>% 
  tibble::column_to_rownames(var = "cell_type") %>% 
  t %>% as.data.frame() %>% 
  tibble::rownames_to_column(var = "sample")
dat.quantiseq <- merge(group, dat.quantiseq, by = "sample")
dat.quantiseq2 <- tidyr::gather(dat.quantiseq, ImmuneCell, Score, -c("sample", "group"))
library(rstatix)
stat_quantiseq <- dat.quantiseq2 %>% 
  group_by(ImmuneCell) %>% 
  wilcox_test(Score ~ group) %>% 
  adjust_pvalue(method = "BH") %>%  # method BH == fdr
  add_significance("p")
stat_quantiseq
write.table(stat_quantiseq,file = 'stat.quantiseq.xls',sep = '\t',row.names = F,quote = F)
colnames(dat.quantiseq2)
violin.quantiseq<-dat.quantiseq2[dat.quantiseq2$ImmuneCell%in%stat_quantiseq$ImmuneCell[which(stat_quantiseq$p<0.05)],]
quantiseq_plot <- ggplot(violin.quantiseq, aes(x=ImmuneCell,
                                      y=Score,
                                      fill=group)) +
  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  #stat_boxplot(geom="errorbar", 
  #             width=0.1,
  #             position = position_dodge(0.9)) +
  #geom_boxplot(width=0.7,
  #             position=position_dodge(0.9),
  #             outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#FF6A6A", "#20B2AA"), name = "Group")+
  labs(title="Immune Cell", x="", y = "Score",size=20) +
  stat_compare_means(data = violin.quantiseq,
                     mapping = aes(group = group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F) +
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18),
        axis.text.x=element_text(angle=0,hjust=,colour="black",face="bold",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+facet_wrap(~ImmuneCell,scales = "free",nrow = 2) 
quantiseq_plot
ggsave(filename = '02.quantiseq.plot.pdf',quantiseq_plot,w=8,h=6)
ggsave(filename = '02.quantiseq.plot.png',quantiseq_plot,w=8,h=6)
## 02 xcell-------
#res.xcell <- deconvolute(as.matrix(dat.tcga), method = "xcell")
library(xCell)
res.xcell <- xCellAnalysis(expr = as.matrix(dat.tcga),rnaseq = T)
write.table(res.xcell,'xcell.txt',sep = '\t',col.names = T,row.names = F,quote = F)
save(file = 'xcell.Rdata',res.xcell)
res.xcell<-res.xcell[c(1:64),]
pdf('03.xcell.box.pdf',w=9,h=8)
res.xcell %>%as.data.frame()%>%rownames_to_column(var = 'cell_type')%>%
  gather(sample, fraction, -cell_type) %>%
  merge(group,by='sample')%>%
  # 绘制堆积条形图
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
  geom_bar(position = 'stack',stat = 'identity')+
  scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  labs(x='',
       y='scores',
       fill='')+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'top') +
  scale_fill_manual(values = mypalette(64))+
  facet_grid(~group,scales= "free",space= "free")
dev.off()

png('03.xcell.box.png',w=700,h=600)
res.xcell %>%as.data.frame()%>%rownames_to_column(var = 'cell_type')%>%
  gather(sample, fraction, -cell_type) %>%
  merge(group,by='sample')%>%
  # 绘制堆积条形图
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
  geom_bar(position = 'stack',stat = 'identity')+
  scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  labs(x='',
       y='scores',
       fill='')+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'top') +
  scale_fill_manual(values = mypalette(64))+
  facet_grid(~group,scales= "free",space= "free")
dev.off()

## 差异-------
dat.xcell <- res.xcell %>% 
 # tibble::column_to_rownames(var = "cell_type") %>% 
  t %>% as.data.frame() %>% 
  tibble::rownames_to_column(var = "sample")
dat.xcell <- merge(group, dat.xcell, by = "sample")
dat.xcell2 <- tidyr::gather(dat.xcell, ImmuneCell, Score, -c("sample", "group"))
library(rstatix)
stat_xcell <- dat.xcell2 %>% 
  group_by(ImmuneCell) %>% 
  wilcox_test(Score ~ group) %>% 
  adjust_pvalue(method = "BH") %>%  # method BH == fdr
  add_significance("p")
write.table(stat_xcell,file = 'stat.xcell.xls',sep = '\t',row.names = F,quote = F)
violin.xcell<-dat.xcell2[dat.xcell2$ImmuneCell%in%stat_xcell$ImmuneCell[which(stat_xcell$p<0.05)],]
xcell_plot <- ggplot(violin.xcell, aes(x=ImmuneCell,
                                      y=Score,
                                      fill=group)) +
  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  #stat_boxplot(geom="errorbar", 
  #             width=0.1,
  #             position = position_dodge(0.9)) +
  #geom_boxplot(width=0.7,
  #             position=position_dodge(0.9),
  #             outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#FF6A6A", "#20B2AA"), name = "Group")+
  labs(title="Immune Cell(Xcell)", x="", y = "Score",size=20) +
  stat_compare_means(data = violin.xcell,
                     mapping = aes(group = group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F) +
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18),
        axis.text.x=element_text(angle=0,hjust=,colour="black",face="bold",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+facet_wrap(~ImmuneCell,scales = "free",nrow = 5) 
xcell_plot
ggsave(filename = '04.xcell.plot.pdf',xcell_plot,w=15,h=12)
ggsave(filename = '04.xcell.plot.png',xcell_plot,w=15,h=12)
## 03 mmcp_counter------
res.mcp_counter <- deconvolute(as.matrix(dat.tcga), method = "mcp_counter")
write.table(res.mcp_counter,'mcp_counter.txt',sep = '\t',col.names = T,row.names = F,quote = F)
save(file = 'mcp_counter.Rdata',res.mcp_counter)
pdf('05.mcp_counter.pdf',w=8,h=6)
res.mcp_counter %>%
  gather(sample, fraction, -cell_type) %>%
  merge(group,by='sample')%>%
  # 绘制堆积条形图
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
  geom_bar(position = 'stack',stat = 'identity')+
  scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  labs(x='',
       y='scores',
       fill='')+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'top') +
  scale_fill_manual(values = mypalette(11))+
  facet_grid(~group,scales= "free",space= "free")
dev.off()

png('05.mcp_counter.png',w=600,h=500)
res.mcp_counter %>%
  gather(sample, fraction, -cell_type) %>%
  merge(group,by='sample')%>%
  # 绘制堆积条形图
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
  geom_bar(position = 'stack',stat = 'identity')+
  scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  labs(x='',
       y='scores',
       fill='')+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'top') +
  scale_fill_manual(values = mypalette(11))+
  facet_grid(~group,scales= "free",space= "free")
dev.off()
## 差异-------
dat.mcp_counter <- res.mcp_counter %>% 
  tibble::column_to_rownames(var = "cell_type") %>% 
  t %>% as.data.frame() %>% 
  tibble::rownames_to_column(var = "sample")
dat.mcp_counter <- merge(group, dat.mcp_counter, by = "sample")
dat.mcp_counter2 <- tidyr::gather(dat.mcp_counter, ImmuneCell, Score, -c("sample", "group"))
library(rstatix)
stat_mcp_counter <- dat.mcp_counter2 %>% 
  group_by(ImmuneCell) %>% 
  wilcox_test(Score ~ group) %>% 
  adjust_pvalue(method = "BH") %>%  # method BH == fdr
  add_significance("p")

write.table(stat_mcp_counter,file = 'stat.mcp_counter.xls',sep = '\t',row.names = F,quote = F)
violin.mcp_counter<-dat.mcp_counter2[dat.mcp_counter2$ImmuneCell%in%stat_mcp_counter$ImmuneCell[which(stat_mcp_counter$p<0.05)],]
mcp_counter_plot <- ggplot(violin.mcp_counter, aes(x=ImmuneCell,
                                      y=Score,
                                      fill=group)) +
  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  #stat_boxplot(geom="errorbar", 
  #             width=0.1,
  #             position = position_dodge(0.9)) +
  #geom_boxplot(width=0.7,
  #             position=position_dodge(0.9),
  #             outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#FF6A6A", "#20B2AA"), name = "Group")+
  labs(title="Immune Cell(mcp_counter)", x="", y = "Score",size=20) +
  stat_compare_means(data = violin.mcp_counter,
                     mapping = aes(group = group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F) +
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18),
        axis.text.x=element_text(angle=0,hjust=,colour="black",face="bold",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+facet_wrap(~ImmuneCell,scales = "free",nrow = 3) 
mcp_counter_plot
ggsave(filename = '06.mcp_counter.plot.pdf',mcp_counter_plot,w=9,h=9)
ggsave(filename = '06.mcp_counter.plot.png',mcp_counter_plot,w=9,h=9)
## 04 timer---------

res.timer <- deconvolute_timer(as.matrix(dat.tcga), indications = rep('ov',378))
write.table(res.timer,'timer.txt',sep = '\t',col.names = T,row.names = F,quote = F)
save(file = 'timer.Rdata',res.timer)
pdf('07.timer.pdf',w=8,h=6)
res.timer %>%as.data.frame()%>%rownames_to_column(var = 'cell_type')%>%
  gather(sample, fraction, -cell_type) %>%
  merge(group,by='sample')%>%
  # 绘制堆积条形图
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
  geom_bar(position = 'stack',stat = 'identity')+
  scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  labs(x='',
       y='scores',
       fill='')+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'top') +
  scale_fill_manual(values = mypalette(7))+
  facet_grid(~group,scales= "free",space= "free")
dev.off()
png('07.timer.png',w=600,h=500)
res.timer %>%as.data.frame()%>%rownames_to_column(var = 'cell_type')%>%
  gather(sample, fraction, -cell_type) %>%
  merge(group,by='sample')%>%
  # 绘制堆积条形图
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
  geom_bar(position = 'stack',stat = 'identity')+
  scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  labs(x='',
       y='scores',
       fill='')+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'top') +
  scale_fill_manual(values = mypalette(7))+
  facet_grid(~group,scales= "free",space= "free")
dev.off()
dat.timer <- res.timer %>% 
#  tibble::column_to_rownames(var = "cell_type") %>% 
  t %>% as.data.frame() %>% 
  tibble::rownames_to_column(var = "sample")
dat.timer <- merge(group, dat.timer, by = "sample")
dat.timer2 <- tidyr::gather(dat.timer, ImmuneCell, Score, -c("sample", "group"))
library(rstatix)
stat_timer <- dat.timer2 %>% 
  group_by(ImmuneCell) %>% 
  wilcox_test(Score ~ group) %>% 
  adjust_pvalue(method = "BH") %>%  # method BH == fdr
  add_significance("p")
write.table(stat_timer,file = 'stat.timer.xls',sep = '\t',row.names = F,quote = F)
violin.timer<-dat.timer2[dat.timer2$ImmuneCell%in%stat_timer$ImmuneCell[which(stat_timer$p<0.05)],]
timer_plot <- ggplot(violin.timer, aes(x=ImmuneCell,
                                      y=Score,
                                      fill=group)) +
  geom_violin(trim=F,color="black") + #绘制小提琴图, “color=”设置小提琴图的轮廓线的颜色(#不要轮廓可以设为white以下设为背景为白色，其实表示不要轮廓线)
  #"trim"如果为TRUE(默认值),则将小提琴的尾部修剪到数据范围。如果为FALSE,不修剪尾部。
  #stat_boxplot(geom="errorbar", 
  #             width=0.1,
  #             position = position_dodge(0.9)) +
  #geom_boxplot(width=0.7,
  #             position=position_dodge(0.9),
  #             outlier.shape = NA)+ #绘制箱线图，此处width=0.1控制小提琴图中箱线图的宽窄
  scale_fill_manual(values= c("#FF6A6A", "#20B2AA"), name = "Group")+
  labs(title="Immune Cell(TIMER)", x="", y = "Score",size=20) +
  stat_compare_means(data = violin.timer,
                     mapping = aes(group = group),
                     label ="p.signif",
                     method = 'wilcox.test',
                     paired = F) +
  theme_bw()+
  theme(plot.title = element_text(hjust =0.5,colour="black",face="bold",size=18),
        axis.text.x=element_text(angle=0,hjust=,colour="black",face="bold",size=10), 
        axis.text.y=element_text(hjust=0.5,colour="black",face="bold",size=12), 
        axis.title.x=element_text(size=16,face="bold"),
        axis.title.y=element_text(size=16,face="bold"),
        legend.text=element_text(face = "bold", hjust = 0.5,colour="black", size=12),
        legend.title = element_text(face = "bold", size = 12),
        legend.position = "top",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+facet_wrap(~ImmuneCell,scales = "free",nrow = 2) 
timer_plot
ggsave(filename = '08.timer.plot.pdf',timer_plot,w=8,h=6)
ggsave(filename = '08.timer.plot.png',timer_plot,w=8,h=6)
# PART B-------
## 01 cancer-immunity cycle--------
library(GSVA)
gene_set <- read.table("cancer-immunity cycle.txt",
                       header = T,
                       sep ="\t")
table(gene_set$Steps)
gene_set$Steps<-c(rep('Step 1:Relase of cancer cell antigenes',40),
                  rep('Step 2:Cancer antigen presentation',36),
                  rep('Step 3:Priming and activation',62),
                  rep('Step 4:Trafficking of immune cells to tumors',76),
                  rep('Step 5:Infiltration of immune cells into tumors',9),
                  rep('Step 6:Recongnition of cancer cells by T cells',25),
                  rep('Step 7:Killing of cancer cells',25))
dat.final <- as.matrix(dat.tcga)
gene_list <- split(as.matrix(gene_set)[,1],
                   gene_set[,2])

ssgsea_score = gsva(dat.final, gene_list, 
                    method = "ssgsea", 
                    ssgsea.norm = TRUE, 
                    verbose = TRUE)
write.table(ssgsea_score,
            file = "ssgsea_result(seven step).xls",
            sep = "\t",
            quote = F)
# PART C-------
## 01 TIDE评分------
# fpkm转TPM
FPKM2TPM <- function(fpkm){
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}
# tide_dat[which(tide_dat<0)] <- 0
tide_dat <- apply(dat.tcga,2,FPKM2TPM)
tide_dat <- log2(tide_dat + 1)
rownmean <- apply(tide_dat,1,mean)
tide_dat2 <- sweep(tide_dat, 1, rownmean)
dim(tide_dat2)
write.table(tide_dat2,
            file ="tide_dat.txt",
            sep = "\t",
            quote = F,
            row.names = T)
tide_result <- read.csv("tide.result.csv",header = T)
colnames(tide_result)

# PART D（heatmap绘图）-----
### 4种免疫浸润结果

timer.heat<-res.timer[stat_timer$ImmuneCell[which(stat_timer$p<0.05)],]
rownames(timer.heat)<-paste0(rownames(timer.heat),'_TIMER')

quantiseq.heat<-res.quantiseq[res.quantiseq$cell_type%in%stat_quantiseq$ImmuneCell[which(stat_quantiseq$p<0.05)],]%>%
  column_to_rownames(var = 'cell_type')
rownames(quantiseq.heat)<-paste0(rownames(quantiseq.heat),'_quantiseq')  
xcell.heat<-res.xcell[rownames(res.xcell)%in%stat_xcell$ImmuneCell[which(stat_xcell$p<0.05)],]
rownames(xcell.heat)<-paste0(rownames(xcell.heat),'_Xcell')
mcp_counter.heat<-res.mcp_counter[res.mcp_counter$cell_type%in%stat_mcp_counter$ImmuneCell[which(stat_mcp_counter$p<0.05)],]%>%
  column_to_rownames(var = 'cell_type')
rownames(mcp_counter.heat)<-paste0(rownames(mcp_counter.heat),'_MCPcounter')
tiic.heat = rbind(timer.heat,quantiseq.heat,xcell.heat,mcp_counter.heat)
### seven step
seven.step.heat<-ssgsea_score
### estimate的4个评分
estimate<-read.delim2('/data/nas1/luchunlin/project/BJTC-258/09_estimate/es_score.xls',row.names = 1)%>%
  dplyr::select(-Description)%>%t
rownames(estimate)<-colname$sample
### TIDE 评分
tide_score <- subset(tide_result, select = c("Patient", "TIDE","Dysfunction","Exclusion","CTL.flag"))%>%
  column_to_rownames(var = 'Patient')
tide_score<-tide_score[rownames(estimate),]
## 热图---------
library(ComplexHeatmap)
library(ggplot2)
library(RColorBrewer)
library(gplots)

rt_group <- cbind(estimate,tide_score)
rt_group$group <- ifelse(rownames(rt_group) %in% high.sample,
                         "High", "Low")
rt_group<-rt_group[order(rt_group$group),]
rt_group$group<-as.factor(rt_group$group)
rt_dat<-rbind(tiic.heat,seven.step.heat)
rt_dat<-rt_dat[,group$sample]
rt_dat<-t(scale(t(rt_dat)))
rt_dat[rt_dat>2]<-2
rt_dat[rt_dat<-2]<- -2
# x[x>2] <- 2
# x[x<-2] <- 2

#rt_dat<-log2(rt_dat+0.01)
rt_dat<-as.matrix(rt_dat)

annotation_col<-rt_group
colnames(annotation_col)

library(circlize)
annotation_col$StromalScore<-as.numeric(annotation_col$StromalScore)
annotation_col$ImmuneScore<-as.numeric(annotation_col$ImmuneScore)
annotation_col$ESTIMATEScore<-as.numeric(annotation_col$ESTIMATEScore)
annotation_col$TumorPurity<-as.numeric(annotation_col$TumorPurity)
top_annotation = HeatmapAnnotation(TIDE = annotation_col$TIDE,
                              Dysfunction = anno_lines(annotation_col$Dysfunction),
                              Exclusion = anno_lines(annotation_col$Exclusion),
                              CTL.flag=annotation_col$CTL.flag,
                              StromalScore=annotation_col$StromalScore,
                              ImmuneScore=annotation_col$ImmuneScore,
                              ESTIMATEScore=annotation_col$ESTIMATEScore,
                              TumorPurity=annotation_col$TumorPurity,
                              group=annotation_col$group,
                            col = list(TIDE=colorRamp2(c(-4, 0, 4), c("lightblue", "white", "tomato")),
                                       StromalScore=colorRamp2(c(-1000,0,1000),c('purple','white','pink')),
                                       CTL.flag=c('False'='blue','True'='red'),
                                       ImmuneScore=colorRamp2(c(-2000,-1000,0,1000,2000),c('powderblue','skyblue','white','pink','deeppink')),
                                       ESTIMATEScore=colorRamp2(c(-4000,-2000,0,2000,4000),c('powderblue','mediumturquoise','white','coral','tomato')),
                                       TumorPurity=colorRamp2(c(0.2,0.4,0.6,0.8,1),c('white','mistyrose','lightpink','hotpink','deeppink')),
                                       group=c('High'='coral','Low'='lightseagreen')
                            ))
pdf(file = '09.heatmap.pdf',w=12,h=10)
ht.list=Heatmap(rt_dat, name = "Score", top_annotation = top_annotation,row_names_gp = gpar(fontsize = 8),
                row_split = factor(c(rep('TIMER',6),rep('quantiseq',6),rep('xCell',34),rep('MCP_counter',8),rep('cancer-immune cycle',7)),
                                   levels = c('TIMER','quantiseq','MCP_counter','xCell','cancer-immune cycle')),
                show_column_names = F,cluster_rows = F,cluster_columns = F,
                row_title_rot = 0
                )

draw(ht.list,heatmap_legend_side='left',annotation_legend_side='left')
dev.off()

png(file = '09.heatmap.png',w=1000,h=800)
ht.list=Heatmap(rt_dat, name = "Score", top_annotation = top_annotation,row_names_gp = gpar(fontsize = 8),
                row_split = factor(c(rep('TIMER',6),rep('quantiseq',6),rep('xCell',34),rep('MCP_counter',8),rep('cancer-immune cycle',7)),
                                   levels = c('TIMER','quantiseq','MCP_counter','xCell','cancer-immune cycle')),
                show_column_names = F,cluster_rows = F,cluster_columns = F,
                row_title_rot = 0
)

draw(ht.list,heatmap_legend_side='left',annotation_legend_side='left')
dev.off()
