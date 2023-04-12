rm(list = ls())
# 01 获取数据集--------------
setwd("/data/nas1/luchunlin/project/BJTC-356/")
if (! dir.exists("./06_quadrant")){
  dir.create("./06_quadrant")
}
setwd("./06_quadrant")
gene <- read.table("../02_DEGs/DEG_sig.xls",sep = "\t",header = T,check.names = F)
gene1 <- gene
gene1$id <- rownames(gene1)
gene2 <- gene1[,c(5,1,2)]

protein <- read.table("../04_DEPs/sig_DEP.xls",sep = "\t",header = T,check.names = F)
protein1 <- protein
protein1$id <- rownames(protein1)
protein2 <- protein1[,c(5,1,2)]


#ID转换
library(org.Hs.eg.db)
library(clusterProfiler)
gene_RNA <- bitr(gene2$id,
                 fromType = "SYMBOL",
                 toType = c("ENSEMBL"),
                 OrgDb = "org.Hs.eg.db")
names(gene_RNA)[1] <- "id"


RNA <- merge(gene_RNA,gene2,by ="id")
RNA <- RNA[,c(2,1,3,4)]
names(RNA)[2] <- "symbol"
names(RNA)[1] <- "id"
library(readxl)
# unmap <- read_xlsx('unmap.xlsx')
# uniprot <- read_xlsx('/data/nas1/luchunlin/project/BJX-309/uniprot-filtered-organism__Homo+sapiens+(Human)+[9606]_.xlsx')
# uniprot.unmap <- uniprot[uniprot$Entry%in%unmap$id,]
# uniprot.unmap <- uniprot.unmap[,c(1,4,5)]%>%as.data.frame()
# class(uniprot.unmap)
# rownames(uniprot.unmap) <- uniprot.unmap$Entry
# uniprot.unmap <- uniprot.unmap[unmap$id,]
# write.table(uniprot.unmap,file = 'unmap.uniprot.xls',sep = '\t',row.names = F,quote = F)

gene_protein <- read_tsv('uniprot-compressed_true_download_true_format_tsv-2023.01.10-05.49.23.42.tsv')

# gene_protein <- bitr(protein2$id,
#                      fromType = "UNIPROT",
#                      toType = c("ENSEMBL"),
#                      OrgDb = "org.Hs.eg.db")
# names(gene_protein)[1] <- "id"

protein <- merge(gene_protein,protein2,by ="id")
protein <- protein[,c(2,1,3,4)]
names(protein)[2] <- "symbol"
names(protein)[1] <- "id"




#整合两个表格================================================================
#合并两个表格。选择表头为"id"的列进行合并(取交集)；
#suffixes：如果两个表列名相同，则会在列名后加入suffixes（后缀）参数中对应的后缀；
#all.x=FALSE,all.y=FALSE，表示输出的是x,y表格的交集；
combine= merge(RNA,protein,by.x="id",by.y="id",suffixes = c("_RNA","_Protein") ,all.x=FALSE,all.y=FALSE)
#查看合并后表格的维度；
dim(combine)
#[1] 93 9
#保存合并后的数据；
write.csv(combine,"RNA_protein.csv",row.names=FALSE)

library(dplyr)
library(ggplot2)
library(ggrepel)

data <- data.frame(combine[,c(1:3,5,6)])

#对数据进行分组；
#生成显著上下调数据标签；
data$part <- case_when(abs(data$logFC_RNA) >= 0.5 & abs(data$logFC_Protein) >= 0 ~ "part1379",
                       abs(data$logFC_RNA) < 0.5 & abs(data$logFC_Protein) > 0 ~ "part28",
                       abs(data$logFC_RNA) > 0.5 & abs(data$logFC_Protein) < 0 ~ "part46",
                       abs(data$logFC_RNA) < 0.5 & abs(data$logFC_Protein) < 0 ~ "part5")

data$part <- case_when(data$logFC_RNA <= 0.5 & data$logFC_Protein >= 0 ~ "part1",
                       data$logFC_RNA >= 0.5 & data$logFC_Protein >= 0 ~ "part3",
                       data$logFC_RNA <= -0.5 & data$logFC_Protein <= 0 ~ "part7",
                       data$logFC_RNA < -0.5 & data$logFC_Protein > 0 ~ "part2",
                       data$logFC_RNA < -0.5 & data$logFC_Protein < 0 ~ "part8",
                       abs(data$logFC_RNA) > 0.5 & abs(data$logFC_Protein) < 0 ~ "part46",
                       abs(data$logFC_RNA) < -0.5 & abs(data$logFC_Protein) < 0 ~ "part5",
                       data$logFC_RNA > 0.5 & data$logFC_Protein < 0 ~ "part9",)



table(data$part) 
# part1 part3 part7 part9 
# 6     8     3     5 
head(data)

#开始尝试绘图；
p0 <-ggplot(data,aes(logFC_RNA,logFC_Protein,color=part))
#添加散点；
p1 <- p0+geom_point(size=3)+guides(color="none")
p1
#改变点颜色
# mycolor <- c("#FF9999","#99CC00","#c77cff","gray80")
# p2 <- p1 + scale_colour_manual(name="",values=alpha(mycolor,0.7))
# p2
#添加辅助线；
p3 <- p1+geom_hline(yintercept = c(0,0),
                    size = 1,
                    color = "grey40",
                    lty = "dashed")+
  geom_vline(xintercept = c(-0.5,0.5),
             size = 1,
             color = "grey40",
             lty = "dashed")
p3
#调整横轴和纵轴绘图区域的范围；
#设置y轴范围（上下两端的空白区域设为1），修改刻度标签；
#expansion函数设置坐标轴范围两端空白区域的大小；mult为“倍数”模式，add为“加性”模式；
p4<-p3+
  scale_y_continuous(expand=expansion(add = c(0.1, 0.1)),
                     limits = c(-1, 1)
                     # breaks = c(-0.3,-3,0,3,6),
                     # label = c("-6","-3","0","3","6")
  )+
  scale_x_continuous(expand=expansion(add = c(0.5, 0.5)),
                     limits = c(-6, 6)
                     # breaks = c(-6,-3,0,3,6),
                     # label = c("-6","-3","0","3","6")
  )
p4
#自定义图表主题，对图表做精细调整；
top.mar=0.2
right.mar=0.2
bottom.mar=0.2
left.mar=0.2
#隐藏纵轴，并对字体样式、坐标轴的粗细、颜色、刻度长度进行限定；
mytheme<-theme_bw()+
  theme(text=element_text(family = "sans",colour ="gray30",size = 17),
        panel.grid = element_blank(),
        panel.border = element_rect(size = 0.8,colour = "gray30"),
        axis.line = element_blank(),
        axis.ticks = element_line(size = 0.6,colour = "gray30"),
        axis.ticks.length = unit(1.5,units = "mm"),
        plot.margin=unit(x=c(top.mar,right.mar,bottom.mar,left.mar),
                         units="inches"))
#应用自定义主题；
p5 <- p4+mytheme
p5
ggsave(filename = "02.nqd.png", width = 7, height = 6)
ggsave(filename = "02.nqd.pdf", width = 7, height = 6)
dev.off()

data1 <- subset(data,data$part == "part3"|data$part=='part7')
dim(data1)
write.table(data1,file = 'RNA_protein2.xls',sep = '\t',row.names = F,quote = F)
symbol <- data.frame(symbol=data1$symbol_RNA)
symbol <- data.frame(symbol=symbol[!duplicated(symbol$symbol),])

write.table(symbol,file = 'final.gene.xls',sep = '\t',row.names = F,quote = F)


