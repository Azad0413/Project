rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/XA-0214-1/")
if (! dir.exists("./14_supplement")){
  dir.create("./14_supplement")
}
setwd("./14_supplement")

library(readxl)
stat3 <- read_xlsx('stat3gene.xlsx')
##同源基因转化
## 将基因转换为大鼠的同源基因
# install.packages('homologene')
library(homologene)
homologene::taxData
homologene <- homologene(stat3$Symbol,inTax = 9606,outTax = 10116)
dim(homologene)
diff <- read.delim2('../01_DEGs/DEG_sig(GSE97537).xls')

intersect <- data.frame(symbol=intersect(homologene$`10116`,rownames(diff)))

write.table(intersect,file = 'intersect.xls',sep = '\t',row.names = T,quote = F)

hubgene <- read.delim2('../04_model/hubgene.xls')
pathway <- rbind(hubgene,intersect)
pathway <- pathway[!duplicated(pathway$symbol),]%>%as.data.frame()
diff.path <- diff[pathway$.,]
write.table(diff.path,file = 'diffpath.xls',sep = '\t',row.names = T,quote = F)
###
library(ggvenn)
mydata<-list('DEGs'=rownames(diff),'Stat3 Pathway'=homologene$`10116`)
pdf('01.venn.pdf',w=6,h=6)
ggvenn(mydata,c('DEGs','Stat3 Pathway'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()
png('01.venn.png',w=500,h=500)
ggvenn(mydata,c('DEGs','Stat3 Pathway'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()
