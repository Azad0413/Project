rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/BJTC-420-1/")
if (! dir.exists("./03_intersect")){
  dir.create("./03_intersect")
}
setwd("./03_intersect")
diff1 <- read.delim2('../01_DEGs(COVID)/DEG_sig.xls',row.names = 1)
diff2 <- read.delim2('../02_DEGs(RA)/DEG_sig(GSE55457).xls')
inter <- data.frame(symbol=intersect(rownames(diff1),rownames(diff2)))
##102
write.table(inter,file = 'intersect.xls',sep = '\t',row.names = F,quote = F)

library(ggvenn)
mydata<-list('DEGs(COVID 19)'=rownames(diff1),'DEGs(RA)'=rownames(diff2))
pdf('01.venn.pdf',w=5,h=5)
ggvenn(mydata,c('DEGs(COVID 19)','DEGs(RA)'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       text_size = 4,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()
png('01.venn.png',w=400,h=400)
ggvenn(mydata,c('DEGs(COVID 19)','DEGs(RA)'),
       fill_color = c("#ffb2b2","#b2e7cb"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       text_size = 4,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb"),
       text_color = 'black')
dev.off()

