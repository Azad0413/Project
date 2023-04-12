rm(list = ls())
setwd("/data/nas1/luchunlin/project/NN-0118-2/")
if (! dir.exists("./05_intersect")){
  dir.create("./05_intersect")
}
setwd("./05_intersect")


diff1 <- read.delim2('../01_DEGs(TCGA)/DEG_all.xls',row.names = 1)
diff2 <- read.delim2('../02_DEGs(GSE40275)/DEG_all(GSE40275).xls')
diff3 <- read.delim2('../03_DEGs(GSE21933)/DEG_all(GSE21933).xls')
diff4 <- read.delim2('../04_DEGs(GSE18842)/DEG_all(GSE18842).xls')
# 
# hubgene <- c('TPX2','KIF14','KIF20A','MKI67','ANLN','TOP2A','SFTPB','SFTPC','SFTPD')
# diff1 <- diff1[hubgene,]
# diff2 <- diff2[hubgene,]
# diff3 <- diff3[hubgene,]
# diff4 <- diff4[hubgene,]
###UP-------
table(diff2$change)
up1 <- diff1[which(diff1$change=='UP'),]
up2 <- diff2[which(diff2$change=='UP'),]
up3 <- diff3[which(diff3$change=='UP'),]
up4 <- diff4[which(diff4$change=='UP'),]

common.up <- intersect(rownames(up1),rownames(up2))
common.up <- intersect(common.up,rownames(up3))
common.up <- data.frame(symbol=intersect(common.up,rownames(up4)))

write.table(common.up,file = 'inter.up.xls',sep = '\t',row.names = F,quote = F)


###down-------
down1 <- diff1[which(diff1$change=='DOWN'),]
down2 <- diff2[which(diff2$change=='DOWN'),]
down3 <- diff3[which(diff3$change=='DOWN'),]
down4 <- diff4[which(diff4$change=='DOWN'),]

common.down <- intersect(rownames(down1),rownames(down2))
common.down <- intersect(common.down,rownames(down3))
common.down <- data.frame(symbol=intersect(common.down,rownames(down4)))

write.table(common.down,file = 'inter.down.xls',sep = '\t',row.names = F,quote = F)


all <- rbind(common.up,common.down)
write.table(all,'inter.all.xls',sep = '\t',row.names = F,quote = F)

mydata1<-list('UP(TCGA)'=rownames(up1),'UP(GSE40275)'=rownames(up2),'UP(GSE21933)'=rownames(up3),'UP(GSE18842)'=rownames(up4))
pdf(file = '01.up.pdf',w=6,h=6)
ggvenn(mydata1,c('UP(TCGA)','UP(GSE40275)','UP(GSE21933)','UP(GSE18842)'),
       fill_color = c("#ffb2b2","#b2e7cb","#F0E442",'#87CEFA'),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb","#F0E442",'#87CEFA'),
       text_color = 'black')
dev.off()
png(file = '01.up.png',w=400,h=400)
ggvenn(mydata1,c('UP(TCGA)','UP(GSE40275)','UP(GSE21933)','UP(GSE18842)'),
       fill_color = c("#ffb2b2","#b2e7cb","#F0E442",'#87CEFA'),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb","#F0E442",'#87CEFA'),
       text_color = 'black')
dev.off()

mydata2<-list('DOWN(TCGA)'=rownames(down1),'DOWN(GSE40275)'=rownames(down2),'DOWN(GSE21933)'=rownames(down3),'DOWN(GSE18842)'=rownames(down4))
pdf(file = '02.DOWN.pdf',w=6,h=6)
ggvenn(mydata2,c('DOWN(TCGA)','DOWN(GSE40275)','DOWN(GSE21933)','DOWN(GSE18842)'),
       fill_color = c("#ffb2b2","#b2e7cb","#F0E442",'#87CEFA'),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb","#F0E442",'#87CEFA'),
       text_color = 'black')
dev.off()
png(file = '02.DOWN.png',w=400,h=400)
ggvenn(mydata2,c('DOWN(TCGA)','DOWN(GSE40275)','DOWN(GSE21933)','DOWN(GSE18842)'),
       fill_color = c("#ffb2b2","#b2e7cb","#F0E442",'#87CEFA'),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb","#F0E442",'#87CEFA'),
       text_color = 'black')
dev.off()
