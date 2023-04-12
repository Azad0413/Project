rm(list = ls())
setwd("/data/nas1/luchunlin/project/JNZK-255-2/")
if (! dir.exists("./02_intersect")){
  dir.create("./02_intersect")
}
setwd("./02_intersect")
diff.all <- read.delim2('../01_DEGs/DEG_all.xls',row.names = 1)
library(readxl)
diff <- read.delim2('../01_DEGs/DEG_sig.xls',row.names = 1)
immport <- readxl::read_xlsx('Immport_IRGs.xlsx')%>%dplyr::select('Symbol')
tisidb <- readxl::read_xlsx('TISIDB.xlsx')%>%dplyr::select('gene')
innatedb <- read_xls('innatedb.xls')%>%dplyr::select('Gene Symbol')


colnames(immport) <- 'symbol'
colnames(tisidb) <- 'symbol'
colnames(innatedb) <- 'symbol'
immugene <- rbind(immport,tisidb,innatedb)
length(unique(immugene$symbol))
immugene <- immugene[!duplicated(immugene$symbol),]

metagene <- readxl::read_xlsx('geneset.xlsx')
# metagene <- metagene[c(1:40),]
length(unique(metagene$Symbol))
metagene <- metagene[!duplicated(metagene$Symbol),]
map <- metagene$Symbol[metagene$Symbol%in%rownames(diff)]

intersect <- intersect(rownames(diff),immugene$symbol)%>%as.data.frame()
##741
immumeta <- intersect(immugene$symbol,metagene$Symbol)%>%as.data.frame()
unmap <- metagene$Symbol[!metagene$Symbol%in%immugene$symbol]%>%as.data.frame()
intersect <- data.frame(symbol=intersect(intersect$.,metagene$Symbol))
##6
write.table(intersect,file = 'intersect.xls',sep = '\t',row.names = F,quote = F)
library(ggvenn)
mydata<-list('DEGs'=rownames(diff),'Immune-Related'=immugene$symbol,'Tryptophan metabolism'= metagene$Symbol)
pdf('01.venn.pdf',w=5,h=5)
ggvenn(mydata,c('DEGs','Immune-Related','Tryptophan metabolism'),
       fill_color = c("#ffb2b2","#b2e7cb","#F0E442"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       text_size = 4,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb","#F0E442"),
       text_color = 'black')
dev.off()
png('01.venn.png',w=400,h=400)
ggvenn(mydata,c('DEGs','Immune-Related','Tryptophan metabolism'),
       fill_color = c("#ffb2b2","#b2e7cb","#F0E442"),
       show_percentage = F,
       stroke_alpha = 0.5,
       stroke_size = 0.3,
       text_size = 4,
       stroke_color="white",
       stroke_linetype="solid",
       set_name_color=c("#ffb2b2","#b2e7cb","#F0E442"),
       text_color = 'black')
dev.off()

