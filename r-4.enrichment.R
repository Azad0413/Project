rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/BJTC-420-1/")
if (! dir.exists("./04_enrichment")){
  dir.create("./04_enrichment")
}
setwd("./04_enrichment")

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(enrichplot)
library(GOplot)
intersect <- read.delim2('../03_intersect/intersect.xls')


library(enrichR)
dbs <- listEnrichrDbs()

##PATHWAY--------
dbs$libraryName
dbs <- c("WikiPathway_2021_Human","KEGG_2021_Human",
         "Reactome_2022","BioPlanet_2019")

symbol <- intersect$symbol
enrichr <- enrichr(symbol,dbs)

result1 <- data.frame(enrichr$WikiPathway_2021_Human)
result2 <- data.frame(enrichr$KEGG_2021_Human)
result3 <- data.frame(enrichr$Reactome_2022)
result4 <- data.frame(enrichr$BioPlanet_2019)
result <- rbind(result1,result2,result3,result4)
result$Method <- c(rep('WikiPathway',213),rep('KEGG',150),rep('Reactome',371),rep('BioPlanet',449))
write.table(result,file = 'Pathway.xls',sep = '\t',row.names = F,quote = F)

trace(plotEnrich,edit = T)
p1 <- plotEnrich(enrichr$WikiPathway_2021_Human,showTerms = 8,numChar = 100,y='Count',orderBy = 'P.value',
                 title = 'WikiPathway',xlab = '')+scale_fill_continuous(low = "#FFDEAD", high = "#FF8C00")

p1

p2 <- plotEnrich(enrichr$KEGG_2021_Human,showTerms = 8,numChar = 80,y='Count',orderBy = 'P.value',
                 title = 'KEGG Pathway',xlab = '')+scale_fill_continuous(low = "#A4D3EE", high = "#7D9EC0")

p2

p3 <- plotEnrich(enrichr$Reactome_2022,showTerms = 8,numChar = 100,y='Count',orderBy = 'P.value',
                 title = 'Reactome',xlab = '')+scale_fill_continuous(low = "#90EE90", high = "#2E8B57")

p3

p4 <- plotEnrich(enrichr$BioPlanet_2019,showTerms = 8,numChar = 80,y='Count',orderBy = 'P.value',
                 title = 'BioPlanet',xlab = '')+scale_fill_continuous(low = "#DDA0DD", high = "#800080")

p4

library(patchwork)
p5 <- p1+p2+p3+p4+plot_layout(ncol = 1)
p5
ggsave(p5,filename = '02.Pathway.pdf',w=10,h=9)
ggsave(p5,filename = '02.Pathway.png',w=10,h=9)



##GO--------
dbs$libraryName
dbs <- c("GO_Biological_Process_2021","GO_Cellular_Component_2021",
         "GO_Molecular_Function_2021")

enrichGO <- enrichr(symbol,dbs)

resultbp <- data.frame(enrichGO$GO_Biological_Process_2021)
write.table(resultbp,file = 'GOBP.xls',sep = '\t',row.names = F,quote = F)

resultcc <- data.frame(enrichGO$GO_Cellular_Component_2021)
write.table(resultcc,file = 'GOCC.xls',sep = '\t',row.names = F,quote = F)

resultmf <- data.frame(enrichGO$GO_Molecular_Function_2021)
write.table(resultmf,file = 'GOMF.xls',sep = '\t',row.names = F,quote = F)

trace(plotEnrich,edit = T)
p6 <- plotEnrich(enrichGO$GO_Biological_Process_2021,showTerms = 8,numChar = 100,y='Count',orderBy = 'P.value',
                 title = 'Biological Process',xlab = '')+scale_fill_continuous(low = "#FFDEAD", high = "#FF8C00")

p6

p7 <- plotEnrich(enrichGO$GO_Cellular_Component_2021,showTerms = 8,numChar = 80,y='Count',orderBy = 'P.value',
                 title = 'Cellular Component',xlab = '')+scale_fill_continuous(low = "#A4D3EE", high = "#7D9EC0")

p7

p8 <- plotEnrich(enrichGO$GO_Molecular_Function_2021,showTerms = 8,numChar = 110,y='Count',orderBy = 'P.value',
                 title = 'Molecular Function',xlab = '')+scale_fill_continuous(low = "#90EE90", high = "#2E8B57")

p8


library(patchwork)
p9 <- p6+p7+p8+plot_layout(ncol = 1)
p9
ggsave(p9,filename = '01.GObar.pdf',w=10,h=9)
ggsave(p9,filename = '01.GObar.png',w=10,h=9)
