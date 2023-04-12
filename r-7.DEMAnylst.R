rm(list = ls())
# 01 获取数据集----------
setwd("/data/nas1/luchunlin/project/GY0324-12/")
if (! dir.exists("./07_DEMAnylst")){
  dir.create("./07_DEMAnylst")
}
setwd("./07_DEMAnylst")

diff.meta1 <- read.delim2('../06_DEMs2/DEM_sig(A1vs.A2).xls')
diff.meta2 <- read.delim2('../06_DEMs2/DEM_sig(A2vs.A3).xls')

intersect <- data.frame(DEM=intersect(rownames(diff.meta1),rownames(diff.meta2)))

write.table(intersect,file = 'InterDEM.xls',sep = '\t',row.names = F,quote = F)                        
library(ggvenn)
mydata<-list('DEMs(A1vs.A2)'=rownames(diff.meta1),'DEMs(A2vs.A3)'=rownames(diff.meta2))
pdf('01.venn.pdf',w=5,h=5)
ggvenn(mydata,c('DEMs(A1vs.A2)','DEMs(A2vs.A3)'),
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
ggvenn(mydata,c('DEMs(A1vs.A2)','DEMs(A2vs.A3)'),
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


### 热图圈图---------

##提取表达矩阵
library(lance)
meta.dat <- read.delim2('../00_rawdata/dat.meta.xls')%>%lc.tableToNum()
meta.dat <- meta.dat[,c(1:30)]
colnames(meta.dat)
meta.dat <- meta.dat[intersect$DEM,]%>%t()
# write.table(meta.dat,'DEMs.dat.txt',sep = '\t',row.names = T,quote = F)
library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

heat.dat <- as.matrix(meta.dat)
cir1 <- t(scale(t(heat.dat)))
head(cir1)
#定义热图颜色梯度：

mycol=colorRamp2(c(-2.5, 0.3, 3.1),c("blue","white","red"))#设置legend颜色
mycol1 = colorRamp2(c(-2, 0, 2), c("#003399", "white", "#cccc00"))
#mycol= colorRamp2(c(-1.45, 0, 2.27),c("#0273c2","white","#efc001"))
#mycol= colorRamp2(c(-1.45, 0, 2.27),c("#0da9ce","white","#e74a32"))
#默认参数绘制普通热图
# Heatmap(cir1)
# Heatmap(cir1,row_names_gp = gpar(fontsize = 6),
#         column_names_gp= gpar(fontsize = 8),
#         name="legend")
# 
# #计算数据大小范围
# 
# range(cir1)
# circos.heatmap(cir1,col = mycol)
# circos.clear()
# #在circos.heatmap()中添加参数进行环形热图的调整和美化：
# 
# circos.par(gap.after=c(30)) #circos.par()调整圆环首尾间的距离，数值越大，距离越宽
# circos.heatmap(cir1,col=mycol,
#                
#                dend.side="inside", #聚类放在环形内测,控制行聚类树的方向，inside为显示在圆环内圈，outside为显示在圆环外圈
#                
#                rownames.side="outside", #基因名放在环形外侧,控制矩阵行名的方向,与dend.side相同；但注意二者不能在同一侧，必须一内一外
#                
#                rownames.col="black",
#                
#                rownames.cex=0.5, #字体大小
#                
#                rownames.font=0.5, #字体粗细
#                
#                bg.border = "black",#背景边缘颜色
#                
#                #show.sector.labels = T，#显示分的区域的标签
#                
#                cluster=TRUE) #cluster=TRUE为对行聚类，cluster=FALSE则不显示聚类
# 
# circos.clear()#画完图结束。一定要运行这个，不然后续画图会叠加
# 
# 
# #聚类树的调整和美化：
# 
# # install.packages("dendextend")#改颜色
# 
# # install.packages("dendsort")#聚类树回调
# 
# library(dendextend)
# 
# library(dendsort)
# 
# circos.par(gap.after=c(30)) #circos.par()调整圆环首尾间的距离，数值越大，距离越宽
# 
# circos.heatmap(cir1,col=mycol,dend.side="inside",#dend.side：控制行聚类树的方向，inside为显示在圆环内圈，outside为显示在圆环外圈
#                
#                rownames.side="outside",#rownames.side：控制矩阵行名的方向,与dend.side相同；但注意二者不能在同一侧，必须一内一外
#                
#                track.height = 0.28, #轨道的高度，数值越大圆环越粗
#                
#                rownames.col="black",
#                
#                rownames.cex=1, #字体大小
#                
#                rownames.font=1, #字体粗细
#                
#                cluster=TRUE, #cluster=TRUE为对行聚类，cluster=FALSE则不显示聚类
#                
#                dend.track.height=0.18, #调整行聚类树的高度
#                
#                dend.callback=function(dend,m,si) {#dend.callback：用于聚类树的回调，当需要对聚类树进行重新排序，或者添加颜色时使用包含的三个参数：dend：当前扇区的树状图；m：当前扇区对应的子矩阵；si：当前扇区的名称
#                  
#                  color_branches(dend,k=15,col=1:15) #color_branches():修改聚类树颜色#聚类树颜色改为1，即单色/黑色
#                  
#                }
#                
# )
# 
# 
# #添加图例标签等
# 
# #library(ComplexHeatmap)
# 
# #install.packages("gridBase")
# 
# library(gridBase)
# 
# lg=Legend(title="Exp",col_fun=mycol,
#           
#           direction = c("vertical"),
#           
#           #title_position= c('topcenter')，
#           
# )
# 
# grid.draw(lg)
# 
# #draw(lg, x = unit(0.9,"npc"), y = unit(0.5,"npc"), just = c("right","center"))#画在右边
# 
# #添加列名：
# 
# circos.track(track.index=get.current.track.index(),panel.fun=function(x,y){
#   
#   if(CELL_META$sector.numeric.index==1){   #if(CELL_META$sector.numeric.index == 3) { # the last sector
#     
#     cn=colnames(cir1)
#     
#     n=length(cn)
#     
#     circos.text(rep(CELL_META$cell.xlim[2],n)+convert_x(0.5,"mm"),#x坐标
#                 
#                 1:n+5,#调整y坐标
#                 
#                 cn,cex=0.6,adj=c(0,0.5),facing="inside")}
#   
# },bg.border=NA)
# 
# circos.clear()


#分组热图绘制#

library(dendextend)

#但如果矩阵数据分组，可用split参数来指定分类变量

ann_row = data.frame(pathway=c(rep("A1",10),rep("A2",10),rep("A3",10)))#对行进行注释，用于后续的热图分裂

row.names(ann_row) = rownames(cir1)

ann_row <- as.matrix(ann_row)#在circlize函数中，需要为matrix

#分组绘图

circos.par(gap.after=c(2,2,30)) #circos.par()调整圆环首尾间的距离，数值越大，距离越宽#让分裂的一个口大一点，可以添加行信息

circos.heatmap(cir1,col=mycol,
               
               dend.side="inside",#dend.side：控制行聚类树的方向，inside为显示在圆环内圈，outside为显示在圆环外圈
               
               rownames.side="outside",#rownames.side：控制矩阵行名的方向,与dend.side相同；但注意二者不能在同一侧，必须一内一外
               
               track.height = 0.28, #轨道的高度，数值越大圆环越粗
               
               rownames.col="black",
               
               bg.border="black", #背景边缘颜色
               
               split = ann_row,#用行注释分裂热图
               
               show.sector.labels = T,
               
               rownames.cex=0.9,#字体大小
               
               rownames.font=1,#字体粗细
               
               cluster=TRUE,#cluster=TRUE为对行聚类，cluster=FALSE则不显示聚类
               
               dend.track.height=0.18,#调整行聚类树的高度
               
               dend.callback=function(dend,m,si) {#dend.callback：用于聚类树的回调，当需要对聚类树进行重新排序，或者添加颜色时使用包含的三个参数：dend：当前扇区的树状图；m：当前扇区对应的子矩阵；si：当前扇区的名称
                 
                 color_branches(dend,k=10,col=1:10)#color_branches():修改聚类树颜色#聚类树颜色改为1，即单色/黑色
                 
               }
               
)

#图例与列名设置
library(gridBase)
lg=Legend(title="",col_fun=mycol,direction = c("vertical"))
circle_size = unit(0.05,"snpc")
draw(lg,x = circle_size,just = 'left')

# grid.draw(lg)

circos.track(track.index=get.current.track.index(),panel.fun=function(x,y){
  
  if(CELL_META$sector.numeric.index==3){# the last sector
    
    cn=colnames(cir1)
    
    n=length(cn)
    
    circos.text(rep(CELL_META$cell.xlim[2],n)+convert_x(1,"mm"),#x坐标

                (1:n)*0.22-1.2,#调整y坐标,行距+距离中心距(1:n)*1.2+5,

                cn,cex=0.8,adj=c(0,1.1),facing="inside")

  }
  
},bg.border=NA)

circos.clear()


# #多轨热图绘制#
# 
# #假设有两个热图的矩阵数据（这里仅为一组重复两次以作示范）
# 
# cir1<-t(scale(t(heat.dat)))
# 
# cir2<-t(scale(t(heat.dat)))
# 
# #但如果矩阵数据分组，可用split参数来指定分类变量
# 
# ann_row = data.frame(pathway=c(rep("A1",10),rep("A2",10),rep("A3",10)))#对行进行注释，用于后续的热图分裂
# 
# ann_row2 = data.frame(pathway=c(rep("A1",10),rep("A2",10),rep("A3",10)))#对行进行注释，用于后续的热图分裂
# 
# row.names(ann_row) = rownames(cir1)
# 
# row.names(ann_row2) = rownames(cir2)
# 
# ann_row <- as.matrix(ann_row)#在circlize函数中，需要为matrix
# 
# ann_row2 <- as.matrix(ann_row2)
# 
# #绘图
# 
# circos.par(gap.after=c(2,2,30)) #circos.par()调整圆环首尾间的距离，数值越大，距离越宽#让分裂的一个口大一点，可以添加行信息
# 
# circos.heatmap(cir1,col=mycol,
#                
#                split=ann_row, #用行注释分裂热图
#                
#                rownames.col="black",
#                
#                show.sector.labels = T,
#                
#                #track.height = 0.28, #轨道的高度，数值越大圆环越粗
#                
#                #rownames.side="inside",控制矩阵行名的方向,与dend.side相同；但注意二者不能在同一侧，必须一内一外
#                
#                rownames.cex=0.2,#字体大小
#                
#                rownames.font=1,#字体粗细
#                
#                bg.border="black", #背景边缘颜色
#                
#                dend.side="outside",#dend.side：控制行聚类树的方向，inside为显示在圆环内圈，outside为显示在圆环外圈
#                
#                cluster=TRUE,#cluster=TRUE为对行聚类，cluster=FALSE则不显示聚类
#                
#                dend.track.height=0.2,#调整行聚类树的高度
#                
#                dend.callback=function(dend,m,si) { #dend.callback：用于聚类树的回调，当需要对聚类树进行重新排序，或者添加颜色时使用包含的三个参数：dend：当前扇区的树状图；m：当前扇区对应的子矩阵；si：当前扇区的名称
#                  
#                  color_branches(dend,k=10,col=1:10) #color_branches():修改聚类树颜色#聚类树颜色改为1，即单色/黑色
#                  
#                }
#                
# )
# 
# circos.heatmap(cir2,
#                
#                col = mycol1,
#                
#                split=ann_row2,
#                
#                rownames.side="inside",
#                
#                bg.border="black", #背景边缘颜色
#                
#                rownames.cex=0.3)#加入第二个热图
# 
# #添加列名#
# 
# #第一个环形列名
# 
# circos.track(track.index=get.current.track.index(),panel.fun=function(x,y){
#   
#   if(CELL_META$sector.numeric.index==3){# the last sector
#     
#     cn=colnames(cir1)
#     
#     n=length(cn)
#     
#     circos.text(rep(CELL_META$cell.xlim[2],n)+convert_x(1,"mm"),#x坐标
#                 
#                 (1:n)*0.4+3.6,#调整y坐标,行距+距离中心距(1:n)*1.2+5,
#                 
#                 cn,cex=0.6,adj=c(0,1),facing="inside")
#     
#   }
#   
# },bg.border=NA)
# 
# #第二个环形列名
# 
# circos.track(track.index=get.current.track.index(),panel.fun=function(x,y){
#   
#   if(CELL_META$sector.numeric.index==3){# the last sector
#     
#     cn=colnames(cir2)
#     
#     n=length(cn)
#     
#     circos.text(rep(CELL_META$cell.xlim[2],n)+convert_x(1,"mm"),#x坐标
#                 
#                 (1:n)*0.4+1,#调整y坐标,行距+距离中心距(1:n)*1.2+5,
#                 
#                 cn,cex=0.6,adj=c(0,1),facing="inside")
#     
#   }
#   
# },bg.border=NA)
# 
# #添加放置在左侧的图例#
# 
# # install.packages("gridBase")
# 
# library(gridBase)
# 
# lg_Exp1=Legend(title="Exp1",col_fun=mycol,direction = c("vertical"))
# 
# lg_Exp2=Legend(title="Exp2",col_fun=mycol1,direction = c("vertical"))
# 
# circle_size= unit(0.07,"snpc")
# 
# h= dev.size()
# 
# lgd_list= packLegend(lg_Exp1,lg_Exp2, max_height = unit(2*h,"inch"))
# 
# draw(lgd_list, x = circle_size, just ="left")
# 
# circos.clear()
