library('dplyr')
library('reshape2')
library(gplots)
library('RColorBrewer')
setwd("./IBD")
sample=read.table("ibd_pop.txt",header=F)
order_list=read.csv('order.csv',header = F)
ibd=fread('plas_ibd')
ibd$length=ibd$V7-ibd$V6
ibd_av=aggregate(ibd$length, by=list(ibd$V1,ibd$V3),sum)
ibd_av$x=log(ibd_av$x)
ibd_av_dup=subset(ibd_av,ibd_av$Group.1!=ibd_av$Group.2)
self=data.frame(Group.1=order_list$V1,Group.2=order_list$V1)
self$x=log(419548847) 
ibd_all=rbind(ibd_av_dup,self)
ind_all=sort(ibd_all$x)

mat=dcast(ibd_all, Group.1~Group.2)
row.names=mat$Group.1
rownames(mat)=mat[,1]
mat=mat[,-1]
mat=as.matrix(mat)
mat[lower.tri(mat)]<-mat[upper.tri(mat)]
mat[lower.tri(mat)]<-t(mat)[lower.tri(mat)]
mat[is.na(mat)] <- 0


colormap=colorRampPalette(rev(brewer.pal(n=10,name="RdBu")))(500)
pdf('ibd_heatmap.pdf',height=20,width=20)
a=heatmap.2(mat,col=colormap,scale=NULL,Colv=NA, Rowv=NA,key=TRUE,symkey=FALSE, density.info="none", trace="none", cexRow=0.3, cexCol=0.3)
dev.off()

write.csv(x=row.names,file='out.csv',row.names = FALSE,col.names = FALSE)
ibd_all=rbind(ibd_av_dup,self)
mat=dcast(ibd_all, Group.1~Group.2)
order1=match(order_list$V1,mat$Group.1)
row.names=order_list$V1
rownames(mat)=mat[,1]
mat=mat[,-1]
mat=as.matrix(mat)
mat[lower.tri(mat)]<-mat[upper.tri(mat)]
mat[lower.tri(mat)]<-t(mat)[lower.tri(mat)]
mat_ordered=mat[,order1]
mat_ordered=mat_ordered[order1,]
mat_ordered[is.na(mat_ordered)] <- 0
pdf('ibd_heatmap_sort.pdf',height=20,width=20)
heatmap.2(mat_ordered,col=colormap,scale=NULL,Colv=NA, Rowv=NA,key=TRUE,symkey=FALSE, density.info="none", trace="none", cexRow=0.3, cexCol=0.3)
dev.off()