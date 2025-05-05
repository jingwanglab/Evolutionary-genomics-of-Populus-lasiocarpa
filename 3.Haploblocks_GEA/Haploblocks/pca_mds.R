dirs=list.files("/usr_storage/lzq/work/dyy/localpca/INV/lostruct_results/",full.names=T)
setwd('/usr_storage/lzq/work/dyy/localpca/INV/')
library(ggpubr)
library(jsonlite)
source('./pca_adaptive.R')
consecutive_windows<-function(x,y){
all_outlier_data=data.frame()
outlier_data=data.frame()
for (i in 1:dim(x)[1]){
  if (x[i,y]>1.5){
    outlier_data=rbind(outlier_data,x[i,])}
  else{
  if (dim(outlier_data)[1]>9){
    reg <- data.frame(chrom=chrom,start=outlier_data[1,2],end=outlier_data[dim(outlier_data)[1],3])
    pca_plot(opt,reg,y)
    write.csv(reg,file=paste('./inv_region/',y,chrom,reg$start,reg$end,".csv",sep='_'))
    all_outlier_data=rbind(all_outlier_data,outlier_data)}
    outlier_data=data.frame()}
  }
return(row.names(all_outlier_data))
}

for (dir in dirs){
  print(dir)
  opt <- fromJSON(paste(dir,"/config.json",sep=""))
  chrom <- opt$chrom.names
  print(chrom)
  region=read.csv(paste(dir,'/',chrom,'.regions.csv',sep=""))
  mds=read.csv(paste(dir,'/mds_coords.csv',sep=""))
  mds=mds[,c(3,4,5)]
  data=cbind(region,mds)
  data$zscore1=abs(data$MDS1-mean(na.omit(data$MDS1)))/sd(na.omit(data$MDS1))
  data[is.na(data$zscore1),'zscore1']<-0
  data$zscore2=abs(data$MDS2-mean(na.omit(data$MDS2)))/sd(na.omit(data$MDS2))
  data[is.na(data$zscore2),'zscore2']<-0
  data$zscore3=abs(data$MDS3-mean(na.omit(data$MDS3)))/sd(na.omit(data$MDS3))
  data[is.na(data$zscore3),'zscore3']<-0
  data[consecutive_windows(data,'zscore1'),'mds1_type']='sv'
  data[is.na(data$mds1_type),'mds1_type']='sv_no'
  data[consecutive_windows(data,'zscore2'),'mds2_type']='sv'
  data[is.na(data$mds2_type),'mds2_type']='sv_no'
  data[consecutive_windows(data,'zscore3'),'mds3_type']='sv'
  data[is.na(data$mds3_type),'mds3_type']='sv_no'
  p1=ggplot(data,aes(x=start,y=MDS1,color=mds1_type))+geom_point(size=1.5)+theme_test()
  p2=ggplot(data,aes(x=start,y=MDS2,color=mds2_type))+geom_point(size=1.5)+theme_test()
  p3=ggplot(data,aes(x=start,y=MDS3,color=mds3_type))+geom_point(size=1.5)+theme_test()
  p4=ggarrange(p1,p2,p3,nrow=3,ncol=1)
  ggsave(file=paste('./mds/',chrom,'.png',sep=''),p4,width=12,height=9)
  ggsave(file=paste('./mds/',chrom,'.pdf',sep=''),p4,width=12,height=9)
}
 warnings()
