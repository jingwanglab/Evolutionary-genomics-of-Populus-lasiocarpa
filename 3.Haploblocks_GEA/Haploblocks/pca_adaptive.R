#!/bin/R
library(lostruct)
library(colorspace)
library(jsonlite)
library(RColorBrewer)
library(dplyr)
library(raster)
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scatterpie) 
library(ggpubr) 

creategroup <- function(tiff){
  colnames(tiff)=c("x","y","bio_value")
  var=max(tiff$bio_value)-min(tiff$bio_value)
  min=floor(min(tiff$bio_value))
  if (min>0){
    group10=paste("T91:",min+floor(var/10*9),'+',sep='')
    group9=paste("T9:",min+floor(var/10*8),'~',min+floor(var/10*9),sep='')
    group8=paste("T8:",min+floor(var/10*7),'~',min+floor(var/10*8),sep='')
    group7=paste("T7:",min+floor(var/10*6),'~',min+floor(var/10*7),sep='')
    group6=paste("T6:",min+floor(var/10*5),'~',min+floor(var/10*6),sep='')
    group5=paste("T5:",min+floor(var/10*4),'~',min+floor(var/10*5),sep='')
    group4=paste("T4:",min+floor(var/10*3),'~',min+floor(var/10*4),sep='')
    group3=paste("T3:",min+floor(var/10*2),'~',min+floor(var/10*3),sep='')
    group2=paste("T2:",min+floor(var/10*1),'~',min+floor(var/10*2),sep='')
    group1=paste("T1:",0,'-',min+floor(var/10*1),sep='')  
    
  }else {
    group10=paste("T91:",min+floor(var/10*9),'+',sep='')
    group9=paste("T9:",min+floor(var/10*8),'~',min+floor(var/10*9),sep='')
    group8=paste("T8:",min+floor(var/10*7),'~',min+floor(var/10*8),sep='')
    group7=paste("T7:",min+floor(var/10*6),'~',min+floor(var/10*7),sep='')
    group6=paste("T6:",min+floor(var/10*5),'~',min+floor(var/10*6),sep='')
    group5=paste("T5:",min+floor(var/10*4),'~',min+floor(var/10*5),sep='')
    group4=paste("T4:",min+floor(var/10*3),'~',min+floor(var/10*4),sep='')
    group3=paste("T3:",min+floor(var/10*2),'~',min+floor(var/10*3),sep='')
    group2=paste("T2:",min+floor(var/10*1),'~',min+floor(var/10*2),sep='')
    group1=paste("T1:",min,'~',min+floor(var/10*1),sep='')
  }
  tiff$level=ifelse(tiff$bio_value>=min+floor(var/10*9),group10,
                    ifelse(tiff$bio_value>=min+floor(var/10*8),group9,
                           ifelse(tiff$bio_value>=min+floor(var/10*7),group8,
                                  ifelse(tiff$bio_value>=min+floor(var/10*6),group7,
                                         ifelse(tiff$bio_value>=min+floor(var/10*5),group6,
                                                ifelse(tiff$bio_value>=min+floor(var/10*4),group5,
                                                       ifelse(tiff$bio_value>=min+floor(var/10*3),group4,
                                                              ifelse(tiff$bio_value>=min+floor(var/10*2),group3,
                                                                     ifelse(tiff$bio_value>=min+floor(var/10*1),group2,group1
                                                                     )))))))))
  
  
  tiff$level<-factor(tiff$level,levels=c(as.character(group1),as.character(group2),as.character(group3),as.character(group4),as.character(group5),as.character(group6),as.character(group7),as.character(group8),as.character(group9),as.character(group10)))
  return(tiff)}
crs.wgs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"  #defines the spatial projection system that the 
mask= read_sf("New_Shapefile.shp")%>%st_transform(, crs.wgs)
jiangyu=c("#E5E8C5","#EFF7B5","#CEECB2","#96D6B9","#5EC0C0","#2DA4C1","#2080B7","#2254A3","#0A3D87","#0E56A5")

running_cov <- function (f, n, normalize.rows=TRUE) {
    if (is.numeric(n) && length(n)==1) { n <- seq_len(n) }
    x <- f(n[1])
    z <- !is.na(x)
    colm <- colMeans(x,na.rm=TRUE)
    colm[is.na(colm)]<-0
    x <- sweep(x,2,colm,"-")
    x[!z] <- 0
    nn <- crossprod(z)       # matrix of number of shared nonmissings
    sums <- crossprod(x,z)       # matrix of conditional sums: sums[i,j] = sum( x[,i] * !is.na(x[,j]) )
    sumsq <- crossprod(x)        # matrix of sum of products of shared nonmissings
    for (k in n[-1]) {
        x <- sweep( f(k), 2, colm, "-" )
        z <- !is.na(x)
        x[!z] <- 0
        nn <- nn + crossprod(z)
        sums <- sums + crossprod(x,z)
        sumsq <- sumsq + crossprod(x)
    }
    out <- ( (1/(nn-1))*sumsq - sums*t(sums)/(nn*(nn-1)) )
    if (normalize.rows){
        Imat <- diag(ncol(x)) - 1/ncol(x)
        out <- Imat %*% out %*% Imat
    }
    return(out)
}

pca_plot<-function(opt,reg,y){
bcf.files <- opt$bcf_files
chroms<-opt$chrom.names
names(bcf.files) <- chroms
corner.npc<-2
sample.ids <- vcf_samples(bcf.files[1])
if (!is.null(opt$sample_info)) {
    samp.file <- opt$sample_info
    samps <- read.table(samp.file,sep="\t",header=TRUE, stringsAsFactors=TRUE)
    names(samps) <- tolower(names(samps))
    # hack for msprime output
    samps <- droplevels( samps[match(sample.ids,samps$id),] )
    samps$population <- factor(samps$population)
} else {
    warning("No population information in the sample file, %s.", samp.file)
    samps <- data.frame( 
            ID=sample.ids,
            population=factor(rep("pop",length(sample.ids))) )
}
pop.pch <- seq_len(nlevels(samps$population))
pop.cols <- rainbow_hcl(nlevels(samps$population))
qfun <- multi_vcf_query_fn( chrom.list=chroms, file=bcf.files, regions=reg )
corner.covmats <- running_cov(qfun,1:nrow(reg), normalize.rows=TRUE)
corner.pca <- cov_pca(covmat=corner.covmats, k=corner.npc, w=opt$weights)
vectors <- matrix( corner.pca[-(1:(1+corner.npc))], ncol=corner.npc )
colnames(vectors) <- paste("PC", c(1,2),c(corner.pca[corner.npc],corner.pca[corner.npc+1]))
pca=data.frame(sample=samps$id,pop=samps$population,PC1=vectors[,1],PC2=vectors[,2])
try_3_clusters <-try(kmeans(pca[,3], 3, centers=c(min(pca[,3]),(min(pca[,3])+max(pca[,3]))/2,max(pca[,3]))))
if("try-error" %in% class(try_3_clusters)){kmeans_cluster <-try(kmeans(pca[,], 2, centers=c(min(pca[,3]),max(pca[,3]))))}else{kmeans_cluster <-kmeans(pca[,3], 3, centers=c(min(pca[,3]),(min(pca[,3])+max(pca[,3]))/2,max(pca[,3])))}
pca$cluster <- kmeans_cluster$cluster - 1
betweenss <- kmeans_cluster$betweenss
pca$ID=gsub("-[^-]+$", "",pca$sample)
if (betweenss>0.95){
pca_plot<-ggplot(pca,aes(x=PC1,y=PC2)) + geom_point(aes(shape=as.factor(cluster),color=as.factor(pop)),size=2) +theme_bw() + ggtitle(paste(" BetweenSS = ",betweenss, sep="")) 
out_name=paste(chrom,":",reg$start,"-",reg$end,sep='')
out_name_het=paste(chrom,reg$start,reg$end,sep='_')
command=paste('bcftools view  --regions ',out_name," ",chrom,"/",chrom,'.DP.GQ.bialleles.bed.mis0.8.maf0.05.180ind.recode.vcf.gz ','|perl vcf2het.pl>/usr_storage/lzq/work/dyy/localpca/INV/het/',out_name_het,sep='')
system(command)
het=read.csv(paste('./het/',out_name_het,sep=''),sep='\t')
plot_het<-pca %>% inner_join(.,het, by = "sample")%>%ggplot(.,aes(x=as.factor(cluster),y=percent_het,color=cluster))+geom_boxplot(aes(x=as.factor(cluster))) +theme_bw() 
head(het)
###########################
location=read.csv('env_gf_input.csv',sep=' ')
cal_freq<-function(x){
  temp=pca[pca$ID==x,]$cluster
  y=sum(as.numeric(temp))/(length(temp)*2)
  return(y)}
location['INV']=sapply(location$ID,cal_freq)
location['no_INV']=1-location$INV
relation_bio=vector()
for (j in 4:31){
    cor=cor.test(location$INV,location[,j],method='spearman',exact=F)
    p=cor$p.value
    roh=cor$estimate
#    if(p<0.01){relation_bio=append(relation_bio,colnames(location)[j])}}
    relation_bio=append(relation_bio,colnames(location)[j])}
  if (length(relation_bio)>1){
  aim_bio=sample(relation_bio,1)
	print(aim_bio)
  bio_raster=raster(paste("/usr_storage/lzq/data/30s/",aim_bio,'.asc',sep=''))
  bio_raster<-mask(bio_raster,mask)
  resdf<-as.data.frame(rasterToPoints(bio_raster),xy=TRUE)%>%drop_na()%>%creategroup()
  head(resdf)
  env<-ggplot()+geom_raster(aes(x=x,y=y,fill=factor(level)),data=resdf,stat = "identity")+scale_fill_manual(values=c("#A242DB","#B8D346",jiangyu))+geom_scatterpie(aes(x=Longitude, y=Latitude), data=location,cols=colnames(location)[33:34],alpha=0.8,size=0.1)+theme_test()
  data=location[,c('ID',aim_bio,'INV')]%>%inner_join(.,pca[,c('ID','pop')])
  colnames(data)=c('ID','BIO','INV','pop')
  rel<-ggplot(data,aes(x=INV,y=BIO))+geom_point(aes(color=pop),size=3)+geom_smooth(method="lm",se=TRUE)+ggtitle(aim_bio)+theme_test()
  p1=ggarrange(pca_plot,plot_het,env,rel,nrow=2,ncol=2)
  ggsave(file=paste('./inv_pca_plot/',y,chrom,reg$start,reg$end,'.png',sep='_'),p1,width=8,height=4.5)
  ggsave(file=paste('./inv_pca_plot/',y,chrom,reg$start,reg$end,'.pdf',sep='_'),p1,width=8,height=4.5)
  write.csv(pca,file=paste('./inv_pca/',y,chrom,reg$start,reg$end,".csv",sep='_'))}
}
}
