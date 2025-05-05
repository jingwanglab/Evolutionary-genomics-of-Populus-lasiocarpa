#!/bin/R
library(lostruct)
library(colorspace)
library(jsonlite)
library(RColorBrewer)
opt <- fromJSON("/usr_storage/lzq/work/dyy/localpca/lostruct_results/type_snp_size_250_weights_none_jobid_653864/config.json")
chroms <- opt$chrom.names
bcf.files <- opt$bcf_files
names(bcf.files) <- chroms
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

setwd("/usr_storage/lzq/work/dyy/localpca/lostruct_results/type_snp_size_250_weights_none_jobid_653864/")
pca.files <- list.files(".",".*.pca.csv")
mds.file <- "mds_coords_20.csv"
regions.files <- list.files(".",".*.regions.csv")
mds.coords <- read.csv(mds.file,header=TRUE, stringsAsFactors=TRUE)
regions <- do.call( rbind, lapply( regions.files, read.csv, header=TRUE, stringsAsFactors=TRUE ) )
mds.cols <- (1:ncol(mds.coords))[-(1:2)]
chroms <- unique(regions$chrom)
chrom.starts <- tapply( regions$start, regions$chrom, min, na.rm=TRUE )
chrom.ends <- tapply( regions$end, regions$chrom, max, na.rm=TRUE )
chrom.spacing <- floor(.05*mean(chrom.ends))
chrom.offsets <- c(0,cumsum(chrom.spacing+chrom.ends))
names(chrom.offsets) <- c(chroms,"end")
chrom.dividers <- c(0,chrom.offsets[-1])-chrom.spacing/2
chrom.mids <- chrom.dividers[-1] - diff(chrom.dividers)/2
regions$pos <- chrom.offsets[match(regions$chrom,chroms)]+(regions$start+regions$end)/2
chrom.cols <- rainbow_hcl(length(chroms), c=90, end=.9*360)[as.numeric(regions$chrom)]
chrom.plot <- function (y,ylab='',main='',chrom.labels=TRUE,...) {
    plot(0, type='n', xlim=range(chrom.offsets/1e6), ylim=range(y,finite=TRUE),
         #xlab='',
         xaxt='n',  ylab=ylab, main=main)
    if (length(chroms)>1) for (k in 1:floor(length(chroms)/2)) {
        rect( xleft=chrom.dividers[2*k-1]/1e6, xright=chrom.dividers[2*k]/1e6,
             ybottom=par("usr")[3], ytop=par("usr")[4],
             border=NA, col=adjustcolor("grey",0.25) )
    }
    abline( v=chrom.dividers/1e6, lty=3, col=adjustcolor("grey",0.5), lwd=2 )
    if (chrom.labels) axis( 1, at=chrom.mids/1e6, labels=paste("chromosome", chroms), las=0, tick=FALSE )
    points( regions$pos/1e6, y,cex=1.2, ...)
}
i=3
top0.01=quantile(mds.coords[,i],0.995,na.rm=TRUE)
low0.01=quantile(mds.coords[,i],0.005,na.rm=TRUE)
ccols <- rep("black",nrow(mds.coords))
top=na.omit(regions[mds.coords[,i] > top0.01,])
top=top[top[,3]-top[,2]<20000,]
low=na.omit(regions[mds.coords[,i] < low0.01,])
low=low[low[,3]-low[,2]<20000,]
ccols[as.numeric(row.names(top))]<-'#CB1F5A'
ccols[as.numeric(row.names(low))]<-'#36CACC'

pdf("mds1.pdf",width=12,height=3.2)
chrom.plot( mds.coords[,i], pch=20,
            xlab= "Position (Mb)" , # main=paste("MDS coordinate",match(k,mds.cols)),
            ylab=colnames(mds.coords)[i],
            col=adjustcolor(ccols,0.75) )
dev.off()
png("mds1.png",width=1200,height=320)
chrom.plot( mds.coords[,i], pch=20,
            xlab= "Position (Mb)" , # main=paste("MDS coordinate",match(k,mds.cols)),
            ylab=colnames(mds.coords)[i],
            col=adjustcolor(ccols,0.75) )
dev.off()
