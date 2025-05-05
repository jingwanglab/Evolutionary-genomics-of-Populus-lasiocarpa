#!/bin/R
library(lostruct)
library(colorspace)
library(jsonlite)
library(RColorBrewer)
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
corner.npc <- 2
chroms <- unique(regions$chrom)
pop.names <- levels(samps$population)
pop.cols <- rainbow_hcl(nlevels(samps$population))
pop.pch <- seq_len(nlevels(samps$population))
i=3
	top0.01=quantile(mds.coords[,i],0.995,na.rm=TRUE)
	low0.01=quantile(mds.coords[,i],0.005,na.rm=TRUE)
	corner.regions <- vector(mode='list', length=3)
	 #print(levels(mds.coords[mds.coords[,i]> top0.01,1]))
	 #print(levels(mds.coords[mds.coords[,i]< low0.01,1]))
	corner.regions[[1]]=na.omit(regions[mds.coords[,i] > top0.01,])
#	corner.regions[[1]]=corner.regions[[1]][(corner.regions[[1]][,3]-corner.regions[[1]][,2]<20000),]
	corner.regions[[2]]=na.omit(regions[mds.coords[,i] < low0.01,])
#	corner.regions[[2]]=corner.regions[[2]][(corner.regions[[2]][,3]-corner.regions[[2]][,2]<20000),]
	temp=regions[-unique(c(as.numeric(row.names(corner.regions[[1]])),as.numeric(row.names(corner.regions[[2]])))),]
	temp=na.omit(temp)
	write.table(temp,"/usr_storage/lzq/work/dyy/localpca/background_0.005",row.names=F,sep="\t",quote=F)
	temp=temp[sample(1:nrow(temp),500,replace=T),]
	write.table(corner.regions[[2]],"/usr_storage/lzq/work/dyy/localpca/divergence_0.005",row.names=F,sep="\t",quote=F)
	write.table(corner.regions[[1]],"/usr_storage/lzq/work/dyy/localpca/balance_0.005",row.names=F,sep="\t",quote=F)
	corner.regions[[3]]=temp
	corner.covmats <- vector(mode='list', length=3)
	for (k in seq_along(corner.covmats)) {
	    reg <- corner.regions[[k]]
	    qfun <- multi_vcf_query_fn( chrom.list=chroms, file=bcf.files, regions=reg )
           corner.covmats[[k]] <- running_cov(qfun,1:nrow(reg), normalize.rows=TRUE)}
       corner.pca <- vector(mode='list', length=3)
	for (k in seq_along(corner.pca)) {
	    	corner.pca[[k]] <- cov_pca(covmat=corner.covmats[[k]], k=corner.npc, w=opt$weights)}
	pdf("mds1_bg_0.005.pdf",width=10,height=3.7) 
	layout(t(1:3))
	for (k in 1:3) {
            vectors <- matrix( corner.pca[[k]][-(1:(1+corner.npc))], ncol=corner.npc )
            colnames(vectors) <- paste("PC", c(1,2),c(corner.pca[[k]][corner.npc],corner.pca[[k]][corner.npc+1]))
            par(mgp=c(0.7,0.7,0), mar=c(2,2,2,0)+.1)
            plot(vectors, pch=pop.pch[samps$population], 
                    col=pop.cols[samps$population],
                    xaxt='n', yaxt='n' )}
	dev.off()


	

