#! /data/apps/R/4.0.3/bin/Rscript --no-save --no-restore
#conda activate /usr_storage/lzq/software/conda_evns/R
#.libPaths("/w/00/u/user201/R/x86_64-pc-linux-gnu-library/4.0")
library(gradientForest)
library(data.table)
require(raster)
require(geosphere)
require(gdm)
require(foreach)
require(parallel)
require(doParallel)
require(gradientForest)
require(fields)
library(sf)
setwd("/usr_storage2/syp/LZQ/Cand2w/GF")
#offset_outdir="/usr_storage2/syp/LZQ/Cand2w/GF/local_offset"
offset_outdir="/usr_storage2/syp/LZQ/GF_3model/offset"
future_envdir="/usr_storage2/syp/LZQ/GF_3model/future_data"
all_SNPs=read.table("/usr_storage2/syp/LZQ/Cand2w/RONA/cand16336_freq.txt",head=T)
all_SNPs=all_SNPs[,-1]
env=read.table("/usr_storage2/syp/LZQ/env_data/now/PL20_NOW.csv",head=T)
future_file <- commandArgs(trailingOnly = TRUE)
envGF=env[, c("bio1", "bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio10","bio12","bio14","bio15","bio16","bio17","bio18","bio19")]
predNames=c("bio1", "bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio10","bio12","bio14","bio15","bio16","bio17","bio18","bio19")
sites=env[, c("lon", "lat")]
Grid=fread("/usr_storage2/syp/LZQ/env_data/now/PL14k_NOW.csv",header=T)
greengrid=Grid[,c("lon","lat","bio1", "bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio10","bio12","bio14","bio15","bio16","bio17","bio18","bio19")]
grid=Grid[, c("lon", "lat")]

pops=env[,1:3]

preds <- colnames(envGF)
specs <- colnames(all_SNPs)
nSites <- dim(envGF)[1]
nSpecs <- dim(all_SNPs)[2]
maxLevel <- log2(0.368*nrow(envGF)/2)
#all_gfmod <- gradientForest(cbind(envGF, all_SNPs), predictor.vars=colnames(envGF),           response.vars=colnames(all_SNPs), ntree=500, compact=T, nbin =1001,maxLevel=maxLevel, trace=T, corr.threshold=0.5)
#save(all_gfmod,file="all_gfmod.data")
load("/usr_storage2/syp/LZQ/Cand2w/GF/all_gfmod.data")

##### present gf #####
all_tgrid=cbind(greengrid[,c("lon","lat")], predict(all_gfmod,greengrid[,c("bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio10","bio12","bio14","bio15","bio16","bio17","bio18","bio19")]))

#########  reverse offset ##########

future_env_file=paste(future_envdir,"/",future_file,sep="")

for (j in 1:length(future_env_file)){
	cl <- makeCluster(2) #### num of cores  
	registerDoParallel(cl)
	a=unlist(strsplit(future_env_file[j],split=".csv")) 
	temp=gsub("/usr_storage2/syp/LZQ/GF_3model/future_data/","",a[1])
	outfile_offset=paste(offset_outdir,"/",temp,"-PL_reverse_offset.csv",sep="")
	fut_cli=fread(future_env_file, header=T)
	futClimDatGF <- data.frame(fut_cli[,c("lon","lat")],predict(all_gfmod,fut_cli[,c("bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio10","bio12","bio14","bio15","bio16","bio17","bio18","bio19")]))
	reverseOffsetGF <- foreach(i = 1:nrow(futClimDatGF), .packages=c("fields","gdm","geosphere")) %dopar%{
       #get the focal population in future climate
       onePopGF <- futClimDatGF[i,]
       #make prediction between focal population and current climate
       combinedDatGF <- data.frame((all_tgrid[,c("lon","lat")]))
       combinedDatGF["gfOffset"] <- c(rdist(onePopGF[,c("bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio10","bio12","bio14","bio15","bio16","bio17","bio18","bio19")], all_tgrid[,c("bio1","bio2","bio3","bio4","bio5","bio6","bio7","bio8","bio10","bio12","bio14","bio15","bio16","bio17","bio18","bio19")]))
       ##Get metrics for the focal population
       #coordinate of focal population
       coordGF <- onePopGF[,c("lon","lat")]
       #choose the pixels with the minimum offset
       minCoordsGF <- combinedDatGF[which(combinedDatGF$gfOffset == min(combinedDatGF$gfOffset)),]
       #calculate the distance to the sites with minimum fst, and selct the one with the shortest distance
       minCoordsGF["dists"] <- distGeo(p1=coordGF, p2=minCoordsGF[,1:2])
       minCoordsGF <- minCoordsGF[which(minCoordsGF$dists == min(minCoordsGF$dists)),]
       #if multiple sites have the same fst, and same distance, one is randomly chosen
       minCoordsGF <- minCoordsGF[sample(1:nrow(minCoordsGF),1),]
       #get local offset
       offsetGF <- combinedDatGF[which(combinedDatGF$lon == coordGF$lon & combinedDatGF$lat == coordGF$lat),"gfOffset"]
       #get the minimum predicted offset - reverse offset in this case
       minValGF <- minCoordsGF$gfOffset
       #get distance and coordinates of site that minimizes fst
       toGoGF <- minCoordsGF$dists
       minPtGF <- minCoordsGF[,c("lon","lat")]
       #get bearing to the site that minimizes fst
       bearGF <- bearing(coordGF, minPtGF)
       #write out
       outGF <- c(x1=coordGF[[1]], y1=coordGF[[2]],local=offsetGF, reverseOffset=minValGF, predDist=toGoGF, bearing=bearGF, x2=minPtGF[[1]],y2=minPtGF[[2]])
       }
    stopCluster(cl)
    reverseOffsetGF <- do.call(rbind, reverseOffsetGF)
    write.csv(reverseOffsetGF,outfile_offset, row.names=FALSE)
}	

#in this resultant dataframe the columns are:
#x1/y1: focal coordinates
#local: local offset
#reverseOffset: reverse offset
#predDist: distance to site of reverse offset
#bearing: bearing to site of reverse offset
#x2/y2: coordinate of site of reverse offset

