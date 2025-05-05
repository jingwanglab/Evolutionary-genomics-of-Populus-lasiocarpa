library(raster)
library(ggplot2)
library(RColorBrewer)
library(reshape2)
library(psych)
library(corrplot)
clim.list <- dir("E:/dyy_env/present/30s", full.names=T) 
clim.layer <-  stack(clim.list) 
sample.coord <-read.table("./new/sample.txt", header=T, stringsAsFactors=F)
crs.wgs <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"  #defines the spatial projection system that the points are in (usually WGS84)
colnames(sample.coord)=c('ID','Longitude','Latitude')
sample.coord.sp <- SpatialPointsDataFrame(sample.coord[,c('Longitude','Latitude')], proj4string=CRS(crs.wgs), data=sample.coord)
clim.points <- raster::extract(clim.layer, sample.coord.sp)  #extracts the data for each point (projection of climate layer and coordinates must match)
clim.points <- cbind(sample.coord, clim.points)  #combines the sample coordinates with the climate data points
ni=raster('E:/dyy_env/Soil_nitrogen_content_5_15cm/nitrogen_5_15.tif')
ni_data<- raster::extract(ni, sample.coord.sp)
clim.points$nitrogen=ni_data
write.table(clim.points,file="./new/culmu_importance/env_gf_input.csv",col.names = T,row.names = F,quote=F)
###
raw_data=clim.points[,5:33]
c=colnames(raw_data)
corrmatrix = cor(raw_data, method = "spearman")
write.csv(corrmatrix  ,file="./new/bio_correlation.csv",row.names = T,quote=