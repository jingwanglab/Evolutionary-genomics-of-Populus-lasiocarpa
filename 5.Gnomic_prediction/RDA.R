.libPaths("D:/R/R_library")
rm(list=ls())
library(vegan)
library(ggplot2)
library(data.table)
library(tidyverse)
library(raster)
library("sf")
library(rgdal)
library(maptools)
library(gridExtra)
library(scales)
library(gstat)
library(sp)
library(RColorBrewer)
setwd("C:/Users/SYP/Desktop/fig6/adaptive_index")
source("E:/Riot_trooper/R/vegan/src/adaptive_index.R")
gen=fread("input/LD_frq.txt",head=F)

gen=gen[,-c(1,1554)]
nowenv=read.table("input/PL20_NOW.txt",head=T)
nowenv=nowenv[,c("bio1","bio4","bio5","bio19")]
RDA_outliers <- rda(gen ~ bio1 + bio4 + bio5 + bio19 ,nowenv)

## RDA biplot
TAB_loci <- as.data.frame(scores(RDA_outliers, choices=c(1:2), display="species", scaling="none"))
TAB_var <- as.data.frame(scores(RDA_outliers, choices=c(1:2), display="bp"))
p=ggplot() +
  geom_hline(yintercept=0, linetype="dashed", color = gray(.80), linewidth=0.6) +
  geom_vline(xintercept=0, linetype="dashed", color = gray(.80), linewidth=0.6) +
  geom_point(data = TAB_loci, aes(x=RDA1*3, y=RDA2*3), colour = "#ea845c", size = 2, alpha = 0.8) + #"#F9A242FF"
  geom_segment(data = TAB_var, aes(xend=RDA1, yend=RDA2, x=0, y=0), colour="black", linewidth=0.15, linetype=1, arrow=arrow(length = unit(0.02, "npc"))) +
  geom_text(data = TAB_var, aes(x=1.1*RDA1, y=1.1*RDA2, label = row.names(TAB_var)), size = 2.5, family = "serif") +
  xlab("RDA 1 (60.46%)") + ylab("RDA 2 (16.64%)") +
  facet_wrap(~"Adaptively enriched RDA space") +
  guides(color=guide_legend(title="Locus type")) +
  theme_bw(base_size = 11, base_family = "serif") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=11))
ggsave(p,filename="picture/RDA_biplot.pdf",width=4,height=4)

nowenv <- scale(nowenv, center=TRUE, scale=TRUE)
scale_env <- attr(nowenv, 'scaled:scale')
center_env <- attr(nowenv, 'scaled:center')
range <- sf::st_read("C:/Users/SYP/Desktop/fig6/adaptive_index/shp/New_Shapefile.shp") 
#crs(range) <- '+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs'

ras_now <- stack(list.files("C:/Users/SYP/Desktop/fig6/adaptive_index/input/now_tif/", pattern = ".tif", full.names = T))
#### K=2: 
res_RDA_proj_current <- adaptive_index(RDA = RDA_outliers, K = 2, env_pres = ras_now, range = range, method = "loadings", scale_env = scale_env, center_env = center_env)

## Vectorization of the climatic rasters for ggplot
RDA_proj <- list(res_RDA_proj_current$RDA1, res_RDA_proj_current$RDA2)
RDA_proj <- lapply(RDA_proj, function(x) rasterToPoints(x))
for(i in 1:length(RDA_proj)){
  RDA_proj[[i]][,3] <- (RDA_proj[[i]][,3]-min(RDA_proj[[i]][,3]))/(max(RDA_proj[[i]][,3])-min(RDA_proj[[i]][,3]))
}

## Adaptive genetic turnover projected across lodgepole pine range for RDA1 and RDA2 indexes
TAB_RDA <- as.data.frame(do.call(rbind, RDA_proj[1:2]))
colnames(TAB_RDA)[3] <- "value"
TAB_RDA$variable <- factor(c(rep("RDA1", nrow(RDA_proj[[1]])), rep("RDA2", nrow(RDA_proj[[2]]))), levels = c("RDA1","RDA2"))
pops=read.csv("C:/Users/SYP/Desktop/RONA/PLnow.csv",head=T)
#province=sf::st_read("D:/GIS/省界/bou2_4l.shp")
#### RDA1 only 
TAB_RDA=subset(TAB_RDA,variable=="RDA2")

p2=ggplot(data = TAB_RDA) + 
  geom_raster(aes(x = x, y = y, fill = cut(value, breaks=seq(0, 1, length.out=10), include.lowest = T))) + 
  scale_fill_viridis_d(alpha = 0.8, direction = -1, option = "A", labels = c("Negative scores","","","","Intermediate scores","","","","Positive scores")) +
  scale_x_continuous(limits = c(100.5, 111.8))+
  scale_y_continuous(limits = c(26, 33.4))+
#  geom_sf(fill="transparent",data=province,linewidth=0.5,color="#1E1E1E")+
  coord_sf(xlim = c(100.5, 111.8), ylim = c(26, 33.4), expand = FALSE) +
  geom_point(aes(x=lon, y=lat),pops,fill='white',shape=21,size=3.5)+
  geom_text(aes(x=lon, y=lat,label=POP),pops,size=2.5,color="black")+
  xlab("Longitude") + ylab("Latitude") +
  guides(fill=guide_legend(title="Adaptive index")) +
  facet_grid(~ "Adaptive index") +
  theme_bw(base_size = 11, base_family = "serif") +
  theme(panel.grid = element_blank(), plot.background = element_blank(), panel.background = element_blank(), strip.text = element_text(size=10))
ggsave(p2,file="picture/adaptive_map2.pdf",width=7,height=5,dpi=1000)

writeRaster(res_RDA_proj_current$RDA1,filename="RDA1.tiff",overwrite=TRUE)
writeRaster(res_RDA_proj_current$RDA2,filename="RDA2.tiff",overwrite=TRUE)
##### extract method1 ###### 

#nowenv_data=read.table("input/PL20_NOW.txt",head=T)
#outdata=nowenv_data[,1:3]
#outdata$adaptive_index1 <- data.frame(extract(res_RDA_proj_current$RDA1, nowenv_data[,2:3]))[,1]
#outdata$adaptive_index2 <- data.frame(extract(res_RDA_proj_current$RDA2, nowenv_data[,2:3]))[,1]
#write.csv(outdata,file="adaptive_index.csv",quote=F,row.names = F)

######### method2 ##########
nowenv_data=read.table("input/PL20_NOW.txt",head=T)
outdata=nowenv_data[,1:3]
PL20_coord=as.data.frame(matrix(nrow=20))
PL20_coord$lon= nowenv_data$lon
PL20_coord$lat= nowenv_data$lat 
coordinates(PL20_coord)<-c("lon","lat")

tif1=stack("RDA1.tiff")
tif2=stack("RDA2.tiff")
value_20 <- extract(tif1, PL20_coord, method='simple',proj4string=CRS(crs.wgs),buffer=1500,fun=mean, df=TRUE)
outdata[,"RDA1"]=value_20[2]
value_20 <- extract(tif2, PL20_coord, method='simple',proj4string=CRS(crs.wgs),buffer=1500,fun=mean, df=TRUE)
outdata[,"RDA2"]=value_20[2]
write.csv(outdata,file="adaptive_index.csv",quote=F,row.names = F)

