.libPaths("F:/360MoveData/Users/dell/Documents/R")
rm(list=ls())
library(ggplot2)
library(cowplot)
setwd("F:/360MoveData/Users/Desktop/IBD&IBE/picture")
geo_dit=read.csv("F:/360MoveData/Users/Desktop/IBD&IBE/neutral/geo_dist.csv",head=F)
env_dist=read.csv("F:/360MoveData/Users/Desktop/IBD&IBE/neutral/scale.csv",head=F)
fst_neutral=read.csv("F:/360MoveData/Users/Desktop/IBD&IBE/neutral/neutral.csv",head=F)
fst_candidate=read.csv("F:/360MoveData/Users/Desktop/IBD&IBE/candidate/candidate.csv",head=F)

trans <- function(raw_data){
  raw_data=raw_data[-1,-1]
  out_data=data.frame(raw_data[,1])
  colnames(out_data)="value"
  for (i in 2:20){
    temp=data.frame(raw_data[i:20,i])
    colnames(temp)="value"
    out_data=rbind(out_data,temp)
  }
  out_data=na.omit(out_data)
  return(out_data)
}

##### all pops #####

plot_data=as.data.frame(matrix(nrow=210,ncol=0))
plot_data$geo_dit=trans(geo_dit)$value
plot_data$env_dist=trans(env_dist)$value
plot_data$fstneutral=trans(fst_neutral)$value
plot_data$fstcandidate=trans(fst_candidate)$value
colnames(plot_data)=c("geo_dit","env_dist","fstneutral","fstcandidate")
plot_data$geo_dit=plot_data$geo_dit/100000
write.csv(plot_data,file="plot_data.csv",quote=F,row.names = F)

plot_data=read.csv("plot_data.csv",head=T)
p1=ggplot(plot_data)+
  geom_point(aes(x=geo_dit,y=fstcandidate),size = 3,alpha=0.7,color="#da6487",shape=21,fill="#f9f2d5")+
  geom_point(aes(x=geo_dit,y=fstneutral),size = 3,alpha=0.7,color="#5d9ac6",shape=21,fill="#d5dde2")+
  geom_smooth(aes(x=geo_dit, y=fstneutral),alpha=0.7,formula = y ~ x, method = lm,se=T,level=0.95,color="#535253", fill="#618cbd",size = 1.5,fullrange = F) +
  geom_smooth(aes(x=geo_dit, y=fstcandidate),alpha=0.7,formula = y ~ x, method = lm,se=T,level=0.95,color="#ca2153",fill="#f5e0ca",fullrange = F, size = 1.5) +
  labs(x = "Geographical Distance (100km)",y = expression(italic(F)[italic(ST)]/(1-italic(F)[italic(ST)])),size = 5.5,colour = "black")+
  panel_border(color = "black", size = 0.7, linetype = 1, remove = FALSE)+
  theme_bw()+
  theme(text=element_text(family="serif"),
        axis.ticks.length = unit(0.25,"lines"),axis.ticks=element_line(colour="black",unit(0.6,"line")),
        axis.text.x=element_text(size=12,colour = "black"),
        axis.text.y=element_text(size=12,colour = "black"), 
        plot.title = element_text(
          size = 15L,
          hjust = 0
        ),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        panel.background=element_rect(fill="white"),
        plot.background = element_rect(fill = "white"),
        #axis.line.x=element_line(colour="black"),
        #axis.line.y=element_line(colour="black"),
        #panel.border=element_blank(),
        panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"mm"))

p2=ggplot(plot_data)+
  geom_point(aes(x=env_dist,y=fstcandidate),size = 3,alpha=0.7,color="#da6487",shape=21,fill="#f9f2d5")+
  geom_point(aes(x=env_dist,y=fstneutral),size = 3,alpha=0.7,color="#5d9ac6",shape=21,fill="#d5dde2")+
  geom_smooth(aes(x=env_dist, y=fstneutral),alpha=0.7,formula = y ~ x, method = lm,se=T,level=0.95,color="#535253", fill="#618cbd",size = 1.5,fullrange = F) +
  geom_smooth(aes(x=env_dist, y=fstcandidate),alpha=0.7,formula = y ~ x, method = lm,se=T,level=0.95,color="#ca2153",fill="#f5e0ca",fullrange = F, size = 1.5) +
  labs(x = "Environment Distance",y = expression(italic(F)[italic(ST)]/(1-italic(F)[italic(ST)])),size = 5.5,colour = "black")+
  panel_border(color = "black", size = 0.7, linetype = 1, remove = FALSE)+
  theme_bw()+
  theme(text=element_text(family="serif"),
        axis.ticks.length = unit(0.25,"lines"),axis.ticks=element_line(colour="black",unit(0.6,"line")),
        axis.text.x=element_text(size=12,colour = "black"),
        axis.text.y=element_text(size=12,colour = "black"), 
        plot.title = element_text(
          size = 15L,
          hjust = 0
        ),
        axis.title.y = element_text(size = 15),
        axis.title.x = element_text(size = 15),
        panel.background=element_rect(fill="white"),
        plot.background = element_rect(fill = "white"),
        #axis.line.x=element_line(colour="black"),
        #axis.line.y=element_line(colour="black"),
        #panel.border=element_blank(),
        panel.grid.major =element_blank(), panel.grid.minor = element_blank(),
        plot.margin=unit(c(0.1,0.1,0.1,0.1),"mm"))
all=plot_grid(p1,p2,align ="v",labels=c("a","b"),label_size = 20,label_fontfamily = "serif",label_fontface = 1,ncol=1)
ggsave(all,file="IBDIBE_allpop.pdf",width=5,height=7.5)