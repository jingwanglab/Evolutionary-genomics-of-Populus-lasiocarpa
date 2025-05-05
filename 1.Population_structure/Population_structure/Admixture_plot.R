#install.packages(c("gridExtra","gtable","label.switching","tidyr"),dependencies=T)
#devtools::install_github('royfrancis/pophelper')
library(gridExtra)
library(gtable)
library(label.switching)
library(tidyr)
library(pophelper)
options(stringsAsFactors = F)
####CV
library("ggplot2")
cv=read.csv("cv.csv",header=FALSE)
colnames(cv)=c("K","CV")
p1=ggplot(cv,aes(x = K, y = CV)) + geom_line() + geom_point(size=3)+
  theme_bw()+scale_x_continuous( labels = as.character(cv$K), breaks = cv$K)
ggsave(p1,filename="CV.pdf",height=5,width=5)

##########################

### INPUT ADMIXTURE RESULT FILES

alist <- readQ(files=list.files(path="./admixture/admixture_out", full.names=T),indlabfromfile=F)
#ԭʼcb_paired=c(4daf4a)
cb_paired=c('#4dbbd5','#fdbf6f','#4daf4a','#ff7f00','#984ea3','#698ed0')

#'#a6611a'��,'#018571'��,'#2166ac'��
inds <- read.delim("plas_0.4.nosex",header=FALSE,stringsAsFactors=F)
twolabset <- read.delim("nj.ind.txt", header=F,stringsAsFactors=F)
if(length(unique(sapply(alist,nrow)))==1) alist <- lapply(alist,"rownames<-",inds$V1)
myFUN<- function(x) {
    x[twolabset$V1, ]
  }
alist<-lapply(alist,myFUN)
plotQ(alignK(alist[2:3]),imgoutput="join",barsize=1,height=2,width=15, clustercol= cb_paired,
      splab=paste0("K=",sapply(alist[2:3],ncol)),showlegend=T,
      showindlab=T,useindlab=T,indlabsize=2,indlabspacer=0,sharedindlab=T,
      barbordersize=0,outputfilename="plas_admixture",imgtype="pdf",exportpath="./admixture")
#barbordercolour="white"
