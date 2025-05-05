ary(gradientForest)
setwd('D:/?./new/culmu_importance/')
env_gf<-read.csv("./new/culmu_importance/env_gf_input.csv",sep=" ")
######### 
library(data.table)
library(vegan)
library(rgdal)
library(dplyr)
candidates=dir("./new/culmu_importance/bio_frq")
reference=fread("reference_frq",  stringsAsFactors=T,sep="\t",fill=TRUE,header = F)
reference=as.data.frame(reference)
del=c()
for (i in seq_along(reference)){
  if(length(unique(reference[,i])) <= 5){del=append(del,i)}
}
reference=reference[,-del]
reference=reference[,-1]


for (file in candidates) {
  file='bio15'
  bio=file
  file_path=paste(".o_frq/",file,sep="")
  if (file.info(file_path)$size == 0){next}
  else {
  candidate=fread(file_path,  stringsAsFactors=T,sep="\t",fill=TRUE,header = F)
  num=as.numeric(dim(candidate)[2])+50
  reference1=reference[,sample(1:nrow(reference),num,replace=T)]
  candidate=as.data.frame(candidate)
  del=c()
  for (i in seq_along(candidate)){
    if(length(unique(candidate[,i])) <= 5){del=append(del,i)}}
  if (length(del)>=1){candidate=candidate[,-del]}

  maxLevel <- log2(0.368*nrow(env_gf)/2)
  env_gf1=as.data.frame(env_gf[,bio])
  colnames(env_gf1)=bio
  gf_candidate <- gradientForest(cbind(env_gf1, candidate),
                               predictor.vars= bio,
                               response.vars=colnames(candidate),
                               ntree=400,maxLevel=maxLevel, trace=T, corr.threshold=0.40)
  gf_reference <- gradientForest(cbind(env_gf1, reference1),
                               predictor.vars= bio,
                               response.vars=colnames(reference1),
                               ntree=400,maxLevel=maxLevel, trace=T, corr.threshold=0.40)
######################################
#bio_cand <- gf_candidate$overall.imp[order(gf_candidate$overall.imp,decreasing = T)]
#pdf("GF_candidate_VariableImportance2.pdf")
#plot(gf_candidate, plot.type = "O")
#dev.off()

#barplot(bio_cand,las=2,cex.names=0.8,col=rep('grey',8),ylab="Weigthed importance (R-sqr)")
  most_cand=bio
  pop_turn <- predict(gf_candidate,env_gf1)
  temp <- data.frame(bio=env_gf[,most_cand],imp=pop_turn[,most_cand]) 
  west <- c(1,2,14:19)
  east <- c(3:13,20)
  categories <- list(east=rownames(pop_turn)[east],west=rownames(pop_turn)[west]) 
  temp_cand_overall <- cumimp(gf_candidate,predictor= most_cand,
                              type=c("Overall"),standardize = T)
  temp_refe_overall <- cumimp(gf_reference,predictor= most_cand,
                              type=c("Overall"),standardize = T)
  pdf(paste('imcu_',most_cand,".pdf",sep=""),height = 5,width=5)
  plot(temp_cand_overall,type="n",
       ylab="Cumulative importance",xlab= most_cand)
  library(MaizePal)
  lines(temp_cand_overall,col=MaizePal::maize_pal("RubyGold")[2],lwd=4,lty=1)
  lines(temp_refe_overall,col=MaizePal::maize_pal("RubyGold")[1],lwd=4,lty=1)
   west_col=rep('#82BCCC',8)
  east_col=rep('#E8B26D',12)
  id_c <- order(temp$bio[east])
  id_ccol <- as.character(cut(1:length(id_c),length(east_col),labels=east_col))
  id_w <- order(temp$bio[west])
  id_wcol <- as.character(cut(1:length(id_w),length(west_col),labels=west_col))
  points(temp$bio[west][id_w],temp$imp[west][id_w],pch=21,bg=rev(id_wcol),cex=1)
  points(temp$bio[east][id_c],temp$imp[east][id_c],pch=21,bg=id_ccol,cex=1)
  dev.off()}}
