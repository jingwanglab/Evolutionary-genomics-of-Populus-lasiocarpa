library(data.table)
library(qvalue)
library(VennDiagram)
library(CMplot)
library(dplyr)
library(ggplot2)
############
for (i in 1:19){
  aim_bio=paste('bio',i,sep="")
  pos=fread('./new/dyy.012.pos',sep='\t',header=F)
  bio=fread(paste('./new/lfmm_out_new/',aim_bio,'_lfmm.env.txt',sep=''),sep='\t')
  bio_adp=fread(paste('./new/asso_loci/lfmm_gemma_5e2/',aim_bio,sep=''),header=F)
  inv_adp=fread('./new_inv/inv_region_95_lfmm_gemma_loci_5e2',header=F)
  id=fread('./new/id.txt',header=F)
  head(inv_adp)
  inv_adp_loc=unique(intersect(unique(inv_adp$V12),bio_adp$V1))
  ma_data=cbind(pos,id,bio[,2])
  head(ma_data)
  colnames(ma_data)=c('CHR','BP','id','pvalue')
  ma_data$TYPE = ifelse(grepl(':',ma_data$id),'1','3')
  ma_data[grepl('INDEL',ma_data$id),"TYPE"]='2'
  
  a=ma_data[ma_data$adp_loc=='adap',]
  ma_data[which(ma_data$id %in% inv_adp_loc),'adp_loc']='adap'
  ma_data[is.na(ma_data$adp_loc),'adp_loc']='netrual'
  ma_data$q=qvalue(ma_data$pvalue)$qvalues
  lfmm=subset(ma_data,q<0.05)
  thre=0.05/6030762
  chr_len <- ma_data %>% 
    group_by(CHR) %>% 
    summarise(chr_len=max(BP))
  ### ????ÿ??chr?ĳ?ʼλ??
  chr_pos <- chr_len  %>%
    mutate(total = cumsum(chr_len) - chr_len) %>%
    dplyr::select(-chr_len)
  ###?????ۼ?SNP??λ??
  Snp_pos <- chr_pos %>%
    left_join(ma_data, ., by="CHR") %>%
    arrange(CHR, BP) %>%
    mutate( BPcum = BP + total)
  X_axis <-  Snp_pos %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  Snp_pos1=subset(Snp_pos,-log10(pvalue)>2)
  #'#E85210',"#F7D092" #3896D6
  #'#'#069EC1','#DCE8BA' jy#DD4D05
  Plot <- ggplot(Snp_pos1, aes(x=BPcum, y=-log10(pvalue))) +
    geom_point(data=subset(Snp_pos1,TYPE=="1"),aes(color=as.factor(CHR)),shape=16,alpha=0.7, size=0.55) +
    geom_point(data=subset(Snp_pos1,TYPE=="2"),aes(color=as.factor(CHR)), shape=15, alpha=0.7,size=1) +
    geom_point(data=subset(Snp_pos1,TYPE=="3"),aes(color=as.factor(CHR)), shape=17,alpha=0.7,size=1)  +
    geom_point(data=subset(Snp_pos1,TYPE=="1"&adp_loc=='adap'),shape=16, size=1,alpha=0.7,color="#DD4D05") +
    geom_point(data=subset(Snp_pos1,TYPE=="2"&adp_loc=='adap'),shape=15, size=1,alpha=0.7,color="#DD4D05") +
    geom_point(data=subset(Snp_pos1,TYPE=="3"&adp_loc=='adap'),shape=17, size=1,alpha=0.7,color="#DD4D05") +
    scale_color_manual(values = rep(c('#069EC1','#DCE8BA' ), 22 )) +
    scale_x_continuous( label = X_axis$CHR, breaks= X_axis$center ) +
    scale_shape_manual(values=c(17,16,19))+
    labs(x="chromosome",y=paste("-log10(P):",aim_bio,sep="")) +
    #ȥ????ͼ????X??֮????gap
    scale_y_continuous(expand = c(0, 0) ) +
    # geom_hline(yintercept = c(min(-log10(lfmm$pvalue)),-log10(thre)), color = c("#5382B1",'grey'), size = 0.8, linetype = c("longdash","longdash")) +
    theme_bw() +
    theme_classic()+
    theme(
      legend.position="none",
      panel.border = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    )
  ggsave(Plot,path='./new_inv/lfmm_plot',filename=paste(aim_bio,'.png',sep=''),width=6,height=1.5,dpi=300)
}
head(Snp_pos)

LG18_INV_data=subset(Snp_pos,CHR=='LG09'&BP>14500000&BP<15300000)
LG18_INV <- ggplot(LG18_INV_data, aes(x=BP, y=-log10(pvalue))) +
  geom_point(data=subset(LG18_INV_data,TYPE=="1"),shape=16, size=1.5,alpha=0.8,color="grey") +
  geom_point(data=subset(LG18_INV_data,TYPE=="2"),shape=15, size=2.5,alpha=0.8,color="#4CB0C8") +
  geom_point(data=subset(LG18_INV_data,TYPE=="3"),shape=17, size=2.5,alpha=0.8,color="#DD5728") +
  #geom_point(data=subset(LG18_INV_data,TYPE=="1"&adp_loc=='adap'),shape=16, size=1.3,alpha=0.7,color="#c64a56") +
  #geom_point(data=subset(LG18_INV_data,TYPE=="2"&adp_loc=='adap'),shape=15, size=1.3,alpha=0.7,color="#c64a56") +
  #geom_point(data=subset(LG18_INV_data,TYPE=="3"&adp_loc=='adap'),shape=17, size=1.3,alpha=0.7,color="#c64a56") +
  scale_shape_manual(values=c(17,16,19))+
  labs(x="chromosome",y=paste("-log10(P):BIO",'4',sep="")) +
  geom_vline(xintercept = c(14863767,14977226),linetype="dashed")+
  geom_hline(yintercept = c(min(-log10(lfmm$pvalue)),-log10(thre)), color = c("#FC7034"), size = 0.5, linetype = c("longdash")) +
  theme_bw() +
  theme_classic()

ggsave(LG18_INV,path='./new_inv/inv/bio4',filename='LG09.pdf',width=6,height=3