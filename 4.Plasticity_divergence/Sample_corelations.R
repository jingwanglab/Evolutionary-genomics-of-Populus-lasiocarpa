setwd('/Users/zhiqinlong/Desktop/computer/work/lasiocarpa/expressssionnew')
library(pheatmap)
es_sc=read.csv("es_sc_all_gene_tpm_filt_lowexp.csv",row.names='gene_id')
es_sc_r = cor(es_sc, method = "spearman")
colnames(es_sc)=gsub('es','E',colnames(es_sc))
colnames(es_sc)=gsub('sc','W',colnames(es_sc))

library(corrplot)
pdf("sc_es_sample_cor.pdf",width=9,height=9)
color=colorRampPalette(rev(c("#EA281F","#F79B69","#F5ECAF",'#DCEEEB',"#2395DB")))
a=pheatmap(es_sc_r,scale='none',clustering_method ='complete',clustering_distance_rows = "euclidean",cluster_cols= T,show_rownames =F,color = color(50))
#corrplot(es_sc_r,col =color(20),method = "color",is.corr=F,tl.col = "black")
dev.off()
library(dplyr)
es_sc_log <- log2(es_sc+1)
colnames(es_sc)
pca1=prcomp(t(es_sc_log),center = T, scale. = F)


da=data.frame(sample=rownames(pca1$x),color = c(rep("ck",3),rep("hs",3),rep("ws",3),
                                                rep("ck",3),rep("hs",3),rep("ws",3),
                                                rep("ck",2),rep("hs",3),rep("ws",3),
                                                rep("ck",3),rep("hs",3),rep("ws",3)),
                                      issues=c(rep("leaf",18),rep("root",17))
                                                ,pca1$x)

#rep("ws-1d-es",3),rep("ws-15d-sc",3),
sum1=summary(pca1)
library(ggplot2)
xlab1 <- paste0("PC1(",round(sum1$importance[2,1]*100,2),"%)")
ylab1 <- paste0("PC2(",round(sum1$importance[2,2]*100,2),"%)")
ylab3 <- paste0("PC3(",round(sum1$importance[2,3]*100,2),"%)")

da$ecotype=sapply(da$sample,function(x){strsplit(x,'_')[[1]][2]})

pdf("pca1_pca2_expression.pdf",width=6,height=4)
ggplot(data = da,aes(x = PC1,y = PC2,color=color,shape=factor(ecotype)))+
  scale_shape_manual(values = c(2,5))+
  geom_point(size = 3.5)+
  labs(x = xlab1,y = ylab1,color = "Condition",title = "PCA Scores Plot")+
  theme_bw()
dev.off()
pdf("pca1_pca3.pdf",width=6,height=4)
 ggplot(data = da,aes(x = PC1,y = PC3,color=color,shape=factor(ecotype)))+
   scale_shape_manual(values = c(2,5))+
  geom_point(size = 3.5)+
  labs(x = xlab1,y = ylab3,color = "Condition",title = "PCA Scores Plot")+
   theme_bw()
dev.off()

