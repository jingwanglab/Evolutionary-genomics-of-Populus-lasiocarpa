library(DESeq2)
library(data.table)
#library(VennDiagram)
library("ggplotify")
library(dplyr)
library(tibble)
library(gridExtra)
library(ggtern)
setwd("/Users/zhiqinlong/Desktop/computer/work/lasiocarpa/expressionnew/")
all_count=as.data.frame(fread('es_sc_all_gene_count'))
row.names(all_count)<-all_count$gene_id
all_count=all_count[,-1]
es_sc=as.data.frame(fread("es_sc_all_gene_tpm_filt_lowexp.csv",sep=","))

all_count=all_count[es_sc$gene_id,]
colnames(all_count)
all_count=all_count[,grep('hs|ws|ck',colnames(all_count))]

colDate<-data.frame(sampleID=colnames(all_count),
                    condition = c(rep('leaf',18),rep('root',17)),
                    population=c(rep(c('es','sc'),each=9),rep('es',8),rep('sc',9)),
                    treat=c(rep(c('ck','hs','ws'),each=3),
                            rep(c('ck','hs','ws'),each=3),
                            rep('ck',2),rep(c('hs','ws'),each=3),
                            rep(c('ck','hs','ws'),each=3)))
tissues=c("leaf","root")
treatments=c("hs","ws")
##################
for (tissue in tissues){
  for (treatment in treatments){
    colDate1=subset(colDate,condition==tissue&(treat==treatment|treat=='ck'))
    all_count1=all_count[,colDate1$sampleID]
    dds <- DESeqDataSetFromMatrix(countData = all_count1, colData =colDate1, design = ~ population + treat + population:treat)
    mod.dds=model.matrix(design(dds),colData(dds))
    es_treatment=colMeans(mod.dds[dds$population=="es"&dds$treat==treatment,])
    es_ck=colMeans(mod.dds[dds$population=="es"&dds$treat=="ck",])
    sc_treatment=colMeans(mod.dds[dds$population=="sc"&dds$treat==treatment,])
    sc_ck=colMeans(mod.dds[dds$population=="sc"&dds$treat=="ck",])
    treatment_dseq=colMeans(mod.dds[dds$treat==treatment,])
    dds2 <- DESeq(dds, test = "Wald")
    save(file=paste(tissue,"_",treatment,"dds2",sep=""),dds2)
    load(paste("./",tissue,"_",treatment,"dds2",sep=""))
    res_ck_es_sc <- results(dds2,pAdjustMethod = "BH",contrast = es_ck - sc_ck)
    res_hs_es_sc <- results(dds2,pAdjustMethod ="BH",contrast = es_hs - sc_hs)
    res_ck_es_sc<-res_ck_es_sc[order(res_ck_es_sc$padj),]
    res_ck_es_sc<-subset(res_ck_es_sc,res_ck_es_sc$baseMean>=10)
    res_hs_es_sc<-res_ck_es_sc[order(res_hs_es_sc$padj),]
    res_hs_es_sc<-subset(res_hs_es_sc,res_hs_es_sc$baseMean>=10)
     }}

