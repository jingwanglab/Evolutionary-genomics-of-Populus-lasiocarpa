# -*- coding: utf-8 -*-
"""
Created on Mon Apr  4 15:23:48 2022

@author: Administrator
"""
import gffutils

path="./new/"
#db = gffutils.create_db(path+'Plas.genome.gff3', dbfn = path+'las.db', force = True, keep_order = True, merge_strategy = 'merge', sort_attribute_values = True)
db = gffutils.FeatureDB('./local_PCA/las.db', keep_order = True)
chr="LG10"
17039910
S=17029910
E=17049910
gffGene = ['molecule,gene,start,end,strand,direction']
gffSubGene = ['molecule,gene,start,end,direction,subgene,from,to']
for gene in list(db.region(region=(chr,S,E),featuretype='gene', completely_within=True)):
    if gene.strand == '+':
        strand = 'forward'
        direction = '1'
    else:
        strand = 'reverse'
        direction = '-1'
    gffGene.append(gene.chrom+ ','+ gene.id+ ',' + str(gene.start) + ',' + str(gene.end) + ',' + strand + ',' + direction)       
    length=[]
    if len(list(db.children(gene.id,featuretype="mRNA")))>1:
        for temp in list(db.children(gene.id,featuretype="mRNA")):
            length.append(temp.bin)
        gene.id=list(db.children(gene.id,featuretype="mRNA"))[length.index(max(length))].id
    for subgene in db.children(gene.id,featuretype=['exon','three_prime_UTR', 'five_prime_UTR']):
        if '.t' in gene.id:
            ids=gene.id.split(".t")
            gene.id=ids[0]
        gffSubGene.append(gene.chrom+ ','+ gene.id + ',' + str(gene.start) + ',' + str(gene.end) + ',' + direction + ',' + subgene.featuretype + ',' + str(subgene.start) + ',' +str(subgene.end))


def write2file(out_path,list_in):
     f1=open(out_path,'w')
     out= '\n'.join(list_in)
     out += '\n'
     f1.write(out+"\n")
     f1.close() 

write2file(path+'distribution/'+chr+"_"+str(S)+'_gene.csv',gffGene)
write2file(path+'distribution/'+chr+"_"+str(S)+'_subgene.csv',gffSubGene)

