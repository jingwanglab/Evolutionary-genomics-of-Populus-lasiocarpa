# -*- coding: utf-8 -*-
"""
Created on Thu May  5 16:57:30 2022

@author: Administrator
"""
import  os
import pandas as pd
import csv
import re
import random
gemma_path='./new/gk1_all/'
lfmm_path='./new/lfmm_5e2/'     
bios=os.listdir(lfmm_path)
out=[]
loc_bio=[]
for bio in bios:
    gemma_bio=gemma_path+bio
    gemma=pd.read_csv(gemma_bio,sep=" ")['rs'].unique().tolist()
    lfmm_bio=lfmm_path+bio
    lfmm=pd.read_csv(lfmm_bio,sep=" ")['id'].unique().tolist()
    intersect=set(gemma)&set(lfmm)
    if len(intersect) > 0 :
        out_path="./new/asso_loci/"+bio
        loc_bio.append(bio)
        out.append(intersect)
dict_out=dict(zip(loc_bio,out))
import itertools 
candidate_loc=set(list(itertools.chain.from_iterable(out)))

ids=pd.read_csv('./new/id.txt',sep='\t',header=None)[0].tolist()
back=set(ids)-set(candidate_loc)
back1=random.sample(back, 50000) 
#chr_pos['chr_pos']=chr_pos[0]+'_'+chr_pos[1].map(str)
frq_dir="./new/culmu_importance/freq/"
files=os.listdir(frq_dir)
pop=[]
frq_dict={}
dict_out['ref']=back1
for bio in dict_out.keys():
    frq_dict[bio]=[]
   


for file in files:    
    if '.frq' in file: 
        print(file)
        list1=file.split(".frq")
        pop.append(list1[0])
        path=frq_dir+file
        frq=pd.read_csv(path,sep="\t",usecols=[0,1,5],skiprows=[0],names=['chr','pos','alt'])
        frq["alt"]=frq["alt"].str.split(":",expand=True)[1]
        frq.index=ids
        for bio in dict_out.keys():
            frq_dict[bio].append(frq.loc[dict_out[bio],'alt'].tolist())
                    
for bio in frq_dict.keys():
        a=pd.DataFrame.from_dict(zip(pop,frq_dict[bio]))
        a.columns=['pop','alt']
        a=pd.concat([a.drop('alt', axis=1), pd.DataFrame(a['alt'].tolist())], axis=1)
        a.set_index(['pop'],inplace=True)
        a.columns=dict_out[bio]
        frq_out='./new/culmu_importance/'+"bio_frq/"+bio
        a.to_csv(frq_out,header=True,sep='\t',index=True)
#########
