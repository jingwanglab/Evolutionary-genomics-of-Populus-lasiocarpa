# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 21:33:15 2021

@author: Administrator
"""
import os
import csv
import sys
if len(sys.argv)!=3:
    sys.stderr.write('Usage:python vcf2sweepfinder2.py vcf_dir outdir')
    sys.exit()

files=os.listdir(sys.argv[1])
for file in files:
    if file[-3:]=='log' or "western" in file:
        continue
    file_path=sys.argv[1]+file
    vcf=open(file_path,'r')
    sweep_out = open(sys.argv[2]+file, 'w')
    sweep_out.write("position\tx\tn\tfolded\n") # write header
    folded = "0"
    sample_names = []
    sample_seqs = []
    vcf_data=csv.reader(vcf,delimiter='\t')
    ###################读取vcf的每一行
    for line in vcf_data:
        if any('##' in strings for strings in line ):
            continue
        if any('#CHROM' in strings for strings in line):
            sample_names = line[9:]
            sample_seqs=[['']]*len(sample_names)    
            continue
        chrom,pos,id,ref,alt,qual,filter,info,format=line[0:9]
        haplotypes = line[9:]
        location=pos
        x = 0
        n = 0
        for index,haplotype in enumerate(haplotypes):
            if haplotype.split(":")[0] != './.':
                n = n + (2-haplotype.split(":")[0].count(".")) # count up the non missing calls ,modified by L-of-IOS
                if haplotype.split(":")[0] != "0/0" and haplotype.split(":")[0] != ".":
                    x = x + (2-haplotype.split(":")[0].count("0")-haplotype.split(":")[0].count("."))  # count derived alleles = not reference and not missing modified by  L-of-IOS
        if x != 0 and x != n:
            sweep_out.write(location+"\t"+str(x)+"\t"+str(n)+"\t"+folded+"\n")
    vcf.close()
    sweep_out.close()
