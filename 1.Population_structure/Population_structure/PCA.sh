#!/bin/bash
#1. performing pca analysis
smartpca -p pk.txt > pca_out.log
#！！！！！！！！！！！！
#2. using the twstats to select the most significant PCs
twstats -t twtable -i shanyang.pca.eval -o shanyang.pca.out

