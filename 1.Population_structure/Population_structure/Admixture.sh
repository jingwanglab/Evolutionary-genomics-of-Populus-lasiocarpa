#!/bin/bash
set -e
file="/export2/home/longzq/lzq/snp/struture/plink/lasiocarpa.bed"
outdir="/export2/home/longzq/lzq/snp/struture/admixture"
if [[ ! -d $outdir ]];then
	mkdir -p $outdir
fi
for K in {1..10} 
do
       echo " cd $outdir && admixture --cv $file $K -j20  |tee log${K}.out  " >command6/ad$K.command
done
