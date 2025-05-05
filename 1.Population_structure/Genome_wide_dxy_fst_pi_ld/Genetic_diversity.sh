#!/bin/bash
bedtools="/usr_storage/software/bedtools2/bin/bedtools"
pixy="/usr_storage/lzq/software/conda_evns/pixy/bin/pixy"
vcf="/usr_storage/lzq/work/mutation/0-4-pai"
workdir="/usr_storage/lzq/work/mutation/plas_0_4_pi"
species=$1
site=$2
outdir_pixy=$workdir/pixy_${site}_$species
if [ ! -d $outdir_pixy ]; then
mkdir -p $outdir_pixy
fi
site_vcf=$workdir/${site}_vcf_$species
if [ ! -d $site_vcf ]; then
mkdir -p $site_vcf
fi
###pai
for file in $vcf/ZSP1912-L.*_rmindel_Qfilter_bed.$species.recode.vcf
do
temp=${file##*/}
LG=${temp%%_*}
/usr_storage/software/bedtools2/bin/bedtools intersect -a $file -b plas_${site}_fold.bed -header > $site_vcf/$LG.new.recode.vcf 
/usr_storage/lzq/software/tabix-0.2.6/bgzip $site_vcf/$LG.new.recode.vcf && /usr_storage/lzq/software/tabix-0.2.6/tabix  $site_vcf/$LG.new.recode.vcf.gz 
echo "nohup $pixy --stats pi --vcf $site_vcf/$LG.new.recode.vcf.gz  --window_size 100000 --n_cores 5  --populations $species.pop  --output_folder pixy_${site}_$species  --output_prefix $LG.pai.out 2>&1&"
done
