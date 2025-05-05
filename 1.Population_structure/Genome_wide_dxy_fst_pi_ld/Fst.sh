#!/bin/bash
pixy="/usr_storage/lzq/software/conda_evns/pixy/bin/pixy"
vcf=" /usr_storage/lzq/work/dyy/strcture_filt/6_dyy"
workdir="/usr_storage/lzq/work/dyy/dxy"
pair=$1
outdir_pixy=fst_out_${pair}
if [ ! -d $outdir_pixy ]; then
mkdir -p $outdir_pixy
fi

for file in $vcf/*rmindel_Qfilter_rmbed_DP5GQ10_dyy_mis0.8.recode.vcf.gz
do
temp=${file##*/}
LG=${temp%%_*}
echo "$pixy --stats fst --vcf $file  --window_size 100000 --n_cores 5  --populations $pair.txt  --output_folder $outdir_pixy  --bypass_invariant_check yes   --output_prefix $LG.pai.out &"
done

