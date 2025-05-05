vcf="/usr_storage/lzq/work/dyy/strcture_filt/dyy201_rmindel_Qfilter_rmbed_DP5GQ10_mis0.8_maf0.0001.recode.vcf.gz"
pop="/usr_storage/lzq/work/dyy"
for chr in LG{01..19}
do
echo "python3 /usr_storage/lzq/software/conda_evns/xp-clr/bin/xpclr --format vcf --input $vcf --samplesA $pop/western.txt --samplesB $pop/eastern.txt --chr $chr --ld 0.7 --maxsnps 200 --size 1000 --step 1000 --out  /usr_storage/lzq/work/dyy/selection/xpclr/$chr.xpclr&"
done
cat LG*xpclr|grep -v 'chrom'|cut -f '2,3,4,12'|awk '{if ($4!="") print $0 }' > all.xpclr
