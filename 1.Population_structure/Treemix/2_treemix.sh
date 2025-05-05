workdir="/usr_storage/lzq/work/dyy/treemix"
vcf_dir="/usr_storage/lzq/work/dyy/strcture_filt/5_GQ10DP5"
pop="$workdir/treemix_clust.txt"
step=$1
if [[ $step == "1" ]];then
less $workdir/freq/ZSP1912-L.LG01.frq.strat|grep 'CHR'>plas.frq.strat
for chr in ZSP1912-L.LG{01..19}
do
	plink --vcf $vcf_dir/${chr}_rmindel_Qfilter_rmbed_DP5GQ10.recode.vcf  --maf 0.05 --geno 0.03 --indep-pairwise 50 10 0.2 --allow-extra-chr --allow-no-sex  --make-bed --double-id --out $workdir/plink_maf0.05_ld/${chr} 
	plink --vcf $vcf_dir/${chr}_rmindel_Qfilter_rmbed_DP5GQ10.recode.vcf --extract $workdir/plink_maf0.05_ld/${chr}.prune.in --allow-extra-chr --allow-no-sex   --double-id  --freq --missing --within $pop --out $workdir/freq/$chr
	less $workdir/freq/$chr.frq.strat |grep -v 'CHR'>>plas.frq.strat
done
python plink2treemix.py plas.frq.strat > plas.treemix.input&
gzip plas.treemix.input
elif [[ $step == "2" ]];then
	for m in {0..5}
   	do
      		echo "treemix  -bootstrap  -i plas.treemix.input.gz  -root yiyang -o treemix_out/TreeMix.${m} -m ${m}  -k 500&"

      	done
fi

