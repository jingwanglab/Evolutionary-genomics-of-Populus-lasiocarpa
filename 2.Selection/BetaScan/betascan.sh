#!/bin/bash
pop="/usr_storage/lzq/work/dyy"
dir="/usr_storage/lzq/work/dyy/selection/clr"
betascan="/usr_storage/lzq/software/BetaScan-master/BetaScan.py"
python vcf2betascan.py /usr_storage/lzq/work/dyy/selection/balance_selection/derived_vcf/ /usr_storage/lzq/work/dyy/selection/balance_selection/beta_input/ 
for juqun in {western,eastern}
do
	for i in LG{01..19}
	do
		vcftools --gzvcf $dir/dyy_outgroup_instersect.vcf.gz --snps  $dir/allele_freq/${juqun}.$i.derived.snps --keep $juqun.15.txt  --recode --recode-INFO-all --out derived_vcf/${juqun}.$i.derived 
		/usr_storage/wjl/software/conda_evns/hicpro/bin/python2 $betascan -i /usr_storage/lzq/work/dyy/selection/balance_selection/beta_input/${juqun}.$i.derived.recode.vcf -w 1000  -o beta_out/${juqun}.$i.beta
	done
done


#for i in LG{01..19};do cat beta_out/eastern.$i.beta|grep -v 'Position'|awk '{print "'$i'""\t"$0}' >> eastern.beta ;done
#for i in LG{01..19};do cat beta_out/western.$i.beta|grep -v 'Position'|awk '{print "'$i'""\t"$0}' >> western.beta;done
