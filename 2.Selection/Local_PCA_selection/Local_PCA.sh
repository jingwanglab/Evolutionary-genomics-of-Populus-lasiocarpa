for chr in LG{01..19}
do
	cp data/sample_info.tsv $chr/data 
	cp sample_info.tsv data
	cp /usr_storage/lzq/work/dyy/local_adaption/env_filt/snp/$chr.recode.vcf   data
	bgzip data/$chr.recode.vcf && bcftools index data/$chr.recode.vcf.gz 
done
	Rscript run_lostruct.R -i data --nmds 10  -t snp -s 250 -I data/sample_info.tsv
	dir=`ls ./lostruct_results/|grep '250_weights'`
paste -d ',' $chr/lostruct_results/$dir/$chr.regions.csv <(cut -d "," -f 3,4,5 $chr/lostruct_results/$dir/mds_coords.csv ) >$chr.csv 
sed -i "s/\"//g" $chr.csv 

