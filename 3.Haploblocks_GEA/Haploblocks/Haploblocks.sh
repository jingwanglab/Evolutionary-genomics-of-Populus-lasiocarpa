for chr in LG{01..19}
do
    mkdir $chr
    vcftools --gzvcf ../data/$chr.recode.vcf.gz --min-alleles 2 --maf 0.05 --recode --recode-INFO-all --out $chr/$chr.DP.GQ.bialleles.bed.mis0.8.maf0.05.180ind
    cp ../data/sample_info.tsv $chr/
    bgzip $chr/$chr.DP.GQ.bialleles.bed.mis0.8.maf0.05.180ind.recode.vcf && bcftools index $chr/$chr.DP.GQ.bialleles.bed.mis0.8.maf0.05.180ind.recode.vcf.gz
    Rscript run_lostruct.R -i $chr/ --nmds 5  -t snp -s 60 -I $chr/sample_info.tsv
done
