step=$1
###extract the cds snp in outgroup 
if [[ $step == "1" ]];then
outgroup="/usr_storage/lzq/work/mutation/dfe_new/outgroup"
cds="/usr_storage/lzq/work/mutation/derived/aim_snp/plas201.cds.all.snps"
out1="ZHY-03-4"
out2="trichocarpa1"
#cp $outgroup/$out1.nomiss.noindel.recode.vcf ./
#cp $outgroup/$out2.nomiss.noindel.recode.vcf ./
#bcftools view -O z -o $out1.nomiss.noindel.recode.vcf.gz $out1.nomiss.noindel.recode.vcf 
#bcftools view -O z -o $out2.nomiss.noindel.recode.vcf.gz  $out2.nomiss.noindel.recode.vcf
#bcftools index $out1.nomiss.noindel.recode.vcf.gz 
#bcftools index $out2.nomiss.noindel.recode.vcf.gz
#bcftools annotate --set-id +'%CHROM\:%POS' -O z -o $out1.nomiss.noindel.annote.recode.vcf.gz $out1.nomiss.noindel.recode.vcf.gz
#bcftools annotate --set-id +'%CHROM\:%POS' -O z -o $out2.nomiss.noindel.annote.recode.vcf.gz $out2.nomiss.noindel.recode.vcf.gz
#vcftools --gzvcf $out1.nomiss.noindel.annote.recode.vcf.gz --snps $cds --recode --recode-INFO-all --out $out1.nomiss.noindel.cds
bgzip $out1.nomiss.noindel.cds.recode.vcf && tabix $out1.nomiss.noindel.cds.recode.vcf.gz
#vcftools --gzvcf $out2.nomiss.noindel.annote.recode.vcf.gz --snps $cds --recode --recode-INFO-all --out $out2.nomiss.noindel.cds
bgzip $out2.nomiss.noindel.cds.recode.vcf && tabix $out2.nomiss.noindel.cds.recode.vcf.gz
elif [[ $step == "2" ]];then
vcf="/usr_storage/lzq/work/mutation/sift/anno/dayeyang.201.rmindel_snp.biallelic.DP5.GQ10.bed.dayeyang.mis0.8.vcf"
snp_type=$2
species=$3
out1="ZHY-03-4"
out2="trichocarpa1"
sift_dir="/usr_storage/lzq/work/mutation/derived/aim_snp"
out_vcf="/usr_storage/lzq/work/mutation/derived/pop_vcf"
pop="/usr_storage/lzq/work/mutation/dfe_new"
#vcftools --vcf $vcf --snps $sift_dir/plas201.cds.${snp_type}.snps --keep $pop/$species.txt --recode --recode-INFO-all --out $out_vcf/${species}.snp.${snp_type}
#bgzip $out_vcf/${species}.snp.${snp_type}.recode.vcf && tabix $out_vcf/${species}.snp.${snp_type}.recode.vcf.gz
#bcftools annotate -x FORMAT -O z -o $out_vcf/${species}.snp.${snp_type}.simple.vcf.gz $out_vcf/${species}.snp.${snp_type}.recode.vcf.gz 
#bcftools index $out_vcf/${species}.snp.${snp_type}.simple.vcf.gz
#bcftools isec -c all -n=3 $out_vcf/${species}.snp.${snp_type}.simple.vcf.gz $out1.nomiss.noindel.cds.recode.vcf.gz $out2.nomiss.noindel.cds.recode.vcf.gz  -w 1 |grep -v "#" | cut -f 3 > $out_vcf/${species}.snp.${snp_type}.intersect.snp.txt 
#vcftools --gzvcf $out_vcf/${species}.snp.${snp_type}.simple.vcf.gz --snps $out_vcf/${species}.snp.${snp_type}.intersect.snp.txt --recode --recode-INFO-all --out $out_vcf/${species}.snp.${snp_type}.simple.intersect
#vcftools --gzvcf $out1.nomiss.noindel.cds.recode.vcf.gz --snps $out_vcf/${species}.snp.${snp_type}.intersect.snp.txt --recode --recode-INFO-all --out $out1.$species.${snp_type}.nomiss.noindel.cds.intersect
#vcftools --gzvcf $out2.nomiss.noindel.cds.recode.vcf.gz --snps $out_vcf/${species}.snp.${snp_type}.intersect.snp.txt --recode --recode-INFO-all --out $out2.$species.${snp_type}.nomiss.noindel.cds.intersect
#bgzip $out_vcf/${species}.snp.${snp_type}.simple.intersect.recode.vcf  &&tabix $out_vcf/${species}.snp.${snp_type}.simple.intersect.recode.vcf.gz
#bgzip $out1.$species.${snp_type}.nomiss.noindel.cds.intersect.recode.vcf && tabix $out1.$species.${snp_type}.nomiss.noindel.cds.intersect.recode.vcf.gz
#bgzip $out2.$species.${snp_type}.nomiss.noindel.cds.intersect.recode.vcf  &&tabix $out2.$species.${snp_type}.nomiss.noindel.cds.intersect.recode.vcf.gz
bcftools merge -m snps $out_vcf/${species}.snp.${snp_type}.simple.intersect.recode.vcf.gz $out1.$species.${snp_type}.nomiss.noindel.cds.intersect.recode.vcf.gz $out2.$species.${snp_type}.nomiss.noindel.cds.intersect.recode.vcf.gz -O z -o ${species}.${snp_type}.outgroup.vcf.gz
elif [[ $step == "3" ]];then
pop="/usr_storage/lzq/work/mutation/dfe_new"
snp_type=$2
species=$3
out1="ZHY-03-4"
out2="trichocarpa1"
#vcftools --gzvcf ${species}.${snp_type}.outgroup.vcf.gz --keep $pop/$species.txt --freq --out est_sfs/${species}.${snp_type}.est_sfs
#vcftools --gzvcf ${species}.${snp_type}.outgroup.vcf.gz --indv $out1 --freq --out est_sfs/${species}.${snp_type}.utgroup1.est_sfs
#vcftools --gzvcf ${species}.${snp_type}.outgroup.vcf.gz --indv $out2 --freq --out est_sfs/${species}.${snp_type}.utgroup2.est_sfs
num=`wc -l $pop/$species.txt`
perl est_sfs/est_sfs_input.dfe.pl est_sfs/${species}.${snp_type}.est_sfs.frq est_sfs/${species}.${snp_type}.utgroup1.est_sfs.frq est_sfs/${species}.${snp_type}.utgroup2.est_sfs.frq $((${num:0:2}*2))
est_sfs_tool="/usr_storage/lzq/software/est-sfs-release-2.03/est-sfs"
config_rate6="/usr_storage/lzq/software/est-sfs-release-2.03/config-rate6.txt"
seed="/usr_storage/lzq/software/est-sfs-release-2.03/seedfile.txt"
config_jc="/usr_storage/lzq/software/est-sfs-release-2.03/config-JC.txt"
config_kimura="/usr_storage/lzq/software/est-sfs-release-2.03/config-kimura.txt"
$est_sfs_tool $config_rate6 est_sfs/${species}.${snp_type}.est_sfs.input.txt $seed est_sfs/${species}.${snp_type}.est_sfs.output.txt est_sfs/${species}.${snp_type}.est_sfs.output.p_anc.txt
$est_sfs_tool $config_jc est_sfs/${species}.${snp_type}.est_sfs.input.txt $seed est_sfs/${species}.${snp_type}.est_sfs.jc.output.txt est_sfs/${species}.${snp_type}.est_sfs.jc.output.p_anc.txt
$est_sfs_tool $config_kimura est_sfs/${species}.${snp_type}.est_sfs.input.txt $seed est_sfs/${species}.${snp_type}.est_sfs.kimura.output.txt est_sfs/${species}.${snp_type}.est_sfs.kimura.output.p_anc.txt
elif [[ $step == "4" ]];then
snp_type=$2
species=$3
outdir="/usr_storage/lzq/work/mutation/derived"
perl allele_freq/est_sfs_ancestral_derived.pl $outdir/est_sfs/${species}.${snp_type}.est_sfs.frq $outdir/est_sfs/${species}.${snp_type}.est_sfs.output.p_anc.txt
fi