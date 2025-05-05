#! /bin/bash -l
beagle="/usr_storage/lzq/software/beagle.27Jan18.7e1.jar"
##-----Input vcf file 
Inputvcf="/usr_storage/lzq/work/dyy/strcture_filt/dyy201_rmindel_Qfilter_rmbed_DP5GQ10_mis0.8_maf0.0001.recode.vcf"

##------BEAGLE imputation

##the first is to impute the missing genotype
java -Xmx100g -jar $beagle gtgl=$Inputvcf out=dyy201_rmindel_Qfilter_rmbed_DP5GQ10_mis0.8_maf0.0001_beagle nthreads=32 window=100000 overlap=10000
java -Xmx100g -jar $beagle gt=dyy201_rmindel_Qfilter_rmbed_DP5GQ10_mis0.8_maf0.0001_beagle.vcf.gz out=dyy201_rmindel_Qfilter_rmbed_DP5GQ10_mis0.8_maf0.0001_beagle_phase nthreads=20 window=100000 overlap=10000 ibd=true ibdtrim=100 ibdlod=5
###estimate IBD shared haplotypes
for chr in LG{01..19}
do
java -Xmx40g -jar $beagle gt=dyy201_rmindel_Qfilter_rmbed_DP5GQ10_mis0.8_maf0.0001_beagle_phase.vcf.gz out=dyy201_rmindel_Qfilter_rmbed_DP5GQ10_mis0.8_maf0.0001_beagle_phase.ibd.$chr   chrom=$chr ibd=true impute=false window=100000 overlap=10000 ibdtrim=100 ibdlod=5 ibdcm=1e-5 
done
