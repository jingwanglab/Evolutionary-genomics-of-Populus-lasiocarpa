##ls ../*txt|while read id;do name=${id##*/};echo " awk '{print \$1}' $id > $name";done
PopLDdecay="/usr_storage/lsy/software/PopLDdecay/bin/PopLDdecay"
vcf="/usr_storage/lzq/work/dyy/strcture_filt/dyy201_rmindel_Qfilter_rmbed_DP5GQ10_mis0.8_maf0.0001.recode.vcf"
workdir="/usr_storage/lzq/work/dyy/LD"
pop=$1
$PopLDdecay -InVCF $vcf -MaxDist 200 -MAF 0.05 -SubPop $pop.txt  -OutStat $pop.200kb.Lddecay.stat 
perl /usr_storage/lsy/software/PopLDdecay/bin/Plot_MultiPop_backup.pl -inList all.list -output all_ld -keepR -bin1 100 -bin2 1000 -break 5000 
perl /usr_storage/lsy/software/PopLDdecay/bin/Plot_OnePop.pl -inFile $pop.200kb.Lddecay.stat.gz -output $pop.200kb.Lddecay 

