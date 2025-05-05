cds_sift="/usr_storage/lzq/work/mutation/sift/anno/result/dayeyang/dayeyang.201.rmindel_snp.biallelic.DP5.GQ10.bed.dayeyang.mis0.8_SIFTannotations.xls"
sift_dir="/usr_storage/lzq/work/mutation/derived/aim_snp"
grep "START\|STOP" $cds_sift > $sift_dir/plas201.cds.loss_of_function.xls
awk '$9=="NONSYNONYMOUS"' $cds_sift |awk '$17=="DELETERIOUS"' > $sift_dir/plas201.cds.deleterious.xls
awk '$9=="NONSYNONYMOUS"' $cds_sift |grep -v "DELETER" |awk '$13!="NA"'> $sift_dir/plas201.cds.tolerated.xls
awk '$9=="SYNONYMOUS"' $cds_sift |grep -v "DELETERIOUS" |awk '$13!="NA"' > $sift_dir/plas201.cds.synonymous.xls
awk '$8=="CDS"' $cds_sift |awk '$13!="NA"' >$sift_dir/all.cds
##extracting snp id
cut -f 1,2 $sift_dir/plas201.cds.loss_of_function.xls |sed 's/\t/:/g'> $sift_dir/plas201.cds.loss_of_function.snps
cut -f 1,2 $sift_dir/plas201.cds.deleterious.xls |sed 's/\t/:/g'> $sift_dir/plas201.cds.deleterious.snps
cut -f 1,2 $sift_dir/plas201.cds.tolerated.xls |sed 's/\t/:/g'> $sift_dir/plas201.cds.tolerated.snps
cut -f 1,2 $sift_dir/plas201.cds.synonymous.xls |sed 's/\t/:/g'> $sift_dir/plas201.cds.synonymous.snps
cut -f 1,2 $sift_dir/all.cds|sed 's/\t/:/g'> $sift_dir/plas201.cds.all.snps
