
sift4g_annotator="/usr_storage/lzq/software/SIFT4G_Annotator-master/SIFT4G_Annotator.jar"
out="/usr_storage/lzq/work/mutation/sift/anno/result/dayeyang"
if [ ! -d "$out" ]; then
mkdir -p $out
fi
vcf="/usr_storage/lzq/work/mutation/sift/anno/dayeyang.201.rmindel_snp.biallelic.DP5.GQ10.bed.dayeyang.mis0.8.vcf"
echo "java -jar $sift4g_annotator -c -i $vcf -d /usr_storage/lzq/work/mutation/sift/anno/ -r $out"



