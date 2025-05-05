#!/bin/bash
# conda activate /usr_storage/cly/.conda/envs/py35/
vcf="/usr_storage/lzq/work/dyy/local_adaption/GF/snp_indel_sv.maf5.mis8.vcf"
plink --vcf $vcf --allow-extra-chr --allow-no-sex  --make-bed  --out gemma_input
paste -d " " <(cut -d " " -f 1,2,3,4,5 gemma_input.fam) <(sed '1d' 180.clim.points |cut -d " " -f 4- ) > gemma_input1.fam
paste -d " " <(cut -d " " -f 1,2,3,4,5 gemma_input.fam) <(less 180.clim.points|grep 
mv gemma_input1.fam gemma_input.fam
gemma -bfile gemma_input -gk 2 -o gemma_matrix
for i in {1..19}
do	
	num=$(($i+3))
	bio=`head -n 1 180.clim.points|cut -d " " -f $num`
	gemma -bfile gemma_input -n $i -k output/gemma_matrix.sXX.txt -lmm -o  $bio.gemma
done
