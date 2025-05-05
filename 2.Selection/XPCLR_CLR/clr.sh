sweepfinder2="/usr_storage/lzq/software/SF2/SweepFinder2"
python vcf2sweepfinder2.py  /usr_storage/lzq/work/dyy/selection/clr/allele_freq/derived_vcf/ /usr_storage/lzq/work/dyy/selection/clr/allele_freq/sweepfinder_out/ &
head -n 1  allele_freq/sweepfinder_out/western.LG19.derived.recode.vcf >header
cat  allele_freq/sweepfinder_out/western* |grep -v 'position'> allele_freq/sweepfinder_out/all_western.temp
cat  allele_freq/sweepfinder_out/eastern*|grep -v 'position'> allele_freq/sweepfinder_out/all_eastern.temp
cat header allele_freq/sweepfinder_out/all_western.temp>allele_freq/sweepfinder_out/all_western 
cat header allele_freq/sweepfinder_out/all_eastern.temp >allele_freq/sweepfinder_out/all_eastern
rm allele_freq/sweepfinder_out/all_western.temp allele_freq/sweepfinder_out/all_eastern.temp
echo "$sweepfinder2 -f allele_freq/sweepfinder_out/all_western  allele_freq/sweepfinder_out/all_western.freq &"
echo "$sweepfinder2 -f allele_freq/sweepfinder_out/all_eastern  allele_freq/sweepfinder_out/all_eastern.freq &"
	for file in allele_freq/sweepfinder_out/*vcf
	do
		temp=${file##*/}
		pop=${temp%%.*}
		echo " cd $workdir/allele_freq/sweepfinder_out/ && $sweepfinder2 -lg 2000 $temp all_${pop}.freq $workdir/sweep_out/$temp.sweep &"
	done

fi

