#!/bin/bash
vcf="dyy201_rmindel_Qfilter_rmbed_DP5GQ10_mis0.8_maf0.0001.recode.vcf.gz"
shfile="/usr_storage2/syp/LZQ/SMCPP_boostrap/run_split.sh"
indir="/usr_storage2/syp/LZQ/SMCPP_boostrap/split_gz"
outdir="/usr_storage2/syp/LZQ/SMCPP_boostrap/split"
rundir="/usr_storage2/syp/LZQ/SMCPP_boostrap"
jasondir="/usr_storage2/syp/LZQ/SMCPP_boostrap/analysis"

if [ -f $shfile ];then
        /bin/rm -rf $shfile
fi

for boostrap in $indir/*
do
	strap=`basename $boostrap`
	if [ ! -d $outdir/east_west/$strap ];then
        	mkdir $outdir/east_west/$strap
        fi
        if [ ! -d $outdir/east_southwest/$strap ];then
                mkdir $outdir/east_southwest/$strap
        fi
        if [ ! -d $outdir/west_southwest/$strap ];then
                mkdir $outdir/west_southwest/$strap
        fi
	echo "cd $rundir;docker run --rm -v \$PWD:/mnt terhorst/smcpp:latest split -o split/east_west/$strap/east_west analysis/east/$strap/model.final.json analysis/west/$strap/model.final.json split_gz/$strap/east_west/*smc.gz" >> $shfile
        echo "cd $rundir;docker run --rm -v \$PWD:/mnt terhorst/smcpp:latest split -o split/east_southwest/$strap/east_southwest analysis/east/$strap/model.final.json analysis/southwest/$strap/model.final.json split_gz/$strap/east_southwest/*smc.gz" >> $shfile
        echo "cd $rundir;docker run --rm -v \$PWD:/mnt terhorst/smcpp:latest split -o split/west_southwest/$strap/west_southwest analysis/west/$strap/model.final.json analysis/southwest/$strap/model.final.json split_gz/$strap/west_southwest/*smc.gz" >> $shfile
done
