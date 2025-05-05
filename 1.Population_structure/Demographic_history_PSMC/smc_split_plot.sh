#!/bin/bash
shfile="/usr_storage2/syp/LZQ/SMCPP_boostrap/plot_split.sh"
rundir="/usr_storage2/syp/LZQ/SMCPP_boostrap/split"

if [ -f $shfile ];then
        /bin/rm -rf $shfile
fi

for pairse in $rundir/*
do
	species=`basename $pairse`
	for strap in boostrap_{1..10}
	do
		echo "cd $pairse/$strap/$species;docker run --rm -v \$PWD:/mnt terhorst/smcpp:latest plot -g 15 -x 1e3 5.5e5 $species.$strap.pdf model.final.json -c" >> $shfile
	done
done
