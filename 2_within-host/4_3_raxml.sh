#!/bin/bash

if [ $# -ne 2 ];
then 
	echo "USAGE: bash 4_3_raxml.sh [input directory] [output directory]"
	exit 0
fi

args=($@)

# save the arguments as variables so they can be edited easily


for idx in $(seq 0 1); do
	
	# check whether the input and output directories end with a slash 
	# add a slash if they do not
	current="${args[idx]}"
	#echo $current	
	if [[ $current != */ ]];
	then
		args[idx]="${args[idx]}/"
		echo "Fixed: ${args[idx]}"	
	fi

	# verify that the provided arguments are valid directories 
	if [ ! -d ${args[idx]} ];then
		echo "Invalid directory provided. Exiting."
		echo "${args[idx]}"
		exit 0
	fi
 
done


for fullname in ${args[0]}*.fasta; do
	seed=$RANDOM
	echo $seed
	echo ${args[1]}
	treename="$(cut -d'.' -f1 <<< `basename $fullname`).tree"
	raxmlHPC -T 8 -f a -p $seed -x $seed -N 100 -s $fullname -n $treename -w ${args[1]} -m GTRGAMMA
done

