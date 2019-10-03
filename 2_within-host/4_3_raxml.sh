#!/bin/bash

if [ $# -ne 1 ];
then 
	echo "USAGE: bash 4_3_raxml.sh [input directory]"
	exit 0
fi

# save the arguments as variables so they can be edited easily
inDir=$1

# check whether the input and output directories end with a slash 
# add a slash if they do not 
if [ ${inDir:(-1)} != "/" ];
then
	inDir="$1/"	
fi


# verify that the provided arguments are valid directories 
if [ ! -d $inDir ]; then
	echo "Invalid first directory provided. Exiting."
	exit 0
fi

for fullname in $inDir*.fasta; do
	echo $fullname
	treename="$(cut -d'.' -f1 <<< `basename $fullname`).tree"
	raxmlHPC -s $fullname -p 123 -n $treename -m GTRGAMMA
done

