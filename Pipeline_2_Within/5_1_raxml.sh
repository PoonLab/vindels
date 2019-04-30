#!/bin/bash

# save the arguments as variables so they can be edited easily
inDir=$1
outDir=$2

# verify that the provided arguments are valid directories 
if [ ! -d $inDir ]; then
	echo "Invalid first directory provided. Exiting."
	exit 0
elif [ ! -d $outDir ]; then 
	echo "Invalid second directory provided. Exiting."
	exit 0
fi

# check whether the input and output directories end with a slash 
# add a slash if they do not 
if [ ${inDir:(-1)} != "/" ];
then
	inDir="$1/"	
fi

if [ ${outDir:(-1)} != "/" ];
then
	outDir="$2/"	
fi

query="RAxML_bestTree"

for fullname in $inDir$query*; do
	echo $fullname
	#name=`basename $fullname`
	echo $outDir$name
	#raxmlHPC -s $fullname -p 123 -n $outDir$name -m GTRGAMMA
done

