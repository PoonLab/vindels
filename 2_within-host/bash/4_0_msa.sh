#!/bin/bash
if [ $# -ne 2 ];
then 
	echo "USAGE: bash 4_0_msa.sh [input directory] [output directory]"
	exit 0
fi

# save the arguments as variables so they can be edited easily
inDir=$1
outDir=$2

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


# verify that the provided arguments are valid directories 
if [ ! -d $inDir ]; then
	echo "Invalid first directory provided. Exiting."
	exit 0
elif [ ! -d $outDir ]; then 
	echo "Invalid second directory provided. Exiting."
	exit 0
fi

for filename in $inDir*.fasta; do
	file=`basename $filename`
	#echo $file
	echo $outDir$file
	mafft --auto --reorder "$filename" > "$outDir$file"
done
q

