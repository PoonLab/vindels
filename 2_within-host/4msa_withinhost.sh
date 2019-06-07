#!/bin/bash
dir=/home/jpalmer/PycharmProjects/hiv-withinhost/4_1Accno/
for filename in /home/jpalmer/PycharmProjects/hiv-withinhost/3_RegionSequences/full_length/*.fasta; do
	echo $filename
	file="$(cut -d'/' -f8 <<<"$filename")2"
	echo $dir$file
	mafft --auto "$filename" > "$dir$file"
done

