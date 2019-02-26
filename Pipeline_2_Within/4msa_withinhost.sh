#!/bin/bash
dir=/home/jpalmer/PycharmProjects/hiv-withinhost/4MSA/
for filename in /home/jpalmer/PycharmProjects/hiv-withinhost/3_RegionSequences/full_length/*.fasta; do
	echo $filename
	file="$(cut -d'/' -f8 <<<"$filename")_"
	echo $dir$file
	mafft --auto "$filename" > "$dir$file"
done

