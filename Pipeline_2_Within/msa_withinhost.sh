#!/bin/bash
dir=/home/jpalmer/PycharmProjects/hiv-withinhost/4MSA/
for filename in /home/jpalmer/PycharmProjects/hiv-withinhost/3RegionSequences/full_length/*.fasta; do
	echo $filename
	file="$(cut -d'/' -f8 <<<"$filename")"
	echo $dir$file
	mafft --auto "$filename" > "$dir$file"
done

