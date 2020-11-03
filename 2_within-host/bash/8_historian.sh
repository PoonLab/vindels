#!/bin/bash
path="/home/jpalmer/PycharmProjects/hiv-withinhost/"
msadir="${path}4MSA/hm-screen/"
outdir="${path}8Historian/mcc/local/"
indir="${path}7_5_MCC/final/"
for filename in $indir*.tree; do
    name="$(cut -d'.' -f1 <<< `basename $filename` )"
    guide="$(cut -d'-' -f1 <<< $name).fasta"
    outfile="${name}_recon.fasta"
    #echo $filename
    #echo $name
    #echo $msadir$guide
    #echo $outdir$outfile
    ~/historian/bin/historian -v -guide $msadir$guide -tree $filename -ancseq -output fasta > "$outdir$outfile"

done