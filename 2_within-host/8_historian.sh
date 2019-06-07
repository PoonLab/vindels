#!/bin/bash
msadir="/home/jpalmer/PycharmProjects/hiv-withinhost/4_1Accno/"
outdir="/home/jpalmer/PycharmProjects/hiv-withinhost/8Historian/"
for foldername in /home/jpalmer/PycharmProjects/hiv-withinhost/7SampleTrees/rescaled_multi/*; do
    folder=`basename $foldername`
    if [ ! -d "$outdir$folder/" ]; then
        mkdir "$outdir$folder/"
    fi

    for filename in $foldername/*.sample; do
        # this name contains the run ID
        name="$(cut -d'.' -f1 <<<"$(basename $filename)")"
        guide="$(cut -d'-' -f1 <<<"$name").fasta"
        outfile=$name"_recon.fasta"
        #echo $filename
        #echo $name
        #echo $guide
        ~/historian/bin/historian r -vv -guide $msadir$guide -tree $filename -ancseq -output fasta > "$outdir$folder/$outfile"
    done
done