#!/bin/bash
msadir="/home/jpalmer/PycharmProjects/hiv-withinhost/4_1Accno/"
outdir="/home/jpalmer/PycharmProjects/hiv-withinhost/8Historian/pat"
for filename in /home/jpalmer/PycharmProjects/hiv-withinhost/7SampleTrees/rescaled_multi/*.tree.sample; do
    # this name contains the run ID
    name="$(cut -d'.' -f1 <<<"$(basename $filename)")"
    guide="$(cut -d'-' -f1 <<<"$name").fasta"
    outfile=$name"_recon.fasta"
    echo $filename
    echo $name
    echo $msadir$guide
    echo $outdir$outfile
    #~/historian/bin/historian r -vv -guide $msadir$guide -tree $filename -ancseq -output fasta > "$outdir$folder/$outfile"

done