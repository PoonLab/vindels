#!/bin/bash
msadir="/home/jpalmer/PycharmProjects/hiv-withinhost/4_1Accno/"
outdir="/home/jpalmer/PycharmProjects/hiv-withinhost/8Historian/"
for filename in /home/jpalmer/PycharmProjects/hiv-withinhost/7SampleTrees/rescaled/*.sample; do
    name="$(cut -d'_' -f1 <<<"$(basename $filename)")"
    guide="$(cut -d'-' -f1 <<<"$name").fasta"
    outfile=$name"_recon.fasta"
    ~/historian/bin/historian r -vvv -guide $msadir$guide -tree $filename -ancseq -output fasta > $outdir$outfile
done
