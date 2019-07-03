# tree annotator bash script 
# used for batch processing the generation of MCC trees from BEAST outputs

folder="/home/jpalmer/PycharmProjects/hiv-withinhost/6BEASTout3/trees/"
outdir="/home/jpalmer/PycharmProjects/hiv-withinhost/7_5_MCC/"

for filename in $folder*trees; do 
    outfile="$(cut -d'.' -f1 <<< `basename $filename`).tree"
    echo $filename
    echo $outdir$outfile
    treeannotator -heights median -burnin 10000000 $filename $outdir$outfile
done 