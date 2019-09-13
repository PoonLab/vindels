# tree annotator bash script 
# used for batch processing the generation of MCC trees from BEAST outputs

if [ $# -ne 1 ]; then
    echo "Usage: bash 6_5_treeannotator.sh [input tree directory]"	    
    exit 0
fi
folder=$1
outdir="/home/jpalmer/PycharmProjects/hiv-withinhost/7_5_MCC/prelim/"

for filename in $folder*trees; do 
    outfile="$(cut -d'.' -f1 <<< `basename $filename`).tree"
    echo $filename
    echo $outdir$outfile
    treeannotator -heights median -burnin 20000000 $filename $outdir$outfile
done 
