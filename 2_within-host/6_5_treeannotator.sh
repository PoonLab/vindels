# tree annotator bash script 
# used for batch processing the generation of MCC trees from BEAST outputs

if [ $# -ne 2 ]; then
    echo "Usage: bash 6_5_treeannotator.sh [input tree directory] [output name]"	    
    exit 0
fi

indir=$1
outdir=$2

if [[ $indir != */ ]];
then
        indir="$indir/"
        echo "Fixed: $indir"      
fi

if [[ $outdir != */ ]];
then
        outdir="$outdir/"
        echo "Fixed: $outdir"      
fi



output="/home/jpalmer/PycharmProjects/hiv-withinhost/7_5_MCC/prelim/$outdir"
if [ ! -d $output ];then
        echo "Making new directory."
        echo "$output"
        mkdir "$output"
fi

for filename in $indir*trees; do 
    outfile="$(cut -d'.' -f1 <<< `basename $filename`).tree"
    echo $filename
    echo $output$outfile
    treeannotator -heights median -burnin 10000000 $filename "$output$outfile"
done 
