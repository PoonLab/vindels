msadir="/home/jpalmer/PycharmProjects/hiv-withinhost/4MSA/hm-screen/"
treedir="/home/jpalmer/PycharmProjects/hiv-withinhost/4_5_Raxml/100BS/"
echo "${treedir}rooted_trees/"

for fullname in ${treedir}rooted_trees/*.tree; do
	filename="$(cut -d'.' -f1 <<< `basename $fullname`)"
	echo "Start $filename"
	output="${treedir}treetime_test2/${filename}/"
	mkdir $output
	treetime --tree $fullname --aln "$msadir$filename.fasta" --dates "${treedir}treetime_test2/dates/${filename}.tsv" --gtr HKY --outdir $output
	
done

