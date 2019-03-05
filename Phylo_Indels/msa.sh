#!/bin/bash
for filename in /home/jpalme56/PycharmProjects/hiv-evolution-master/4_1_Edited/*.fasta; do
	mafft --auto "$filename" > "${filename}.msa"
done
