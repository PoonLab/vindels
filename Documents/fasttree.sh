#!/bin/bash
for filename in /home/jpalme56/PycharmProjects/hiv-evolution-master/5_1_final/*.fasta; do
	fasttree -nt "$filename" > "${filename}.tree"
done
