#!/bin/bash
for filename in /home/jpalme56/PycharmProjects/hiv-evolution-master/5_3_AccHeaders/*.fasta; do
	fasttree -nt "$filename" > "${filename}.out"
done
