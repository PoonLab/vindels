#!/bin/bash
for filename in /home/jpalme56/PycharmProjects/hiv-evolution-master/PhylipFiles6/*.phy; do
	phyml -i "$filename" 
done
