#!/bin/bash
inDir=$1

if [ ${inDir:(-1)} != "/"];
then
    inDir="$1/"
fi

for filename in $inDir*.xml; do
    beast $filename
done 