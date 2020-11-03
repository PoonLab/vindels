#!/bin/bash
inDir=$1

if [[ $inDir != */ ]];
then
    inDir="$inDir/"
fi

for filename in $inDir*.xml; do
    beast $filename
done 
