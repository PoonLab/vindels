#!/bin/bash
for foldername in /home/jpalmer/PycharmProjects/hiv-withinhost/5BEAST/*; do
	for filename in $foldername/*.xml; do
		echo $filename
		beast -overwrite $filename
	done
done

