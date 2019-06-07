#!/bin/bash
dir=/home/jpalmer/PycharmProjects/hiv-withinhost/5Trees/
for filename in /home/jpalmer/PycharmProjects/hiv-withinhost/4MSA/*.fastamsa; do
	name="$(cut -d'.' -f1 <<<"$(cut -d'/' -f7 <<<"$filename")")"
	fasttree -nt $filename > "$dir$name.tree"
done

