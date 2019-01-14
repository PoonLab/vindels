from glob import glob
from beauti import parse
import os
import re

p = re.compile('.+_([0-9]+)\..+')
infile = open("template.xml")


template_file = 'template.xml'

#parse(template_file=args.template, fasta_file=args.fasta, stem=args.stem, outfile=args.out, nreps=args.nreps)
files = glob('/home/jpalmer/PycharmProjects/hiv-withinhost/4MSA/*.fasta')

for f in files:
    name = os.path.basename(f)
    path = os.path.dirname(f)
    #os.mkdir('/home/jpalmer/PycharmProjects/hiv-withinhost/5BEAST/'+name.split(".fasta")[0]+'/')
    stem = re.sub('\.fasta', '.xml', name)
    print(stem)

    ofn = '/home/jpalmer/PycharmProjects/hiv-withinhost/5BEAST/'+stem
    parse(template_file, f, stem, ofn, 'days', '1')
