import re
from glob import glob
import os
from seqUtils import *



folder = glob("/home/jpalmer/PycharmProjects/hiv-withinhost/3RegionSequences/full_length/*.fasta")
outdir = '/home/jpalmer/PycharmProjects/hiv-withinhost/3RegionSequences/full_length_accno/'
for infile in folder:
    with open(infile) as handle:
        data = parse_fasta(handle)

    filename = os.path.basename(infile)
    outfile = open(outdir+filename,'w')
    for header in data.keys():
        accno = header.split(".")[4]  # isolate just the accession number in the fifth position 

        outfile.write( '>' + accno +'\n'+ data[header] +'\n')
    