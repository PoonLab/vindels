import re
from glob import glob
import os
from seqUtils import *



folder = glob("/home/jpalmer/PycharmProjects/hiv-withinhost/4MSA/*.fasta")
outdir = '/home/jpalmer/PycharmProjects/hiv-withinhost/4_1Accno/'
for infile in folder:
    with open(infile) as handle:
       data = parse_fasta(handle)

    #newfilename = re.sub("_\d{3}.", ".", infile)
    #os.rename(infile,newfilename)
    filename = os.path.basename(infile)
    outfile = open(outdir+filename,'w')
    if filename != "VN_Data.fasta":
        for header in data.keys():
            accno = header.split(".")[4]  # isolate just the accession number in the fifth position 
            date = header.strip("\n").split("_")[1]

            outfile.write('>' + accno+"_"+date+'\n'+ data[header] +'\n')
    else:
        for header in data.keys():
            outfile.write('>' + header+'\n'+ data[header] +'\n')
    