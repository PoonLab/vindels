from glob import glob
import regex
import os
from seqUtils import *


# used the Perl file provided by Abayomi to convert fasta files to phylip format

list = []

for file in glob("/home/jpalme56/PycharmProjects/hiv-evolution-master/MSAlignments5/*.fasta"):

    name = file.split(".fasta")

    list.append(name[0])


for filename in list:
    path1 = filename + ".fasta"
    path2 = filename + ".phy"
    phylip = "/home/jpalme56/Downloads/fasta2phylip.pl"

    cmd = "perl "+ phylip + " " + path1 +" "+ path2
    print(filename)
    os.system(cmd)
