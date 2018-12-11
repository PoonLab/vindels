import re
from glob import glob
import os
from gotoh2 import *
from seqUtils import *



def get_boundaries(str):
    gap_prefix = re.compile('^[-]+')
    gap_suffix = re.compile('[-]+$')

    # returns indices of without gap prefix and suffix
    res = [0, len(str)]
    left = gap_prefix.findall(str)
    right = gap_suffix.findall(str)
    if left:
        res[0] = len(left[0])

    if right:
        res[1] = len(str) - len(right[0])

    return res


gp120 = open("/home/jpalmer/PycharmProjects/hiv-withinhost/gp120.txt",'r')
ntref = ""
for line in gp120:
    line = line.strip("\n")
    ntref += line
aaref = translate_nuc(ntref, 0)

print(ntref)

folder = glob("/home/jpalmer/PycharmProjects/hiv-withinhost/2FilteredSeqs/*.fasta")

for file in folder:

    input = open(file, "r")

    data = parse_fasta(input)

    filename = file.split("/")[-1]

    output = open("/home/jpalmer/PycharmProjects/hiv-withinhost/3_1_pairwise/"+filename, 'w')

    unequal = []
    for header in data:

        nt_pair = Aligner()
        nt_pair.set_model('HYPHY_NUC')
        nt_pair.is_global = False
        nt_pair.gap_open_penalty = 30
        nt_pair.gap_extend_penalty = 10

        result = nt_pair.align(ntref, data[header])

        #use the reference sequence to cut the query sequence down
        left, right = get_boundaries(result[0])
        ntqry = result[1][left:right].replace("-","")

        result2 = nt_pair.align(ntref, ntqry)

        if "*" in ntqry:
            unequal.append(header)
            print(header)
            print(result2[0])
            print(result2[1])
            #continue



        output.write(">" + header + '\n')
        output.write(">ref\n" + result2[0] + "\n>query\n" + result2[1] + '\n')


    print(unequal)





