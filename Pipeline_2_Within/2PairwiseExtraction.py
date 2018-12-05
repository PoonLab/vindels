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

        left, right = get_boundaries(result[0])

        ntqry = result[1][left:right].replace("-","")

        aaqry = translate_nuc(ntqry,0)


        aa_pair = Aligner()
        aa_pair.set_model('EmpHIV25')
        aa_pair.gap_extend_penalty = 10
        aa_pair.gap_open_penalty = 30
        aa_pair.is_global = True

        result2 = aa_pair.align(aaref,aaqry)

        newRef = list(ntref)
        newQry = list(ntqry)

        # reads through the amino acid alignment and adds codon gaps to the proper locations
        for i in range(len(result2[0])):
            if result2[0][i] == '-':
                newRef[i * 3:i * 3] = ['-', '-', '-']

            if result2[1][i] == '-':
                newQry[i * 3:i * 3] = ['-', '-', '-']

        if len(newRef) != len(newQry):
            unequal.append(header)
            print(header)

            continue

        finalRef = "".join(newRef)
        finalQry = "".join(newQry)

        output.write(">" + header + '\n')
        output.write(">ref\n" + finalRef + "\n>query\n" + finalQry + '\n')


    print(unequal)




