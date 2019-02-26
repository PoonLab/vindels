import os
from glob import glob
import re
import sys
from seqUtils import *

#GP120 Reference sequence file
rfile = open("./Phylo_Indels/hxb2_gp120_sequence.txt", 'r')

gp120 = ''
for i in rfile:
    gp120 += i.strip('\n').upper()

#Conserved regions
c_regions = [(0,390),  (588,885) , (993, 1152), (1254, 1377), (1410, 1532)]

#normal v regions
#v_regions = [(390, 468), (468, 588), (885, 993), (1152, 1254), (1377, 1410)]
#v regions with modified V5
v_regions = [(390, 468), (468, 588), (885, 993), (1152, 1254), (1368, 1419)]

#input folder of pairwise alignments from 2_pairwise_codon.py 
alignments = glob('/home/jpalmer/PycharmProjects/hiv-evolution-master/2_AAPairwise/*.fasta')

gap_prefix = re.compile('^[-]+')
gap_suffix = re.compile('[-]+$')

def get_boundaries (str):
    # return a tuple giving indices of subsequence without gap prefix and suffix
    res = [0,len(str)]
    left = gap_prefix.findall(str)
    right = gap_suffix.findall(str)
    if left:
        res[0] = len(left[0])

    if right:
        res[1] = len(str) - len(right[0])

    return res

for infile in alignments:
    count = 0
    print(infile)

    with open(infile) as handle:
        data = parse_fasta2(handle)

    subtype = os.path.basename(infile).split(".")[0]
    print(subtype)
    incorrect =[]

    outputv = open("/home/jpalmer/PycharmProjects/hiv-evolution-master/3_RegionSequences/VRegions-pre/"+ subtype + ".csv", "w")
    outputc = open("/home/jpalmer/PycharmProjects/hiv-evolution-master/4_Conserved/" + subtype + ".fasta", "w")

    for header, seq in data.items():
        #Extract the reference and query sequences
        ref, query = seq

        #Extract the accession number (last in the header)
        accno = header[-8:]
        if "." in accno:
            accno = accno.split(".")[1]

        ri = 0  
        index = {} # generate an alignment index
        for ai, x in enumerate(ref):
            if x != '-':
                index.update({ri: ai})
                ri += 1

        #Output
        #_____________________________________________________

        outputv.write(accno + ",")
        outputc.write(">" + header + "\n")

        for n1, n2 in v_regions:
            seq= query[index[n1]:index[n2]]
            gapcount = float(seq.count("-")) / len(seq)
            qcount = float(seq.count("?")) / len(seq)
            if n1 != 1368:
                outputv.write(seq.replace("-","") + ",")  # commas separators to make csv files
            else:
                outputv.write(seq.replace("-","") + "\n")

        for n1, n2 in c_regions:
            seq = query[index[n1]:index[n2]]
            if n1 != 1410:
                outputc.write(seq.replace("-",""))    # no separators so that the conserved regions are concatenated together
            else:
                outputc.write(seq.replace("-","") + "\n")
        
    print(count)
    print(len(incorrect))
