import os
from glob import glob
import re
import sys
from seqUtils import *


#GP120 Reference sequence file
ref_file = open("hxb2_gp120_sequence.txt", 'r')

gp120 = ''
for i in ref_file:
    gp120 += i.strip('\n').upper()

#Variable regions
#print(len(v1), len(v2), len(v3), len(v4), len(v5))


c_regions = [(0,390),  (588,885) , (993, 1152), (1254, 1377), (1410, 1532)]
#v_regions = [(390, 468), (468, 588), (885, 993), (1152, 1254), (1377, 1410)]
#modified
v_regions = [(390, 468), (468, 588), (885, 993), (1152, 1254), (1377, 1410)]
alignments = glob('/home/jpalmer/PycharmProjects/hiv-evolution-master/2_AAPairwise/*++.fasta')

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



#Read and parse all subtype alignments
for file in alignments:
    count = 0
    print(file)

    with open(file) as handle:
        data = parse_fasta2(handle)


    filename = file.split("/")[-1].split("++")[0]

    incorrect =[]

    outputv = open("/home/jpalmer/PycharmProjects/hiv-evolution-master/3_RegionSequences/VRegions_mod/"+ filename + "_VR.csv", "w")
    outputc = open("/home/jpalmer/PycharmProjects/hiv-evolution-master/4_2_Conserved/" + filename + "_CR.fasta", "w")

    for header, seq in data.items():
        #Extract the reference and query sequences


        ref, query = seq

        accno = header[-8:]

        if "." in accno:
            accno = accno.split(".")[1]



        #Determines the start and end positions of the reference sequence in the alignment
        #REMOVED : was moved to the pairwise alignment step earlier #2
        #left, right = get_boundaries(ref)

        #Work only with relevant sections
        r120 = ref   # cut down to gp120
        q120 = query
        #print(header)


        # generate map from alignment to reference coordinates
        ri = 0  # reference index

        #{nucleotide #: actual position}
        index = {}

        #Scan the reference for the variable region location
        for ai, x in enumerate(r120):
            # ai = alignment index, literal position

            if x != '-':
                # otherwise alignment has a gap, do not increment reference index
                index.update({ri: ai})
                ri += 1

        #Output

        #FILTERING STEP
        #remove any sequences that contain a problematic variable region (more than 70% gaps)


        outputv.write(accno + ",")
        outputc.write(">" + header + "\n")


        for n1, n2 in v_regions:
            seq= q120[index[n1]:index[n2]]
            gapcount = float(seq.count("-")) / len(seq)
            qcount = float(seq.count("?")) / len(seq)


            if n1 != 1368:
                outputv.write(seq.replace("-","") + ",")  # commas separators to make csv files
            else:
                outputv.write(seq.replace("-","") + "\n")

        for n1, n2 in c_regions:
            seq = q120[index[n1]:index[n2]]
            if n1 != 1410:
                outputc.write(seq.replace("-",""))    # no separators so that the conserved regions are concatenated together
            else:
                outputc.write(seq.replace("-","") + "\n")
    print(count)
    print(len(incorrect))
