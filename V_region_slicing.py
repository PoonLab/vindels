import os
from glob import glob
import re
import sys

def parse_fasta2(handle):
    # Modified parse fasta to return a dictionary of lists containing the reference [0] and the query [1]
    res = {}
    nt = ''
    h = ''
    for i in handle:
        if i[0] == '$': # skip h info
            continue
        elif i[0] == '>' or i[0] == '#':
            if ">ref" in i:
                continue
            elif ">query" in i:                    #Append the reference sequence
                res.update({h:[]})
                res[h].append(nt)
                nt = ''                             # reset containers
            elif len(nt) > 0 and ">query" not in i:             #Append the query sequence // Occurs when reaching the header of a new sequence
                res[h].append(nt)
                nt = ''
                #Parse the new header
                h = i.strip('\n')[1:]

            else:
                h = i.strip('\n')[1:]         #Occurs on the first header
        else:
            nt += i.strip('\n').upper()

    res[h].append(nt)
    return res

    

#GP120 Reference sequence file
ref_file = open("data/hxb2_gp120_sequence.txt", 'r')

gp120 = ''
for i in ref_file:
    gp120 += i.strip('\n').upper()

#Variable regions
v1 = gp120[390:469]
v2 = gp120[469:588]
v3 = gp120[885:993]
v4 = gp120[1152:1254]
v5 = gp120[1377:1410]

v_region = {0:v1,1:v2,2:v3,3:v4,4:v5}
regions = [(0,390), (390,469), (469,588), (885,993), (1152,1254), (1377,1410)]

print(v1, len(v1))
print(v2, len(v2))
print(v3, len(v3))
print(v4, len(v4))
print(v5, len(v5))

'''
#Single file testing
fasta_in = open("/home/jpalme56/PycharmProjects/hiv-evolution-master/Alignments/CD_pairwise.fasta", 'r')
data = parse_fasta(fasta_in)
output = open("/home/jpalme56/PycharmProjects/hiv-evolution-master/VRegions/V1-sequences.txt",'w')
'''


alignments = glob('data/Alignments/*.fasta')
print(alignments)



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

    with open(file) as handle:
        data = parse_fasta2(handle)

    for header, seq in data.items():
        #Extract the reference and query sequences
        #ref = data[i][0]
        #query = data[i][1]
        ref, query = seq
        left, right = get_boundaries(ref)
        r120 = ref[left:right]  # cut down to gp120
        q120 = query[left:right]
        
        # generate map from alignment to reference coordinates
        ri = 0  # reference index
        index = {}
        #Scan the reference for the variable region location
        for ai, x in enumerate(r120):
            # ai = alignment index
            if x != '-':
                # otherwise alignment has a gap, do not increment reference index
                index.update({ri: ai})
                ri += 1
        
        for left, right in regions:
            print(index[left])
            print(index[right])
            print(q120[index[left]:index[right]])
            print(gp120[left:right])
        sys.exit()
            




