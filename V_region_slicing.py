import os

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
ref_file = open("hxb2_gp120_sequence.txt", 'r')

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

'''
#File testing
print(v1, len(v1))
print(v2, len(v2))
print(v3, len(v3))
print(v4, len(v4))
print(v5, len(v5))

fasta_in = open("/home/jpalme56/PycharmProjects/hiv-evolution-master/Alignments/CD_pairwise.fasta", 'r')
data = parse_fasta(fasta_in)
output = open("/home/jpalme56/PycharmProjects/hiv-evolution-master/VRegions/V1-sequences.txt",'w')
'''


#Slice one variable region at a time
for t in v_region:
    output = open("/home/jpalme56/PycharmProjects/hiv-evolution-master/VRegions/V" + str(t + 1) + "-sequences.txt",'w')

    #Read and parse all subtype alignments
    for filename in os.listdir("/home/jpalme56/PycharmProjects/hiv-evolution-master/Alignments"):

        fasta_in = open("/home/jpalme56/PycharmProjects/hiv-evolution-master/Alignments/" + filename, 'r')

        data = parse_fasta2(fasta_in)

        for i in data:
            #Extract the reference and query sequences
            ref = data[i][0]
            query = data[i][1]

            nt = ''
            once = 0

            #Scan the reference for the variable region location
            for n, x in enumerate(ref):

                if x.isalpha():
                    nt += x

                    #Record the end and start locations of the V region (one iteration)
                    if v_region[t] in nt and once == 0:
                        once = 1
                        end = n + 1

                        back_count = 0
                        nt_count = 0
                        #Count  backwards on the query to find the start position
                        for a in reversed(query[0:end]):
                            back_count += 1
                            if a.isalpha():
                                nt_count += 1
                            if nt_count == len(v_region[t]):
                                break

            output.write(query[end-back_count:end])
            output.write("\n")
