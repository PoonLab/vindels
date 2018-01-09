import sys
from glob import glob
from seqUtils import *



input = open("/home/jpalme56/PycharmProjects/hiv-evolution-master/5MSAlignments/C_MSA_edit.fasta",'r')

#list format
fasta = convert_fasta(input)

transposed = transpose_fasta(fasta)

pos = 0

for x in transposed:
    pos += 1

    length = len(x)
    gaps = x.count("-")

    freq = float(gaps)/length
    print(freq, pos)
    if freq < 0.2:

        print('Cut Point:'+ str(pos))
        break
print(pos)

input.close()

input = open("/home/jpalme56/PycharmProjects/hiv-evolution-master/5MSAlignments/C_MSA_edit.fasta",'r')
data = parse_fasta(input)

output = open("/home/jpalme56/PycharmProjects/hiv-evolution-master/5MSAlignments/C_MSA_shortened.fasta",'w')


#dictionary format


for i in data.keys():
    data[i] = data[i][pos-1:]

    output.write(">" + i + "\n")
    output.write(data[i] + "\n")

input.close()
