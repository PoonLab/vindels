from seqUtils import *
from glob import glob


#listed sequences hinder the accurate formation of alignments
#delete the sequences in the given lists, and re-print the Conserved Region .fasta file


file = open('/home/jpalme56/PycharmProjects/hiv-evolution-master/ConservedRegions4/B_CRegions.fasta', 'r')


blacklistC = ['AY772698','MF373131', 'KU319530', 'DQ275642', 'JN681248', 'KF725932', 'AB254149','AY878069']

blacklistB =['AF042101','AF042100', 'DQ295193','KY658703','MF373191','JF932483','AB485638','KF384814','KC935958',
             'KC935959','KF384811','KF384812','KU168257','KP109512','KY778594','JF932486','AY561239','KX960974',
             'KR182195','KR182182']

data = parse_fasta(file)

for i in data.keys():
    for x in blacklistB:
        if x in i:
            print(i)
            del data[i]

output_file = open('/home/jpalme56/PycharmProjects/hiv-evolution-master/ConservedRegions4/B_CRegions_edit.fasta', 'w')

for n in data.keys():
    output_file.write(">"+ n)
    output_file.write("\n")
    output_file.write(data[n])
    output_file.write("\n")

output_file.close()
