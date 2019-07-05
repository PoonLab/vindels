
gene = '''atgaaagtgaaggggatacaaatgaattgtcagcacttattgagatgggg
gactatgatcttgggattggtaataatctgtagtgctgcaaacaacttgt
gggttactgtctactatggggtacctgtgtggagagatgcagagaccacc
ttattttgtgcatcagatgctaaagcatatgtcacagaaagacataatgt
ctgggctacacatgcttgtgtacccacagaccccaacccacaagaaatac
atttggcaaatgtgacagaaaagtttgacatgtggaaaaacagcatggta
gagcagatgcatacagatataatcagtctatgggacgaaagcctaaagcc
atgtgtaaaattaacccctctttgcattactttaaattgtactaacatca
tca----caa-----ata-agaacaccacaggggggaacctcacagagga
gggcaaggaagaattaagaaactgctctttcaatgcaaccacagaactaa
aggataaggaacaaaaagtacattcacttttttatagacttgatctagta
gaacttaatgagg---------gt--------------------------
----aatggtaatagtagt------a---------a------tactagta
tgtatagattaataaattgtaatacctcagcaattacacaggcttgtcca
aaggtatcctttgagccaattcccatacattattgtgccccagctggttt
tgcgattctaaagtgtagggagaaggattttaatggaacagggctatgca
acaatgtcagcacagtacaatgcacacatgggatcaagccagtagtatca
actcagctgctgttaaatggcagtctagcagaaggaaaggtaaatattag
atgtgaaaatatctcaaacaatgccaaaaccatactagtccaacttacta
agcctgtgagaattaattgtaccaggcctaacaacaatacaagaaaaagt
atacgtataggaccaggacaatcgttctatgcaacaggtgacataatagg
ggatataagaaaagcacattgtaatgtttcaggatcagaatggatggaag
ctttaaggaatgtaagtgaacaattaagaaaacactttaata---ataaa
acaataacctttaatagctcctcaggaggggatttagaaatcacaacaca
tagttttaattgtggaggagaatttttctattgcaatacatcaagcctgt
ttaatagtacttg------gaatgggactagcattactaacagcacgcag
gagtcaaataaaaatataactctcccatgcaggataaagcaaattataaa
tatgtggcaaagaacagggcaagcaatgtatgcccctcccatcccaggaa
aaataaggtgtgactcaaacatcacaggactaatattaacaagag-----
-------atggtgg------------------------aat---------
---taataacaatga---------cagtgaaacctttagacctg
gaggaggagacatgaggaacaattggagaagtgagttatataagtataaa
gtagtaaagattgacccactaggagtagcacccaccagggcaaagagaag
agtggtggagagagaaaaaaga'''



idx = {}
ai = 0

for n, i in enumerate(gene):
    if i != "-":
        idx.update({ai:n})
        ai += 1

print(gene[:1549].replace("-",""))

print(gene[idx[1391-6]:idx[1391]])