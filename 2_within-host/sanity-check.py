import re


gene = '''atgaaagtgaaggggatacaaatgaattgtcagcacttattgagatgggg
gactatgatcttgggattgataataatctgtagtgctgcgaacaacttgt
gggttactgtctactatggggtacctgtgtggagagatgcagagaccacc
ttattttgtgcatcagatgctaaagcatatgtcacagaaaggcataatgt
ctgggctacacatgcttgtgtacccacagaccccaacccacaagaaatac
atttggcaaatgtgacagaaaagtttgacatgtggaaaaacagcatggta
gagcagatgcatacagatataatcagtctatgggacgaaaacctaaagcc
atgtgtaaaattaacccctctttgcattactttaaattgtactaacatca
tca----caa-----ata-agaacgccacaggggggaacctcacagagga
gggcaaggaagaattaagaaactgctctttcaatgcaaccacagaactaa
aggataaggaacaaaaagtacattcacttttttatagacttgatctagtt
gaacttaatgag---------ggt--------------------------
----aatagtagagatagt------a---------a------taatagta
tgtatagattaataaattgtaatacctcagcaattacacaggcttgtcca
aaggtatcctttgagccaattcccatacattattgtgccccagctggttt
tgcgattctaaagtgtagggagaaggagtttaatggaacagggccatgca
caaatgtcagcacagtacaatgcacacatgggatcaagccagtagtatca
actcagctgctgttaaatggcagtctagcagaagaaagggtacaaattag
atgtgaaaatatctcaaacaatgccaaaaccatactagtccaacttacta
tgcctgtgaaaattaattgtaccaggcctaacaacaatacaagaaaaagt
atacgtataggaccaggacaatcgttctatgcaacaggtgatataatagg
ggatataagaaaagcacattgtaatgtttcagaatcagaatggctggaag
ctttagggaaggtagctgaacaattaagaagacactttaata---ataaa
acaataacctttaatagctcctcaggaggggatttagaaatcacaacaca
tagttttaattgtggaggagaatttttctattgcaatacatcaagcctgt
ttaatagtactt------------ggaatagcactactaacagcacgcag
gagtcaaatgaaactataactctcccatgcaggataaagcaaattataaa
tatgtggcagagaacagggcaagcaatgtatgcccctcccatcccaggaa
aaataaggtgtgactcaaacatcacaggactaatattaacaagag-----
-------atggtgg------------------------aat---------
---taataacaat------ga---------cagtgaaacctttagacctg
gaggaggagacatgaggaacaattggagaagtgagttatataagtataaa
gtagtaaagattaacccactaggagtagcacccaccagggcaaagagaag
agtggtggagagagaaaaaaga'''
anc = '''atgaaagtgaaggggatacaaatgaattgtcagcacttattgagatgggg
gactatgatcttgggattgataataatctgtagtgctgcaaacaacttgt
gggttactgtctactatggggtacctgtgtggagagatgcagagaccacc
ttattttgtgcatcagatgctaaagcatatgtcacagaaaggcataatgt
ctgggctacacatgcttgtgtacccacagaccccagcccacaagaaatac
atttggcaaatgtgacagaaaagtttgacatgtggaaaaacagcatggta
gagcagatgcatacagatataatcagtctatgggacgaaagcctaaagcc
atgtgtaaaattaacccctctttgcattactttaaattgtactaacatca
tca----caa-----ata-agaacgccacaggggggaacctcacagagga
gggcaaggaagaattaagaaactgctctttcaatgcaaccacagaactaa
aggataaggaacaaaaagtacattcacttttttatagacttgatctagta
gaacttaatgag---------agt--------------------------
----aatagtagaaatagt------a---------a------tactagta
tgtatagattaataaattgtaatacctcagcaattacacaggcttgtcca
aaggtatcctttgagccaattcccatacattattgtgccccagctggttt
tgcgattctaaagtgtagggagaaggagtttaatggaacagggccatgca
caaatgtcagcacagtacaatgcacacatggcatcaggccagtagtatca
actcagctgctgttaaatggcagtctagcagaagaaaaggtacaaattag
atgtgaaaatatctcaaacaatgccaaaaccatactagtccaacttacta
cgcctgtgaaaattaattgtaccaggcctaacaacaatacaagaaaaagt
atacgtataggaccaggacaatcgttctatgcaacaggtgacataatagg
ggatataagaaaagcacattgtaatgtttcagaatcagaatggaaggaag
ctttagggaaggtagttgaacaattaagaagacactttaata---ataaa
acaataacctttaatagctcctcaggaggggatttagaaatcacaacaca
tagttttaattgtggaggagaatttttctattgcaatacatcaagcctgt
ttaatagtacttg------gaatgggaatagcactattaacagcacgcag
gagtcaaatgaaactataactctcccatgcaggataaagcaaattataaa
tatgtggcagagaacagggcaagcaatgtatgcccctcccatcccaggaa
aaataaggtgtgactcaaacatcacaggactaatattaacaagag-----
-------atggtgg------------------------aat---------
---taataacaat------ga---------cagtgaaacctttagacctg
gaggaggagacatgaggaacaattggagaagtgagttatataagtataaa
gtagtaaagattgacccactaggagtagcacccaccagggcaaagagaag
agtggtggagagagaaaaaaga'''


gene = gene.replace("\n","")
anc = anc.replace("\n","")
gene2 = gene.replace("-","")
anc2 = anc.replace("-","")
print(anc)
print("")
print(gene)

print(len(anc))
print(len(gene))

gidx = {}
ai = 0
for n, i in enumerate(gene):
    if i != "-":
        gidx.update({ai:n})
        ai += 1

aidx = {}
ai = 0
for n, i in enumerate(gene):
    if i != "-":
        aidx.update({ai:n})
        ai += 1

print(gene[gidx[540]:gidx[560]])
print(anc[aidx[540]:aidx[560]])

#gv5 = gene2.find(v5.lower())
#av5 = anc2.find(v5.lower())

#print(gene2[gv5-12:gv5+len(v5)])
#print(anc2[av5-12:av5+len(v5)])



#print(idx)
#print(gene.replace("\n","").replace("-",""))

#print(gene[:idx[574]])