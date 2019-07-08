
gene = '''atgagagtgatggggatattgaggaattgtcaacgatggtggatatgggg
catcttaggcttttggatgttaatgatttgtagtgggagcttgtgggtaa
cagtctattatggggtacctgtgtggaaagaagcaaagactactctattt
tgtgcatcagatgctaaagcatatgagagagaggtgcataatgtctgggc
tacacatgcctgtgtacccacagaccctgacccacaagaattagttttgg
aaaatgtaacagaaaattttaacatgtggaagaatgatatggtggatcag
atgcatgaggatataatcagtttatgggatgaaagcctaaagccatgcgt
aaagctgaccccactctgtgtcactttaaactg-----------------
-------ta------------------------------gtcataatgtt
accatca------at------------------a----------------
-----at------aa---tgata---------------------------
--------------------------------------------------
--------------------------------------------------
--------------------------------------------------
--------------------------------c-----------------
-t------------------a------------------ccaa-------
-----------t------------------g------------------g
gaa------------------tagtagca---------------------
------tc---------------a---atgg---------g---------
---------------a------------------------cgataaagga
aggaatgaaaaattgctctttcactatacccacagaactgaaagacaaga
caa---agaaggtatatacactttttaatgaacttgatgtagtgaaactt
------agtggaaatgacag------------------------------
----------------t------------------------g------ag
gatgagtacagattgatacattgtaatacctcagccataaaacaagcctg
tccaaagatctcttttgaaccaattcctatacacttttgtgctccagctg
gttatgcgattctaaagtgtaataataagacattcaatggaacaggacca
tgcaacaatgtcagcacagtacaatgtacacatggaattaagccagtagt
atcaactcaactactgttaaatggtagtctagcagaagaaaagataataa
ttagatctgaaaatataacaaataatgccaaaacaataatagtacacctt
aaagaccctgtagaaattgtgtgtacaaggcctggcaataatacaaggac
aagtgtga------------g------gataggaccaggacaaacattct
atgcaataggagacataataggagatataagaaaagcatattgtaacatt
agtgaagcaaaatggaatgaaactttaaagcaggtagctggagaactaca
aaaatactataacaca------aacacaaccataatctttaaccaaccct
caggaggggacctagaaattacaacacatagctttaattgtggaggagaa
tttttctattgcaatacatcaaaactgtttaa---------------tag
tacatacatgca------------------------------taatcata
caata------a------gtaata------gta---------ca------
---agt------------ag---------tt-------------------
------------------caaattcgaccatcaccatccactgcagaata
aaacaaattataaacatgtggcagggggtaggacaagcaatatatgcccc
tcccattaaaggaaacattacatgtagatcaaatatcacaggactactat
tgacacgtgatggaggtaataata---atactaatgagacattcagacct
ggaggaggagatatgagagacaattggagaagtgaattatataagtataa
agtagtagaaattaagccactaggagtagcacccaccaaggcgaaaagac
gagtggtggagaaagaaaaaaga'''

nogaps = '''atgagagtgatggggatattgaggaattgtcaacgatggtggatatggggcatcttaggcttttggatgttaatgatttgtagtgggagcttgtgggtaacagtctattatggggtacctgtgtggaaagaagcaaagactactctattttgtgcatcagatgctaaagcatatgagagagaggtgcataatgtctgggctacacatgcctgtgtacccacagaccctgacccacaagaattagttttggaaaatgtaacagaaaattttaacatgtggaagaatgatatggtggatcagatgcatgaggatataatcagtttatgggatgaaagcctaaagccatgcgtaaagctgaccccactctgtgtcactttaaactgtagtcataatgttaccatcaataataatgatactaccaatgggaatagtagcatcaatgggacgataaaggaaggaatgaaaaattgctctttcactatacccacagaactgaaagacaagacaaagaaggtatatacactttttaatgaacttgatgtagtgaaacttagtggaaatgacagtgaggatgagtacagattgatacattgtaatacctcagccataaaacaagcctgtccaaagatctcttttgaaccaattcctatacacttttgtgctccagctggttatgcgattctaaagtgtaataataagacattcaatggaacaggaccatgcaacaatgtcagcacagtacaatgtacacatggaattaagccagtagtatcaactcaactactgttaaatggtagtctagcagaagaaaagataataattagatctgaaaatataacaaataatgccaaaacaataatagtacaccttaaagaccctgtagaaattgtgtgtacaaggcctggcaataatacaaggacaagtgtgaggataggaccaggacaaacattctatgcaataggagacataataggagatataagaaaagcatattgtaacattagtgaagcaaaatggaatgaaactttaaagcaggtagctggagaactacaaaaatactataacacaaacacaaccataatctttaaccaaccctcaggaggggacctagaaattacaacacatagctttaattgtggaggagaatttttctattgcaatacatcaaaactgtttaatagtacatacatgcataatcatacaataagtaatagtacaagtagttcaaattcgaccatcaccatccactgcagaataaaacaaattataaacatgtggcagggggtaggacaagcaatatatgcccctcccattaaaggaaacattacatgtagatcaaatatcacaggactactattgacacgtgatggaggtaataataatactaatgagacattcagacctggaggaggagatatgagagacaattggagaagtgaattatataagtataaagtagtagaaattaagccactaggagtagcacccaccaaggcgaaaagacgagtggtggagaaagaaaaaaga'''

idx = {}
ai = 0

for n, i in enumerate(gene):
    if i != "-":
        idx.update({ai:n})
        ai += 1
#print(idx)
#print(gene.replace("\n","").replace("-",""))

print(gene[:idx[574]])