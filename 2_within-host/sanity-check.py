import re


gene = '''atgagagtgatggggatacagaggaattatccacaatggtggatatggag
catgttaggcttttggatgctaatgatttgtaatgggatgtgggtcacag
tctactatggggtacctgtgtggaaagaagcaaaaactactctattttgt
gcatc---------------------------------------------
--------------------------------------------------
-------agatgctaaagcatatgagaaagaagtgcataatgtctgggct
acacatgcctgtgtacccacagaccccaatccacaagaaatggttttaaa
aaatgtaacagaaaatttcaacatgtggaaaaatgacatggtggatcaga
tgcatgaagacgtaattagtttatgggatcaaagcctcaagccatgtgta
aagttgaccccactctgtgtcactctaaactgtaccaatgctactac---
--------------------------------------------------
--------------------------------------------------
--------------------------------------------------
--------------------------------------------------
------------------c------a------------------------
----a------------------------------------------t--
-------------------gc--------------------tact-----
----------------------------g---c---------------ca
------------------------------g---c---------------
---------a------------------a---t-----------------
-------agc------------------------------------agta
------ta------------------------a---t-------------
-----------------------------------agagggaatgaaaaa
ttgctctttcaatataaccacagaattaagagataagagagagaaaaaga
atgcacttttttataaacttgatatagtacaactag------------at
------------ggca------actctagtcagtatagattaataaattg
taatacctcagtcataacacaagcctgtccaaaggtctcttttgacccaa
ttcctatacattattgtgctccagctggttatgcgattctaaagtgtaat
aataagacattcactggaacaggaccgtgtaataatgtcagcacagtaca
atgtacacatggaattaagccagtggtttcaactcaactattgttaaatg
gtagcctagcagaaggagagataataattagatctgaaaatataacaaac
aatgccaaaacaataatagtacatctcaatgaatctgtaaagattgagtg
tacgagacccagtaataaaacaagaacaagtataagaataggaccaggac
aagcattttatgcaacaggacaagtaataggagacataagagaagcatat
tgtaacattagtgaaagtaaatggaatgaaactttacaaagggtaagtaa
aaaattaaaagaatacttccctcataagaatataacatttcaaccatcct
caggaggggacctagaaattacaacacatagctttaattgtggaggagaa
tttttctattgcaatacatcaagcctgtttaataggacat----------
--------at------------------atggctaatagtacagatatgg
ctaatagtacagaaactaacagtacacgaaccatcacaatccgctgcaga
ataaaacaaattataaacatgtggcaggaggtgggacgagcaatgtatgc
ccctcccattgcaggaaacataacatgtatatcaaatatcacaggactac
tattg------------------------------------------a--
-------------------ca---------------------aggga---
------------tggagg------------------a---aaa---acca
--------------------------------------------------
----------------------------------atacgga------ga-
--------ca---------------------ttcagacctggaggaggaa
atatgaaggacaattggagaagtgaattatataaatataaagtggtagaa
gttaagccattaggagtagcacccactaatgcaagaaggagagtggtgga
gagagaaaaaaga'''
anc = '''atgagagtgatggggatacagaggaattatccacaatggtggatatggag
catgttaggcttttggatgctaatgatttgtaatgggatgtgggtcacag
tctactatggggtacctgtgtggaaagaagcaaaaactactctattttgt
gcatc---------------------------------------------
--------------------------------------------------
-------agatgctaaagcatatgagaaagaagtgcataatgtctgggct
acacatgcctgtgtacccacagaccccaatccacaagaaatggttttaaa
aaatgtaacagaaaatttcaacatgtggaaaaatgacatggtggatcaga
tgcatgaagacgtaattagtttatgggatcaaagcctcaagccatgtgta
aagttgaccccactctgtgtcactctaaactgtaccaatgctactac---
--------------------------------------------------
--------------------------------------------------
--------------------------------------------------
--------------------------------------------------
------------------c------a------------------------
----a------------------------------------------t--
-------------------gc--------------------tact-----
----------------------------g---c---------------ca
------------------------------g---c---------------
---------a------------------a---t-----------------
-------agc------------------------------------agta
------ta------------------------a---t-------------
-----------------------------------agagggaatgaaaaa
ttgctctttcaatataaccacagaattaagagataagagagagaaaaaga
atgcacttttttataaacttgatatagtacaactag------------at
------------ggca------actctagtcagtatagattaataaattg
taatacctcagtcataacacaagcctgtccaaaggtctcttttgacccaa
ttcctatacattattgtgctccagctggttatgcgattctaaagtgtaat
aataagacattcactggaacaggaccgtgtaataatgtcagcacagtaca
atgtacacatggaattaagccagtggtttcaactcaactattgttaaatg
gtagcctagcagaaggagagataataattagatctgaaaatataacaaac
aatgccaaaacaataatagtacatctcaatgaatctgtaaagattgagtg
tacgagacccagtaataaaacaagaacaagtataagaataggaccaggac
aagcattttatgcaacaggacaagtaataggagacataagagaagcatat
tgtaacattagtgaaagtaaatggaatgaaactttacaaagggtaagtaa
aaaattaaaagaatacttccctcataagaatataacatttcaaccatcct
caggaggggacctagaaattacaacacatagctttaattgtggaggagaa
tttttctattgcaatacatcaagcctgtttaataggacat----------
--------at------------------atggctaatagtacagatatgg
ctaatagtacagaaactaacagtacacgaaccatcacaatccgctgcaga
ataaaacaaattataaacatgtggcaggaggtgggacgagcaatgtatgc
ccctcccattgcaggaaacataacatgtatatcaaatatcacaggactac
tattg------------------------------------------a--
-------------------ca---------------------aggga---
------------tggagg------------------a---aaa---aaca
---------------------a------t---------------a-----
-cg------------------------g------atacgga------gaplot(tre)-
--------ca---------------------ttcagacctggaggaggaa
atatgaaggacaattggagaagtgaattatataaatataaagtggtagaa
gttaagccattaggagtagcacccactaatgcaagaaggagagtggtgga
gagagaaaaaaga'''

#TAGTACATACATGCA
#TGCAATACATCAAAACTGTTTAATAGTACATACATGCATAGCACATACATGAATAATCATACAAACAGTAATAGTACAAGTAATCCAAATTCGACCATCACCATCCAATGC

#TGCAATACATCAAAACTGTTTAATAGTACATACATGCATAGCACATACATGAATAATCATACAAACAGTAATAGTACAAGTAATCCAAATTCGACCATCACCATCCAATGC
#TGCAATACATCAAAACTGTTTAA               TAGCACATACATGAATAATCATACAAACAGTAATAGTACAAGTAATCCAAATTCGACCATCACCATCCAATGC

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
for n, i in enumerate(anc):
    if i != "-":
        aidx.update({ai:n})
        ai += 1

print(gene[gidx[1348]:gidx[1354]])
print(anc[gidx[1348]:gidx[1354]])


#gv5 = gene2.find(v5.lower())
#av5 = anc2.find(v5.lower())
#TGTACCGATGCTAATGCTACTGCCAGCAATAGCAGTATAATAAAGGGAATGAATAGCAGTATGATAGAGGAAATGAAAAAT

#print(gene2[gv5-12:gv5+len(v5)])
#print(anc2[av5-12:av5+len(v5)])
#TGTAGTCATAAGGTTAATGCCATCAATGGGAGTATTAGCATCAATGGGACGATAAAGGAAGGAATGAAAAAT

'''
> ins[86,]
                                  Accno Vloop Vlength Subtype Count    Seq Pos                                                                              Vseq     Pat Vpos before after
19110 C.MW.2009.56552.KC247387.24.3_673     1      81       C     1 GATGCT   6 TGTACCGATGCTAATGCTACTGCCAGCAATAGCAGTATAATAAAGGGAATGAATAGCAGTATGATAGAGGAAATGAAAAAT 56552-a  384     NA     0
> ins[87,]
                                  Accno Vloop Vlength Subtype Count                      Seq Pos                                                                              Vseq     Pat
191.1 C.MW.2009.56552.KC247387.24.3_673     1      81       C     1 AAAGGGAATGAATAGCAGTATGAT  41 TGTACCGATGCTAATGCTACTGCCAGCAATAGCAGTATAATAAAGGGAATGAATAGCAGTATGATAGAGGAAATGAAAAAT 56552-a
      Vpos before after
191.1  419     NA     0

TGTACCGATGCTAATGCTACTGCCAGCAATAGCAGTATAATAAAGGGAATGAATAGCAGTATGATAGAGGAAATGAAAAAT
TGTACC      AATGCTACTGCCAGCAATAGCAGTATAAT                        AGAGGAAATGAAAAAT
GAATG
TGCAATACATCAAGCCTGTTTAATAGTACTTGGAATGGGAATAGCACTACTAACAGCACGCAGGAGTCAAATGAAACTATAACTCTCCCATGC
'''

#print(idx)
#print(gene.replace("\n","").replace("-",""))

#print(gene[:idx[574]])