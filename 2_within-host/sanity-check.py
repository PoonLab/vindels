import re


gene = '''atgagagtgatggggatacagaggaattatccacaatggtggatatggag
catgttaggcttgtggatgctaatgatttgtaatgggatgtgggtcacag
tctactatggggtacctgtgtggaaaaaagcaaaaactactctattttgt
gcatc---------------------------------------------
--------------------------------------------------
-------agatgctaaagcatatgagaaagaagtgcatagtgtctgggct
acacatgcctgtgtacccacagaccccaatccacaagaaatggttttaaa
aaatgtaacagaaaatttcaacatgtggaaaaatgacatggtagatcaga
tgcatgaagatgtaattagtttatgggatcaaagcctcaagccatgtgta
aagttgaccccactctgtgtcactctaaactg------------taccga
tgct------aa---------t---------------------------g
c------------------tactgccag------------caa-------
-----t------agc------agtata------at---------------
---------------------a---------------aagg---------
------------------g------a------------------------
----a------------------------------------------t--
-------------------g------------------------------
--------------------------------------------------
--------------------------------------------------
---------a------------------a---t-----------------
-------agc------------------------------------agta
------tg------------------------a---t-------------
-----------------------------------agaggaaatgaaaaa
ttgctctttcaatataaccacagaattaagagataagagagagaaaaaga
atgcacttttttataaacttgatatagtacaactag------------at
------------ggca------actctagtcagtatagattaataaattg
taatacctcagtcataacacaagcctgtccaaaggtctcttttgacccaa
ttcctatacattattgtgctccagctggttatgcgattctaaagtgtaat
aataagacattcaatggaacaggaccgtgtaataatgtcagcacagtaca
atgtacacatggaattaagccagtggtttcaactcaactattgttaaatg
gcagcctagcagaaggagagataataattagatctgaaaatataacaaac
agtgccaaaacaataatagtacatctcaatgaatctgtaaagattgagtg
tacgagacccagtaataacacaagaacaagtataagaataggaccaggac
aagcattttatgcaacaggacaagtaataggagacataagaaaagcacat
tgtaacattagtgaaagtaaatggaatgaaactttacaaagggtaagtaa
aaaattaaaagaatacttccctgataggaatataacatttcaaccatcct
caggaggggacccagaaattacaacacatagctttaattgtggaggaaaa
tttttctattgcaatacatcaagcctgtttaatagaacat----------
--------at------------------atggctaatagtacagatatgg
ctaatagtacagaaactaacagtacacgaatcatcacaatccgctgcaga
ataaaacaaattataaacatgtggcaggaggtgggacgagcaatgtatgc
ccctcccattgcaggaaacataacatgtatatcaaatatcacaggactac
tattg------------------------------------------a--
-------------------ca---------------------aggga---
------------tggagg------------------a---aat---aaca
---------------------a------t---------------a-----
-cg------------------------g------------a------ga-
--------ca---------------------ttcagacctgtaggaggaa
atatgaaggacaattggagaagtaaattatataaatataaagtggtagaa
gttaagccattaggagtagcacccactaaggcaagaaggagaatggtgga
gagagaaaaaaga'''
anc = '''atgagagtgatggggatacagaggaattatccacaatggtggatatggag
catgttaggcttgtggatgctaatgatttgtaatgggatgtgggtcacag
tctactatggggtacctgtgtggaaagaagcaaaaactactctattttgt
gcatc---------------------------------------------
--------------------------------------------------
-------agatgctaaagcatatgagaaagaagtgcataatgtctgggct
acacatgcctgtgtacccacagaccccaatccacaagaaatggttttaaa
aaatgtaacagaaaatttcaacatgtggaaaaatgacatggtagatcaga
tgcatgaagatgtaattagtttatgggatcaaagcctcaagccatgtgta
aagttgaccccactctgtgtcactctaaactg------------taccga
tgct------aa---------t---------------------------g
c------------------tactgccag------------caa-------
-----t-------gc----------ta------ct---------------
---------------------g---------------ccag---------
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
aataagacattcaatggaacaggaccgtgtaataatgtcagcacagtaca
atgtacacatggaattaagccagtggtttcaactcaactattgttaaatg
gtagcctagcagaaggagagataataattagatctgaaaatataacaaac
agtgccaaaacaataatagtacatctcaatgaatctgtaaagattgagtg
tacgagacccagtaataacacaagaacaagtataagaataggaccaggac
aagcattttatgcaacaggacaagtaataggagacataagaaaagcacat
tgtaacattagtgaaagtaaatggaatgaaactttacaaagggtaagtaa
aaaattaaaagaatacttccctcataagaatataacatttcaaccatcct
caggaggggacctagaaattacaacacatagctttaattgtggaggagaa
tttttctattgcaatacatcaagcctgtttaataggacat----------
--------at------------------atggctaatagtacagatatgg
ctaatagtacagaaactaacagtacacgaatcatcacaatccgctgcaga
ataaaacaaattataaacatgtggcaggaggtgggacgagcaatgtatgc
ccctcccattgcaggaaacataacatgtatatcaaatatcacaggactac
tattg------------------------------------------a--
-------------------ca---------------------aggga---
------------tggagg------------------a---aat---aaca
---------------------a------t---------------a-----
-cg------------------------g------------a------ga-
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

print(gene[gidx[378]:gidx[460]])
print(anc[gidx[378]:gidx[460]])


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