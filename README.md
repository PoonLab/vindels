# vindels
Developing an empirical model of sequence insertion and deletion in virus genomes

1) collected gp120 sequence data from patient samples on LANL (60,000 sequences)
2) Filtered sequences by >1400 nt; subtype present ; date present (8957 sequences) [subtype_output.py]
3) Performed pairwise alignments on the gp120 reference sequence (HXB2) using MAFFT [pairwise.py]
4) Cut and output the pairwise-aligned Variable and Conserved Regions by reading in frame [V_region_slicing2.py]
5) Concatenated and output the conserved regions by subtype [V_region_slicing3.py]
6) Multiple sequence alignment in MAFFT by subtype [MSA.py]

      
