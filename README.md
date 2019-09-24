# Variable Loop Indels
Project aimed at determining the rates of insertions and deletions in the five variable regions (V1-V5) of the HIV-1 gp120 surface envelope glycoprotein 

Overview:
1) parsed over 26,000 HIV-1 gp120 sequences from the Los Alamos National Laboratory (LANL) HIV Database and sorted them into their respective group M subtypes and circulating recombinant forms (CRFs)
2) filtered sequences to ensure sufficient coverage of gp120 (>1,400 nt) and availability of collection dates
3) performed a pairwise alignments between each sequence and the HXB2 reference genome to locate and extract the five variable and five conserved regions of gp120
4) performed multiple sequence alignments (MSAs) among concatenated conserved regions within each group M clade 
5) reconstructed phylogenetic trees from these MSAs, and rescaled the trees in time using sequence collection dates
6) extracted cherries of the phylogenetic trees and checked for length differences in their variable regions to detect indels
7) applied a binomial-Poisson model to these data to determine indel rates for each variable loop within each group M clade





# kaitlyns idea
- look at the flanking regions left and right of the indel event
- see if +/- 5 nucleotides might contain higher proportions of A or T 
- 