## Within Host Analysis

Pipeline of Python, Bash, and R scripts used to perform a within-host phylogenetic analysis of gp120 insertions and deletions using longitudinally-sampled HIV-1 sequence data.





Methodology:
1) Downloaded gp120 sequence data from LANL from within-host HIV-1 sequence data (1SequenceProcessing.py)

    1. Parameters for LANL search (29 sets and 11350 sequences)

        * longitudinal only
        * known timeline
        * 100 seqs per patient minimum
        * HIV-1, any  project or subtype
        * genomic region: gp120 

    2. FASTA labels used

        1) subtype
        2) country
        3) sampling year
        4) patient ID 
        5) accession no 
        6) day 1st sample
        7) days since treatment start
        8) days of infection
        9) days since seroconversion
        10) number of timepoints 

    3. sequence filtering

        * length > 1400 nt 
        * at least one of the four day fields (6-9) has been populated
        * subtype and collection year present


2) Perform pairwise alignment of each patient sequence with HXB2 (reference) to remove any irrelevant nucleotides outside gp120 

3) Use full length gp120 sequences to create  a MSA between all sequences within each host 
4) Generate phylogenetic trees using Bayesian Inference (implemented using BEAST software package)
5) Apply an ancestral reconstruction method to extract insertion and deletion sequences from data 
6) Estimate the rates, lengths, and compositions of indels in each gp120 variable region using past pipeline




