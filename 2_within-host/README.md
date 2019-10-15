## Within Host Analysis

Pipeline of Python, Bash, and R scripts used to perform a within-host phylogenetic analysis of gp120 insertions and deletions using longitudinally-sampled HIV-1 sequence data.





Methodology:
1) Downloaded gp120 sequence data from LANL from within-host HIV-1 sequence data (1_SequenceProcessing.py)

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


2) Perform pairwise alignment of each patient sequence with HXB2 (reference) (3Extraction.py) 
    
    1. to remove any flanking nucleotides outside of the gp120 gene   

    2. to extract the Variable Loops 1-5 for later analysis


3) Apply MAFFT software to all full length gp120 sequences within each patient to create within-host MSAs (4MSA_withinhost.sh)

4) Generate time-scaled phylogenetic trees using Bayesian Inference (implemented using BEAST) (5_Beauti.py, 6)

    1. Substitution model: TN93
    2. Site heterogeneity model: Gamma
    3. Clock model: Uncorrelated relaxed clock
    4. Coalescent model: Bayesian Skyline
    5. MCMC: 10e8 iterations

5) Use ancestral reconstruction software (Historian) to extract insertion and deletions from all tips (6_sample-beast-trees.py, 7_tree_mod.r, 8_historian.sh, 9_ancestors.py)



6) Estimate the rates, lengths, and compositions of indels in each gp120 variable region using past pipeline







