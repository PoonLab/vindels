## Notes 


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


3) Apply MAFFT software to all full-length gp120 sequences within each patient to create within-host MSAs (4MSA_withinhost.sh)

4) Perform a hypermutation screen on all patient data sets (4_1_hm-screen)

    - IMPORTANT: patient 30651 was NOT submit to this hypermutation screen
        - we are choosing to NOT use the hypermutation screen for this one patient data set because there is a second T/F virus that essentially TAKES OVER the population
        - hypermutation essentially treats all sequences found at later timepoints to be hypermutants, which is problematic for our analysis

5) Perform a manual screen for recombination (111848)
    - looked at sequences displayed in RDP4
    - looked at the phylogenetic tree of superinfection patient 111848 
        - removed many of the most divergent branches on the tree; the ones found along the longest middle branch length 
    - split the patient data set into two 
    

4) Generate time-scaled phylogenetic trees using Bayesian Inference (implemented using BEAST) (5_Beauti.py, 6)

    1. Substitution model: TN93
    2. Site heterogeneity model: Gamma
    3. Clock model: Uncorrelated relaxed clock
    4. Coalescent model: Bayesian Skyline
    5. MCMC: 10e8 iterations

5) Use ancestral reconstruction software (Historian) to extract insertion and deletions from all tips (6_sample-beast-trees.py, 7_tree_mod.r, 8_historian.sh, 9_ancestors.py)



6) Estimate the rates, lengths, and compositions of indels in each gp120 variable region using past pipeline


FILTERING 

- 30651 performed poorly with 