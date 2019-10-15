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
    
    * to remove any flanking nucleotides outside of the gp120 gene   
    * to extract the Variable Loops 1-5 for later analysis
    * to sort experiment-based files into 
    * LANL-derived sequence data was C1-C5 (full gp120) 
    * Vlad Novitsky's sequence data was V1-C5 (missing C1) and therefore, required different handling 
    * FILTER: removed redundant patient data sets (of the identical patient data sets, only the one containing the maximum number of sequences was chosen)
    * FILTER: removed patient data sets that did not contain at least 5 unique time points 

3) Apply MAFFT software to all full-length gp120 sequences within each patient to create within-host MSAs (4MSA_withinhost.sh)
    * `mafft --auto --reorder $filename` 
    * 
4) Perform a hypermutation screen on all patient data sets (4_1_hm-screen)

    - this step takes alignment file inputs from the `4MSA/prelim` folder, finds the 1st time point, generates a consensus sequence from this first timepoint, uses this consensus to detect hypermutants, and creates a filtered patient dataset from the original file in `3RegionSequences/full_length/hm-screen/` 
    - multiple sequence alignments were then re-generated using these screened files in `3RegionSequences/full_length/hm-screen/` to produce the final output in `/4MSA/hm-screen/`
    - SUPERINFECTION (111848)
        - this data set had its outlier removed and was split apart in `/4MSA/prelim/`, generating two separate alignment files that contained 447 and 616 sequences, respectively.
        - the hypermutation screen was performed on these two files separately using each of their 1st time points
        - the result was two filtered data sets in `3RegionSequences/full_length/hm-screen/` which were then combined to generate a combined, filtered data set
        - MSA was re-generated using the combined data set to produce a single result 
    - SPECIAL CASE (30651):
        - this patient was NOT submit to this hypermutation screen
        - we are choosing to NOT use the hypermutation screen for this one patient data set because there is a second T/F virus that essentially takes over the population
        - hypermutation essentially treats all sequences found at later timepoints to be hypermutants, which is problematic for our analysis
    

5) Perform a manual screen for recombination (111848)
    - looked at the phylogenetic tree of superinfection patient 111848 
        - there is a single long branch separating the two distinct populations within this patient 
        - made manual records of the sequences falling along this longest branch (and treated them as recombinants to be removed)
    - split the patient data set into two 
    - looked at sequences displayed in RDP4
        - manually looked up these sequences in the alignment to check whether they are in fact notable recombinants 
    

6) Generate time-scaled phylogenetic trees using Bayesian Inference (implemented using BEAST) (5_Beauti.py, 5_1_xml_edit.py)

    - testing multiple iterations of BEAST
        - constant size coalescent model has equal likelihood to the Bayesian Skygrid coalescent model in many cases
        - shows evidence of using the constant size instead of the Skygrid (due to it being the simpler model containing fewer parameters)
    - clock prior 
        - still need to test this to see whether this is beneficial 

5) Use ancestral reconstruction software (Historian) to extract insertion and deletions from all tips (6_sample-beast-trees.py, 7_tree_mod.r, 8_historian.sh, 9_ancestors.py)



6) Estimate the rates, lengths, and compositions of indels in each gp120 variable region using past pipeline


FILTERING 

- 30651 performed poorly with 