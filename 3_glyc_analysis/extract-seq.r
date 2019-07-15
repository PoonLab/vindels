# extract sequences 
require(ape)
tre <- read.tree("~/PycharmProjects/glyc-analysis/6_pruned/pruned.tree")

full.len <- read.csv("~/PycharmProjects/glyc-analysis/3_sequences/full/gp120.csv", stringsAsFactors = F)
# retrieve the sequences from the tree 


idx <- match(tre$tip.label, full.len$header)
# look up and retrieve their FULL LENGTH sequences 

fasta <- full.len[idx,1:2]

x <- as.vector(fasta$seq)
names(x) <- fasta$header
y <-as.DNAbin.DNAStringSet(x)

write.FASTA(y, "~/PycharmProjects/glyc-analysis/7_prnseq/gp120.fasta")
