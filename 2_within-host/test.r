require(ape)


msa <- read.FASTA("~/PycharmProjects/hiv-withinhost/4MSA/12255.fasta")

seqs <- read.FASTA("~/PycharmProjects/hiv-withinhost/3RegionSequences/full_length/fullname/12255.fasta")

u <- unique(c(names(msa),names(seqs)))
c(names(msa),names(seqs))
