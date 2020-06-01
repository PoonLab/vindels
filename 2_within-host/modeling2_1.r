require(bbmle)
require(ape)
# used to fill in deletion gaps found in the tip sequences 
nt <- c("A", "C", "G", "T")
require(stringr)
source("~/vindels/2_within-host/utils.r")
source("~/vindels/2_within-host/slippage-model.r")


# ------  DATA IMPORT AND MANIPULATION -----

path <- "~/PycharmProjects/hiv-withinhost/"
#path <- "~/Lio/"
insertions <- read.csv(paste0(path,"10_nucleotide/all/ins-sep-all.csv"),row.names=1, stringsAsFactors = F)

# FIX HEADERS
#insertions$Header <- gsub("_\\d$","",insertions$Header)

# CASE: remove instances missing ancestor and tip 
insertions <- insertions[-c(which(insertions$anc == "")),]

# CASE: remove instances with gaps in the ancestor but NO INSERTION
insertions <- insertions[-c(which(grepl("-",insertions$anc) & insertions$indel=="")),]

# CASE: remove instances with insertion position 0
insertions <- insertions[-c(which(!is.na(insertions$pos) & insertions$pos==0)),]

# Restore all deletions found in tip sequences and adjust the POS values accordingly
insertions$tip <- unname(mapply(restoreTipDel,insertions$tip, insertions$anc, insertions$indel, insertions$pos))
# insertions$tip <- as.character(res[,1])
# insertions$pos <- as.numeric(as.character(res[,2]))

# CASE: restore all gaps from OTHER insertions
insertions$anc <- mapply(removeOtherGaps, insertions$anc, insertions$tip, insertions$indel, insertions$pos)

# SANITY CHECK: to make sure all seqs are equal
a <- nchar(insertions$tip) - nchar(insertions$indel)
b <- nchar(gsub("-","",insertions$anc))
sum(a!=b)==0              # should be TRUE
insertions <- insertions[-which(a!=b),]  # should be nrow = 0

# CASE: replace "R" nucleotides with the corresponding one found in the ancestor
cases <- which(grepl("[RYSWKMBDHVN]", insertions$tip))
for (idx in cases){
  toEdit <- insertions[idx, "tip"]
  pos <- gregexpr("[RYSWKMBDHVN]", toEdit)[[1]]
  toUse <- insertions[idx,"anc"]
  insertions[idx, "tip"] <- paste0(substr(toEdit, 1, pos-1),substr(toUse,pos,pos), substr(toEdit, pos+1,nchar(toEdit)))
}

setup(insertions$tip, insertions$anc, nchar(insertions$indel), insertions$pos, insertions$Date)
# add the a/b replicate label to the headers
#insertions$Header <- unname(mapply(patLabel, insertions$Header, insertions$Pat))
names(slip.list) <- insertions$Header

# # C.-.-.QT.10R.-.-_289_1_b


# needed for use in the CHANGESLIP function
idx <- which(unname(lapply(slip.list,sum))>0)

nt <- c("A", "C", "G", "T")




# RUN MCMC
startvalue <- c(0.001, 0.8, 0.001)
chain <- runMCMC(startvalue, 100000, slip.list)




