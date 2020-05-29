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
insertions <- read.csv(paste0(path,"10_nucleotide/all/ins-sep-all2.csv"),row.names=1, stringsAsFactors = F)

# FIX HEADERS
#insertions$Header <- gsub("_\\d$","",insertions$Header)

# CASE: remove instances missing ancestor and tip 
insertions <- insertions[-c(which(insertions$Anc == "")),]

# CASE: remove instances with gaps in the ancestor but NO INSERTION
insertions <- insertions[-c(which(grepl("-",insertions$Anc) & insertions$Seq=="")),]

# CASE: remove instances with insertion position 0
insertions <- insertions[-c(which(insertions$Pos==0)),]

# Restore all deletions found in tip sequences and adjust the POS values accordinaly
res <- as.data.frame(t(unname(mapply(restoreTipDel,insertions$Vseq, insertions$Anc, insertions$Seq, insertions$Pos))))
insertions$Vseq <- as.character(res[,1])
insertions$Pos <- as.numeric(as.character(res[,2]))

# CASE: restore all gaps from OTHER insertions
insertions$Anc <- mapply(removeOtherGaps, insertions$Anc, insertions$Vseq, insertions$Seq, insertions$Pos)

# SANITY CHECK: to make sure all seqs are equal
a <- nchar(insertions$Vseq) - nchar(insertions$Seq)
b <- nchar(gsub("-","",insertions$Anc))
sum(a!=b)==0              # should be TRUE
insertions[which(a!=b),]  # should be nrow = 0

# CASE: replace "R" nucleotides with the corresponding one found in the ancestor
cases <- which(grepl("[RYSWKMBDHVN]", insertions$Vseq))
for (idx in cases){
  toEdit <- insertions[idx, "Vseq"]
  pos <- gregexpr("[RYSWKMBDHVN]", toEdit)[[1]]
  toUse <- insertions[idx,"Anc"]
  insertions[idx, "Vseq"] <- paste0(substr(toEdit, 1, pos-1),substr(toUse,pos,pos), substr(toEdit, pos+1,nchar(toEdit)))
}

setup(insertions$Vseq, insertions$Anc, nchar(insertions$Seq), insertions$Pos, insertions$Date)
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




